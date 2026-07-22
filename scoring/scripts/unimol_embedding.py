"""Uni-Mol representation extraction through the supported unimol-tools API."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Sequence

import numpy as np
import torch

try:
    from ..asset_integrity import verify_asset
except ImportError:  # PYTHONPATH=scoring compatibility.
    from asset_integrity import verify_asset


def _stable_padding_length(n_tokens: int, *, max_tokens: int = 258) -> int:
    """Return a molecule-only padding length, independent of batch composition."""
    n_tokens = int(n_tokens)
    max_tokens = int(max_tokens)
    if n_tokens < 1:
        raise ValueError("n_tokens must be positive")
    if max_tokens < 1:
        raise ValueError("max_tokens must be positive")
    if n_tokens > max_tokens:
        return n_tokens
    for bucket in (32, 64, 128):
        if bucket < max_tokens and n_tokens <= bucket:
            return bucket
    return max_tokens


def _deterministic_crop_indices(n_atoms: int, max_atoms: int) -> np.ndarray:
    """Select stable endpoints when Uni-Mol must crop an oversized molecule."""
    n_atoms = int(n_atoms)
    max_atoms = int(max_atoms)
    if n_atoms < 1 or max_atoms < 1:
        raise ValueError("n_atoms and max_atoms must be positive")
    if n_atoms <= max_atoms:
        return np.arange(n_atoms, dtype=np.int64)
    return np.linspace(0, n_atoms - 1, max_atoms, dtype=np.int64)


def _build_unimol_features(smiles: list[str], params: dict) -> list[dict]:
    """Build features while making Uni-Mol's oversized-atom crop deterministic."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from unimol_tools.data.datahub import DataHub
    import unimol_tools.data.conformer as conformer

    max_atoms = int(params.get("max_atoms", 256))
    oversized = []
    for index, smi in enumerate(smiles):
        mol = Chem.MolFromSmiles(smi)
        if mol is not None and Chem.AddHs(mol).GetNumAtoms() > max_atoms:
            oversized.append(index)

    def load(indices: list[int], *, multi_process: bool) -> list[dict]:
        if not indices:
            return []
        local_params = dict(params)
        local_params["multi_process"] = bool(multi_process)
        datahub = DataHub(
            data=np.asarray([smiles[index] for index in indices]),
            task="repr",
            is_train=False,
            **local_params,
        )
        return datahub.data["unimol_input"]

    if not oversized:
        return load(list(range(len(smiles))), multi_process=bool(params.get("multi_process", False)))

    features: list[dict | None] = [None] * len(smiles)
    normal = [index for index in range(len(smiles)) if index not in set(oversized)]
    for index, feature in zip(normal, load(normal, multi_process=bool(params.get("multi_process", False)))):
        features[index] = feature

    original_choice = conformer.np.random.choice

    def deterministic_choice(a, size=None, replace=True, p=None):
        n_items = int(a) if isinstance(a, (int, np.integer)) else len(a)
        if not replace and isinstance(size, (int, np.integer)) and int(size) < n_items:
            return _deterministic_crop_indices(n_items, int(size))
        return original_choice(a, size=size, replace=replace, p=p)

    conformer.np.random.choice = deterministic_choice
    try:
        generator_params = dict(params)
        generator_params["multi_process"] = False
        generator = conformer.ConformerGen(**generator_params)
        for index in oversized:
            molecule = Chem.MolFromSmiles(smiles[index])
            if molecule is None:
                raise ValueError(f"invalid SMILES: {smiles[index]}")
            molecule = Chem.AddHs(molecule)
            AllChem.Compute2DCoords(molecule)
            atoms = np.asarray([atom.GetSymbol() for atom in molecule.GetAtoms()])
            coordinates = molecule.GetConformer().GetPositions().astype(np.float32)
            features[index] = conformer.coords2unimol(
                atoms,
                coordinates,
                generator.dictionary,
                generator.max_atoms,
                remove_hs=generator.remove_hs,
            )
    finally:
        conformer.np.random.choice = original_choice
    if any(feature is None for feature in features):
        raise RuntimeError("Uni-Mol feature construction lost an input row")
    return [feature for feature in features if feature is not None]


class UniMolEmbedding:
    """Extract deterministic 512-dimensional Uni-Mol CLS representations."""

    def __init__(
        self,
        weights_path: str,
        dict_path: str,
        config_path: str | None = None,
        device: str = "cuda",
        batch_size: int = 32,
        multi_process: bool = False,
    ):
        del config_path  # The public API owns the architecture matching these weights.
        weights = Path(weights_path).resolve()
        dictionary = Path(dict_path).resolve()
        if not weights.is_file():
            raise FileNotFoundError(f"Uni-Mol weights not found: {weights}")
        if not dictionary.is_file():
            raise FileNotFoundError(f"Uni-Mol dictionary not found: {dictionary}")
        if weights.parent != dictionary.parent:
            raise ValueError("Uni-Mol weights and dictionary must share one frozen asset directory")
        verify_asset(weights)
        verify_asset(dictionary)

        self.weight_dir = weights.parent
        if device.startswith("mps") and torch.backends.mps.is_available():
            self.device = "mps"
        else:
            self.device = "cuda" if device.startswith("cuda") and torch.cuda.is_available() else "cpu"
        self.batch_size = int(batch_size)
        self.multi_process = bool(multi_process)
        if self.batch_size < 1:
            raise ValueError("batch_size must be positive")

        self._configure_weight_dir(self.weight_dir)
        from unimol_tools import UniMolRepr

        self.repr_model = UniMolRepr(
            data_type="molecule",
            batch_size=self.batch_size,
            remove_hs=False,
            use_cuda=self.device == "cuda",
        )
        # unimol-tools otherwise spawns a process pool on POSIX for every call,
        # which is unsafe inside pytest/CLI workers and slower for our shards.
        self.repr_model.params["multi_process"] = self.multi_process

    @staticmethod
    def _configure_weight_dir(weight_dir: Path) -> None:
        """Point already-imported and future unimol-tools modules at frozen local assets."""
        path = str(weight_dir)
        os.environ["UNIMOL_WEIGHT_DIR"] = path
        import unimol_tools.data.conformer as conformer
        import unimol_tools.models.unimol as unimol_model
        import unimol_tools.weights as weights_package
        import unimol_tools.weights.weighthub as weighthub

        conformer.WEIGHT_DIR = path
        unimol_model.WEIGHT_DIR = path
        weights_package.WEIGHT_DIR = path
        weighthub.WEIGHT_DIR = path

    def extract_one(self, smiles: str) -> np.ndarray:
        return self.extract([smiles], batch_size=1, verbose=False)[0]

    def extract(
        self,
        smiles_list: Sequence[str],
        batch_size: int | None = None,
        verbose: bool = False,
    ) -> np.ndarray:
        del verbose
        smiles = list(smiles_list)
        if not smiles:
            return np.empty((0, 512), dtype=np.float32)
        if not all(isinstance(item, str) and item.strip() for item in smiles):
            raise ValueError("all Uni-Mol inputs must be non-empty SMILES strings")
        requested_batch_size = int(batch_size or self.batch_size)
        if requested_batch_size < 1:
            raise ValueError("batch_size must be positive")
        self.repr_model.params["batch_size"] = requested_batch_size
        if self.device == "mps":
            output = self._extract_mps(smiles, requested_batch_size)
        else:
            output = self.repr_model.get_repr(smiles)
        embeddings = np.asarray(output["cls_repr"], dtype=np.float32)
        if embeddings.shape != (len(smiles), 512):
            raise RuntimeError(
                f"unexpected Uni-Mol representation shape {embeddings.shape}; "
                f"expected {(len(smiles), 512)}"
            )
        if not np.isfinite(embeddings).all():
            raise RuntimeError("Uni-Mol returned non-finite representations")
        return embeddings

    def _extract_mps(self, smiles: list[str], batch_size: int) -> dict:
        """Run the public Uni-Mol model with its official collator on Apple MPS."""
        from unimol_tools.predictor import MolDataset
        from unimol_tools.tasks.trainer import NNDataLoader

        params = dict(self.repr_model.params)
        params["batch_size"] = batch_size
        params["use_cuda"] = False
        params["multi_process"] = self.multi_process
        features = _build_unimol_features(smiles, params)
        # UniMol's stock collator pads to the largest molecule *in the current
        # call*.  On MPS this changes the numerical path through the attention
        # kernels, so the same molecule can receive a different CLS vector in
        # a different batch.  Use a frozen maximum length to make the scorer a
        # function of the molecule and model assets only.
        configured_max_atoms = int(params.get("max_atoms", 256))
        max_tokens = configured_max_atoms + 2
        padding_idx = int(self.repr_model.model.padding_idx)
        bucketed_features: dict[int, list[tuple[int, dict]]] = {}
        for original_index, feature in enumerate(features):
            n_tokens = len(feature["src_tokens"])
            fixed_length = _stable_padding_length(n_tokens, max_tokens=max_tokens)
            padded = dict(feature)
            tokens = np.full(fixed_length, padding_idx, dtype=np.int64)
            tokens[:n_tokens] = np.asarray(feature["src_tokens"], dtype=np.int64)
            coords = np.zeros((fixed_length, 3), dtype=np.float32)
            coords[:n_tokens] = np.asarray(feature["src_coord"], dtype=np.float32)
            distance = np.zeros((fixed_length, fixed_length), dtype=np.float32)
            distance[:n_tokens, :n_tokens] = np.asarray(feature["src_distance"], dtype=np.float32)
            edge_type = np.full((fixed_length, fixed_length), padding_idx, dtype=np.int64)
            edge_type[:n_tokens, :n_tokens] = np.asarray(feature["src_edge_type"], dtype=np.int64)
            padded.update({
                "src_tokens": tokens,
                "src_coord": coords,
                "src_distance": distance,
                "src_edge_type": edge_type,
            })
            bucketed_features.setdefault(fixed_length, []).append((original_index, padded))

        model = self.repr_model.model.to(torch.device("mps")).eval()
        restored = [None] * len(features)
        for fixed_length in sorted(bucketed_features):
            group = bucketed_features[fixed_length]
            dataset = MolDataset([feature for _, feature in group])
            loader = NNDataLoader(
                feature_name=None,
                dataset=dataset,
                batch_size=batch_size,
                shuffle=False,
                collate_fn=model.batch_collate_fn,
                distributed=False,
                valid_mode=True,
            )
            group_offset = 0
            for batch in loader:
                net_input, _ = batch
                net_input = {key: value.to(torch.device("mps")) for key, value in net_input.items()}
                captured = []

                def capture_cls(_module, inputs):
                    captured.append(inputs[0])

                hook = model.classification_head.register_forward_pre_hook(capture_cls)
                with torch.no_grad():
                    model(**net_input, return_repr=False)
                hook.remove()
                if len(captured) != 1:
                    raise RuntimeError("Uni-Mol classification-head hook did not capture CLS representation")
                for item in captured[0]:
                    original_index = group[group_offset][0]
                    restored[original_index] = item.detach().cpu().numpy()
                    group_offset += 1
        return {"cls_repr": restored}


if __name__ == "__main__":
    base = Path(__file__).resolve().parent.parent
    assets = (
        base
        / "models"
        / "unimol"
        / "models--dptech--Uni-Mol-Models"
        / "snapshots"
        / "9f19c45c718192888a1c8a1c905f69f0755ea502"
    )
    extractor = UniMolEmbedding(
        weights_path=str(assets / "mol_pre_all_h_220816.pt"),
        dict_path=str(assets / "mol.dict.txt"),
        device="cpu",
    )
    values = extractor.extract(["c1ccccc1", "CCO"], batch_size=2)
    print(values.shape, float(values.std()))

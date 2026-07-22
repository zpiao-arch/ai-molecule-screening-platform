from __future__ import annotations

from dataclasses import dataclass
import os
from pathlib import Path


@dataclass(frozen=True)
class BenchmarkConfig:
    root: Path
    n_targets: int = 1000
    pool_size: int = 10000
    seed: int = 42

    def __post_init__(self):
        object.__setattr__(self, "root", Path(self.root).resolve())

    @property
    def scoring_dir(self) -> Path:
        return Path(os.environ.get("FOUR_LEVEL_SCORING_DIR", self.root / "scoring")).resolve()

    @property
    def bindingdb_examples(self) -> Path:
        return Path(os.environ.get(
            "FOUR_LEVEL_BINDINGDB_EXAMPLES",
            self.root / "data" / "bindingdb_202606_target_match_examples.csv",
        )).resolve()

    @property
    def chembl_examples(self) -> Path:
        return Path(os.environ.get(
            "FOUR_LEVEL_CHEMBL_EXAMPLES",
            self.root / "data" / "chembl_37_target_match_examples.csv",
        )).resolve()

    @property
    def chembl_uniprot_mapping(self) -> Path:
        return Path(os.environ.get(
            "FOUR_LEVEL_CHEMBL_UNIPROT_MAPPING",
            self.root / "data" / "chembl_uniprot_mapping.txt",
        )).resolve()

    @property
    def chembl_library(self) -> Path:
        return Path(os.environ.get(
            "FOUR_LEVEL_CHEMBL_LIBRARY",
            self.root / "data" / "chembl_37_chemreps.txt.gz",
        )).resolve()

    @property
    def frozen_library_smiles(self) -> Path:
        return Path(os.environ.get(
            "FOUR_LEVEL_LIBRARY_SMILES",
            self.root / "data" / "lib_smiles.txt",
        )).resolve()

    @property
    def frozen_library_features(self) -> Path:
        return Path(os.environ.get(
            "FOUR_LEVEL_LIBRARY_FEATURES",
            self.root / "data" / "lib_feats.npy",
        )).resolve()

    @property
    def l2_model_path(self) -> Path:
        compatible = self.scoring_dir / "models" / "bindingdb_l2" / "l2_model_sklearn_1_7_2.joblib"
        return compatible if compatible.is_file() else self.scoring_dir / "models" / "bindingdb_l2" / "l2_model.joblib"

"""
Uni-Mol Zero-Shot 分子评分
用预训练嵌入 + 余弦相似度，无需训练

用法:
    scorer = UniMolScorer()
    score = scorer.score("CC(=O)Oc1ccccc1C(=O)O")
    # → {"unimol_score": 0.78, "top_refs": [...]}
"""

import numpy as np
from pathlib import Path
from sklearn.metrics.pairwise import cosine_similarity

try:
    from .unimol_embedding import UniMolEmbedding
except ImportError:  # PYTHONPATH=scoring compatibility.
    from scripts.unimol_embedding import UniMolEmbedding

try:
    from ..asset_integrity import verify_asset
except ImportError:  # PYTHONPATH=scoring compatibility.
    from asset_integrity import verify_asset


# 从文件加载FDA批准药物参考集
def _load_drugbank_refs():
    """加载FDA批准药物参考集"""
    ref_path = Path(__file__).parent.parent / "models" / "drugbank_ref.txt"
    drugs = []
    if ref_path.exists():
        for line in ref_path.read_text(encoding="utf-8").split("\n"):
            line = line.strip()
            if line and not line.startswith("#"):
                smi = line.split(",")[0].strip()
                if smi and len(smi) > 2:
                    drugs.append((f"FDA_{len(drugs)+1}", smi))
        if len(drugs) > 20:
            return drugs
    return _default_refs()

def _default_refs():
    """内置回退参考集"""
    return [
        ("Aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
        ("Ibuprofen", "CC(C)Cc1ccc(C(C)C(=O)O)cc1"),
        ("Acetaminophen", "CC(=O)Nc1ccc(O)cc1"),
        ("Caffeine", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"),
        ("Diazepam", "CN1C(=O)CN=C(c2ccccc2Cl)c2ccccc21"),
        ("Naproxen", "COc1ccc2cc(C(C)C(=O)O)ccc2c1"),
        ("Oseltamivir", "CCC(CC)OC1C=C(CC(N)=O)C(O)CC1NC(C)=O"),
        ("Metformin", "CN(C)C(=N)NC(N)=N"),
        ("Omeprazole", "COc1ccnc(CS(=O)c2nc3ccccc3[nH]2)c1C"),
        ("Warfarin", "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O"),
        ("Darunavir", "CC(C)CN(CC(C)C)S(=O)(=O)C1=CC=C(NC(=O)OC2CCOC2)C=C1"),
        ("Ritonavir", "CC(C)C1=NC(=S)SC1=C(C#N)C(=O)NC(CC1=CC=CC=C1)CC1=CC=CC=C1"),
        ("Doxorubicin", "COc1cccc2C(=O)c3c(O)c4CC(OC5CC(N)C(O)C(C)O5)(C(=O)CO)c4c(O)c3C(=O)c12"),
        ("Penicillin_G", "CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O"),
        ("Ciprofloxacin", "OC(=O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O"),
    ]


# 已知有毒分子 (阴性参考)
TOXIC_REFS = [
    ("Aflatoxin_B1", "COc1cc2c(c3c1O[C@@H]1C=CO[C@@H]1[C@@H]3O)C(=O)CC2"),
    ("Thalidomide", "O=C1NC(=O)C2=C1CCN(C2=O)C1=CC=CC=C1"),
    ("Paraquat", "C[N+]1=CC=C(C=C1)c2cc[n+](C)cc2"),
]


class UniMolScorer:
    """Uni-Mol零-shot评分器"""

    def __init__(self, device: str = "cuda", multi_process: bool = False):
        BASE = Path(__file__).parent.parent  # scripts/ → 评分/

        # 尝试加载缓存
        cache_path = BASE / "models/ref_embeddings.npz"
        smiles_cache = BASE / "models/ref_smiles.pkl"

        if cache_path.exists() and smiles_cache.exists():
            # 从缓存加载 (秒级)
            verify_asset(cache_path)
            verify_asset(smiles_cache)
            data = np.load(cache_path)
            self.ref_embeddings = data["embeddings"]
            import pickle
            with open(smiles_cache, "rb") as f:
                cached_smiles = pickle.load(f)

            # 构建参考集元数据 (全部作为FDA药物)
            self.ref_smiles = {}
            self.ref_labels = {}
            self.ref_names = []
            n_fda = len(cached_smiles) - 3  # 排除最后3个TOXIC
            if n_fda < 100:
                n_fda = len(cached_smiles)  # 兼容旧缓存
            for i, smi in enumerate(cached_smiles[:n_fda]):
                name = f"FDA_{i+1}"
                self.ref_names.append(name)
                self.ref_smiles[name] = smi
                self.ref_labels[name] = "positive"
            self.ref_embeddings = self.ref_embeddings[:n_fda]  # 裁剪TOXIC

            print(f"从缓存加载参考集: {len(self.ref_names)} 个FDA药物")

        else:
            # 实时计算 (首次)
            WEIGHTS = BASE / "models/unimol/models--dptech--Uni-Mol-Models/snapshots/9f19c45c718192888a1c8a1c905f69f0755ea502"
            self.extractor = UniMolEmbedding(
                weights_path=str(WEIGHTS / "mol_pre_all_h_220816.pt"),
                dict_path=str(WEIGHTS / "mol.dict.txt"),
                device=device,
            )

            fda_drugs = _load_drugbank_refs()
            self.ref_smiles = {}
            self.ref_labels = {}
            for name, smi in fda_drugs:
                self.ref_smiles[name] = smi
                self.ref_labels[name] = "positive"

            print(f"预计算 {len(self.ref_smiles)} 个FDA药物嵌入...")
            ref_names = list(self.ref_smiles.keys())
            ref_smiles_list = [self.ref_smiles[n] for n in ref_names]
            self.ref_embeddings = self.extractor.extract(ref_smiles_list, verbose=False)
            self.ref_names = ref_names

        # 归一化
        self.ref_norms = np.linalg.norm(self.ref_embeddings, axis=1, keepdims=True)
        self.ref_normalized = self.ref_embeddings / (self.ref_norms + 1e-8)

        self._cache = {}
        self.device = device
        self.multi_process = bool(multi_process)
        self.extractor = None  # 候选分子实时计算时懒加载

        print(f"参考集就绪: {len(self.ref_names)} 个分子")

    def _get_extractor(self):
        if self.extractor is None:
            BASE = Path(__file__).parent.parent
            WEIGHTS = BASE / "models/unimol/models--dptech--Uni-Mol-Models/snapshots/9f19c45c718192888a1c8a1c905f69f0755ea502"
            self.extractor = UniMolEmbedding(
                weights_path=str(WEIGHTS / "mol_pre_all_h_220816.pt"),
                dict_path=str(WEIGHTS / "mol.dict.txt"),
                device=self.device,
                multi_process=self.multi_process,
            )
        return self.extractor

    def score(self, smiles: str) -> dict:
        """对单个分子打分"""
        try:
            emb = self._get_embedding(smiles)
        except Exception as exc:
            return {"unimol_score": 0.0, "pos_similarity": 0.0, "neg_similarity": 0.0,
                    "top_refs": [], "status": f"error:{type(exc).__name__}",
                    "error": str(exc)}

        return self._score_embedding(emb)

    def _score_embedding(self, emb: np.ndarray) -> dict:
        """Convert one verified representation into the frozen L4 similarity score."""

        emb_norm = emb / (np.linalg.norm(emb) + 1e-8)

        # 计算与FDA药物的余弦相似度 (只看positive)
        similarities = cosine_similarity([emb_norm], self.ref_normalized)[0]

        # Top-10最相似FDA药物
        K = 10
        pos_similarities = [(similarities[i], i) for i in range(len(similarities))
                            if self.ref_labels[self.ref_names[i]] == "positive"]
        pos_similarities.sort(key=lambda x: x[0], reverse=True)
        top_K = pos_similarities[:K]

        # 平均相似度 → L4分数
        if top_K:
            avg_sim = sum(s[0] for s in top_K) / len(top_K)
        else:
            avg_sim = 0.0

        unimol_score = round(avg_sim, 4)
        top_refs = [(self.ref_names[i], round(float(similarities[i]), 4)) for _, i in top_K[:5]]

        return {
            "unimol_score": unimol_score,
            "pos_similarity": round(float(avg_sim), 4),
            "neg_similarity": 0.0,
            "top_refs": top_refs,
            "status": "ok",
        }

    def _get_embedding(self, smiles: str) -> np.ndarray:
        if smiles in self._cache:
            return self._cache[smiles]
        ext = self._get_extractor()
        emb = ext.extract_one(smiles)
        self._cache[smiles] = emb
        return emb

    def score_many(self, smiles_list: list, batch_size: int = 32) -> list:
        """Score SMILES in real Uni-Mol batches while retaining per-row failures."""
        smiles = list(smiles_list)
        missing = list(dict.fromkeys(smi for smi in smiles if smi not in self._cache))
        if missing:
            try:
                embeddings = self._get_extractor().extract(
                    missing, batch_size=batch_size, verbose=False
                )
                self._cache.update(zip(missing, embeddings))
            except Exception:
                # Isolate malformed molecules without turning their failures into zeros.
                for smi in missing:
                    try:
                        self._cache[smi] = self._get_extractor().extract_one(smi)
                    except Exception:
                        continue

        rows = []
        for smi in smiles:
            emb = self._cache.get(smi)
            if emb is None:
                rows.append({
                    "unimol_score": None,
                    "pos_similarity": None,
                    "neg_similarity": None,
                    "top_refs": [],
                    "status": "error:embedding_failed",
                    "error": "embedding_failed",
                })
            else:
                rows.append(self._score_embedding(emb))
        return rows

    def score_batch(self, molecules: list) -> list:
        """批量打分: [(id, smiles), ...]"""
        scores = self.score_many([smi for _, smi in molecules])
        results = []
        for (mid, smi), r in zip(molecules, scores):
            r = dict(r)
            r["id"] = mid
            r["smiles"] = smi
            results.append(r)
        return results


if __name__ == "__main__":
    scorer = UniMolScorer(device="cpu")  # CPU模式更稳定

    # 测试几个分子
    tests = [
        ("test_oseltamivir", "CCC(CC)OC1C=C(CC(N)=O)C(O)CC1NC(C)=O"),
        ("test_aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
        ("test_benzene", "c1ccccc1"),
        ("test_thalidomide", "O=C1NC(=O)C2=C1CCN(C2=O)C1=CC=CC=C1"),
    ]

    for mid, smi in tests:
        r = scorer.score(smi)
        print(f"\n{mid}:")
        print(f"  Uni-Mol分: {r['unimol_score']:.4f}")
        print(f"  阳性相似: {r['pos_similarity']:.4f}  阴性相似: {r['neg_similarity']:.4f}")
        print(f"  Top参考: {r['top_refs'][:3]}")

"""
预计算FDA药物参考集嵌入缓存
每100个保存一次，支持断点续传
"""

import pickle
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))

try:
    from .unimol_embedding import UniMolEmbedding
    from .unimol_scorer import _load_drugbank_refs, TOXIC_REFS
except ImportError:  # Direct script/PYTHONPATH=scoring compatibility.
    from scripts.unimol_embedding import UniMolEmbedding
    from scripts.unimol_scorer import _load_drugbank_refs, TOXIC_REFS

BASE = Path(__file__).parent.parent
WEIGHTS = BASE / "models/unimol/models--dptech--Uni-Mol-Models/snapshots/9f19c45c718192888a1c8a1c905f69f0755ea502"
CACHE = BASE / "models/ref_embeddings.npz"
REF_SMILES = BASE / "models/ref_smiles.pkl"


def main():
    print("加载参考药物列表...")
    fda_drugs = _load_drugbank_refs()
    refs = [(name, smi, "positive") for name, smi in fda_drugs]
    refs += [(name, smi, "negative") for name, smi in TOXIC_REFS]

    all_smiles = [smi for _, smi, _ in refs]
    total = len(all_smiles)
    print(f"共计 {total} 个参考分子 (正:{len(fda_drugs)}, 负:{len(TOXIC_REFS)})")

    # 检查断点
    if CACHE.exists() and REF_SMILES.exists():
        existing = np.load(CACHE)["embeddings"]
        done = existing.shape[0]
        print(f"已有缓存: {done}/{total}")
        if done >= total:
            print("缓存已完整，无需重新计算")
            return
    else:
        done = 0
        existing = np.empty((0, 512), dtype=np.float32)

    print("加载 Uni-Mol 模型...")
    extractor = UniMolEmbedding(
        weights_path=str(WEIGHTS / "mol_pre_all_h_220816.pt"),
        dict_path=str(WEIGHTS / "mol.dict.txt"),
        device="cpu",
    )

    batch = 100
    for start in range(done, total, batch):
        end = min(start + batch, total)
        print(f"[{start+1}-{end}/{total}] 计算中...", end=" ", flush=True)
        batch_smiles = all_smiles[start:end]
        embs = extractor.extract(batch_smiles, verbose=False)

        existing = np.vstack([existing, embs])
        np.savez_compressed(CACHE, embeddings=existing)
        with open(REF_SMILES, "wb") as f:
            pickle.dump(all_smiles[:len(existing)], f)

        print(f"已保存 ({existing.shape[0]}/{total})")

    print(f"\n完成! 缓存: {CACHE} ({CACHE.stat().st_size/1024:.0f} KB)")


if __name__ == "__main__":
    main()

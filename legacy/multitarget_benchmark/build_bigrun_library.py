# -*- coding: utf-8 -*-
"""
大通量闭环 Stage A — 构建 100,000 分子库 + 标注映射
================================================================================
库 = 49,947 已缓存标注分子(含全部已知活性) + 从 chembl_37_chemreps 采样的干扰分子
     补齐到 100,000。干扰分子无需标签(是"草垛"), 只让检索更难更真实。
产出 (checkpoint 到 <external-library-cache>/):
  - lib_feats.npy      (100000, 520) float32 分子特征
  - lib_smiles.txt     100000 行 SMILES (与 feats 行对齐)
  - label_map.json     {chembl_target_id: {"pos":[smiles...], "neg":[...]}}  (来自 aligned 标注)
"""
import os, sys, csv, gzip, json, pickle, time, random
import numpy as np
from pathlib import Path

REPO = Path("<validated-workspace>")
SC = REPO / "评分_work_package/评分"
OUT = Path("<external-library-cache>"); OUT.mkdir(exist_ok=True)
sys.path.insert(0, str(SC))
csv.field_size_limit(10**7)

from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")
from l2_bindingdb import BindingDBFeature
from target_seq_embedding import extract_chembl_id

TARGET_N = 100_000
SEED = 42
t0 = time.time()
def log(*a): print(f"[{time.time()-t0:6.0f}s]", *a, flush=True)

feat = BindingDBFeature()

# ── 1. 载入已缓存标注分子特征 (49,947) ──
cache = pickle.load(open(SC / "models/bindingdb_l2_seq/molfeat_cache.pkl", "rb"))
labeled_smiles = list(cache.keys())
labeled_feats = np.stack([np.asarray(cache[s], dtype=np.float32) for s in labeled_smiles])
log(f"标注分子(缓存): {len(labeled_smiles)}  feat {labeled_feats.shape}")
have = set(labeled_smiles)

# ── 2. 构建标注 target->actives 映射 (ChEMBL 命名空间) ──
CH = REPO / "data_lake/chembl/latest/aligned_model_input/chembl_37_target_match_examples.csv"
BD = REPO / "data_lake/bindingdb/aligned_model_input/bindingdb_202606_target_match_examples.csv"
label_map = {}
def add(tid, smi, lab):
    if not tid or not smi: return
    d = label_map.setdefault(tid, {"pos": set(), "neg": set()})
    (d["pos"] if lab == 1 else d["neg"]).add(smi)

with open(CH) as f:
    for r in csv.DictReader(f):
        tid = r.get("target_chembl_id", "").strip()
        smi = r.get("canonical_smiles", "").strip()
        try: lab = int(float(r.get("label", "0")))
        except: continue
        add(tid, smi, lab)
with open(BD) as f:
    for r in csv.DictReader(f):
        tid = extract_chembl_id(r.get("target_text", "") or "")
        smi = r.get("canonical_smiles", "").strip()
        try: lab = int(float(r.get("label", "0")))
        except: continue
        add(tid, smi, lab)
# 序列化 (set->list)
label_map_ser = {t: {"pos": sorted(d["pos"]), "neg": sorted(d["neg"])} for t, d in label_map.items()}
json.dump(label_map_ser, open(OUT / "label_map.json", "w"))
n_tgt_lab = sum(1 for d in label_map.values() if len(d["pos"]) >= 5)
log(f"标注靶点: {len(label_map)} (其中 >=5 活性: {n_tgt_lab})")

# ── 3. 采样干扰分子补齐到 100k ──
need = TARGET_N - len(labeled_smiles)
log(f"需干扰分子 {need} 个")
rng = random.Random(SEED)
CR = REPO / "data_lake/chembl/latest/chembl_37_chemreps.txt.gz"
# 蓄水池: 先读全部 chembl id->smiles (仅 smiles 列), 跳过已在库中的
distract = []
seen = 0
with gzip.open(CR, "rt") as f:
    header = f.readline()
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 2: continue
        smi = parts[1].strip()
        seen += 1
        if not smi or smi in have: continue
        # 蓄水池抽样
        if len(distract) < need * 3:   # 先多收 3x 备用(部分特征会失败)
            distract.append(smi)
        elif rng.random() < (need * 3) / seen:
            distract[rng.randrange(len(distract))] = smi
log(f"chemreps 扫描 {seen} 行, 候选干扰 {len(distract)}")

# ── 4. 计算干扰分子特征 (直到补足 need) ──
dist_smiles, dist_feats = [], []
rng.shuffle(distract)
for i, smi in enumerate(distract):
    if len(dist_feats) >= need: break
    mf = feat.mol_features(smi)
    if mf is None: continue
    dist_smiles.append(smi); dist_feats.append(np.asarray(mf, dtype=np.float32))
    if len(dist_feats) % 5000 == 0:
        log(f"干扰特征 {len(dist_feats)}/{need}")
log(f"干扰分子最终 {len(dist_feats)}")

# ── 5. 合并并保存 ──
all_smiles = labeled_smiles + dist_smiles
all_feats = np.vstack([labeled_feats, np.stack(dist_feats)])
np.save(OUT / "lib_feats.npy", all_feats)
with open(OUT / "lib_smiles.txt", "w") as f:
    f.write("\n".join(all_smiles))
# 记录标注分子在库中的行号 (前 len(labeled) 行)
meta = {"n_total": len(all_smiles), "n_labeled": len(labeled_smiles),
        "n_distract": len(dist_smiles), "feat_dim": int(all_feats.shape[1])}
json.dump(meta, open(OUT / "lib_meta.json", "w"), indent=1)
log(f"库完成: {all_feats.shape} 保存到 {OUT}")
print("STAGE_A_DONE", json.dumps(meta))

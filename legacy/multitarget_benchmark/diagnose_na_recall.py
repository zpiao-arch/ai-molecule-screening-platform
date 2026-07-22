#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""NA 活性分子 L2 召回失败 — 根因归因诊断(可复现).
把 18 个 NA 活性在 100k 库里的真实 L2 排名/分位挖出来, 并量化 L2 分数塌缩,
定位"罪魁": 校准崩坏 / 靶点条件化失效 / 库本身全为药物样 / 标签-库 SMILES 对不齐.
"""
import os, sys, json, time
import numpy as np
sys.path.insert(0, "评分_work_package/评分")
from l2_bindingdb import BindingDBFeature, Layer2BindingDB
from target_resolver import TargetResolver
from rdkit import Chem

BASE = "<validated-workspace>"
LIB_FEATS = "<external-library-cache>/lib_feats.npy"
LIB_SMI = "<external-library-cache>/lib_smiles.txt"
LABEL = "<external-library-cache>/label_map.json"
TARGET = "CHEMBL2051"

def canon(s):
    try:
        return Chem.MolToSmiles(Chem.MolFromSmiles(s))
    except Exception:
        return None

t0 = time.time()
feats = np.load(LIB_FEATS)
smiles = [l.strip() for l in open(LIB_SMI) if l.strip()]
lib_canon = {}
for i, s in enumerate(smiles):
    c = canon(s)
    if c: lib_canon[c] = i
N = len(smiles)
print(f"[load] 库 {N} 分子, 特征矩阵 {feats.shape}, 用时 {time.time()-t0:.1f}s")

# ── L2 打分 ──
f = BindingDBFeature(); r = TargetResolver()
tt, method = r.resolve(name="neuraminidase", chembl_id=TARGET)
th = f.target_features(tt)
X = np.hstack([feats, np.tile(th, (feats.shape[0], 1))]).astype(np.float32)
l2 = Layer2BindingDB(); l2._ensure_model()
proba = l2.predict_proba(X)
if proba.ndim > 1: proba = proba[:, 1]
proba = proba.astype(float)
med = float(np.median(proba))
print(f"[L2] 打分完成, 背景中位={med:.3f}, 用时 {time.time()-t0:.1f}s")

# ── 标签-库对齐 ──
lm = json.load(open(LABEL))[TARGET]
pos = lm["pos"]; neg = lm["neg"]
pos_exact = [s for s in pos if s in set(smiles)]
pos_canon = [s for s in pos if canon(s) in lib_canon]
neg_exact = [s for s in neg if s in set(smiles)]
neg_canon = [s for s in neg if canon(s) in lib_canon]
print(f"[align] 活性 pos: 精确匹配 {len(pos_exact)}/{len(pos)}, 规范化匹配 {len(pos_canon)}/{len(pos)}")
print(f"[align] 阴性 neg: 精确匹配 {len(neg_exact)}/{len(neg)}, 规范化匹配 {len(neg_canon)}/{len(neg)}")

# ── 18 活性真实排名 ──
order = np.argsort(-proba)            # 降序
rank_of = {int(j): int(k) for k, j in enumerate(order)}  # idx -> 名次(0-based)
rows = []
for s in pos:
    c = canon(s)
    if c and c in lib_canon:
        idx = lib_canon[c]
        rank = rank_of[idx] + 1
        rows.append({"smiles": s, "idx": idx, "l2": round(float(proba[idx]), 4),
                     "rank": rank, "pct": round(100 * rank / N, 2)})
rows.sort(key=lambda d: d["rank"])
print(f"\n[actives] {len(rows)} 个活性在库中的 L2 排名:")
for d in rows:
    print(f"  名次 {d['rank']:>6}/{N} ({d['pct']:>5.2f}%)  L2={d['l2']:.3f}  {d['smiles'][:40]}")
ranks = [d["rank"] for d in rows]
if ranks:
    print(f"[actives] 中位排名={int(np.median(ranks))}, 最靠前={min(ranks)}, 最靠后={max(ranks)}")
    # 假设 top-N 截断为 250, 落在窗口内的活性数
    for topN in (250, 500, 1000, 2000, 5000):
        print(f"  若 L2 top-{topN} 接对接: 落入窗口的活性数 = {sum(1 for x in ranks if x <= topN)}/{len(ranks)}")

# ── 分数分布(塌缩量化) ──
counts, edges = np.histogram(proba, bins=20, range=(0, 1))
gt09 = int((proba > 0.9).sum()); gt095 = int((proba > 0.95).sum())
gt_med = int((proba > med).sum())
print(f"\n[dist] 分数>背景中位({med:.3f}) 的分子 = {gt_med} ({100*gt_med/N:.1f}%)")
print(f"[dist] 分数>0.90 的分子 = {gt09} ({100*gt09/N:.1f}%)")
print(f"[dist] 分数>0.95 的分子 = {gt095} ({100*gt095/N:.1f}%)")

# ── 阴性(诱饵)对照: 18活性 vs 随机阴性 的 L2 分布 ──
neg_idx = [lib_canon[canon(s)] for s in neg_canon if canon(s) in lib_canon]
if neg_idx:
    print(f"\n[contrast] 活性 L2 中位={np.median([d['l2'] for d in rows]):.3f}  vs  阴性 L2 中位={np.median(proba[neg_idx]):.3f}")
    print(f"[contrast] 活性 L2 均值={np.mean([d['l2'] for d in rows]):.3f}  vs  阴性 L2 均值={proba[neg_idx].mean():.3f}")

# ── 落盘 ──
out = {
    "target": TARGET, "library_size": N,
    "l2_median": round(med, 4),
    "frac_gt_median": round(gt_med / N, 4),
    "frac_gt_0_90": round(gt09 / N, 4),
    "frac_gt_0_95": round(gt095 / N, 4),
    "n_pos_in_lib": len(rows), "n_pos_label": len(pos),
    "pos_exact": len(pos_exact), "pos_canon": len(pos_canon),
    "n_neg_in_lib": len(neg_idx),
    "active_ranks": [d["rank"] for d in rows],
    "active_rank_median": int(np.median(ranks)) if ranks else None,
    "active_l2_median": round(float(np.median([d["l2"] for d in rows])), 4) if rows else None,
    "neg_l2_median": round(float(np.median(proba[neg_idx])), 4) if neg_idx else None,
    "topN_window_recall": {str(k): sum(1 for x in ranks if x <= k) for k in (250, 500, 1000, 2000, 5000)},
    "hist_edges": [round(float(e), 3) for e in edges],
    "hist_counts": [int(c) for c in counts],
}
json.dump(out, open("<external-library-cache>/diagnose_na_recall.json", "w"),
          indent=2, ensure_ascii=False)
print(f"\n[done] 用时 {time.time()-t0:.1f}s, 结果 -> <external-library-cache>/diagnose_na_recall.json")

# -*- coding: utf-8 -*-
"""
高通量测试闭环 (序列嵌入 L2 版) — 对比升级
============================================
复用 hts_closed_loop.py 已落盘的 L1/L3/MOLF 缓存, 仅把 L2 从
"靶点哈希 Layer2BindingDB(256维)" 替换为 "序列嵌入 Layer2BindingDBSeq(128维CNN)"。
对每个靶点用 predict_matrix(molfeat, chembl_id) 全库向量化打分, 同口径复算:
  recall@1/5/10%, EF@5%, BEDROC, 全库AUC, 级联剪枝。
末尾与旧哈希结果逐靶点对比, 重点看 "28 个反向靶点(AUC<0.5)" 是否修复。
"""
import json, csv, time, sys, pickle, math
from pathlib import Path
import numpy as np
from sklearn.metrics import roc_auc_score

REPO = Path("<validated-workspace>")
SC = REPO / "评分_work_package/评分"
MTB = REPO / "scientific_validation/multitarget_benchmark"
sys.path.insert(0, str(SC))
from l2_bindingdb import BindingDBFeature, Layer2BindingDBSeq

t0 = time.time()
def log(*a): print(f"[{time.time()-t0:6.0f}s]", *a, flush=True)

# ── 载入库 + 面板 ──
lib_ids, lib_smiles = [], []
with open(MTB / "candidates_10k.csv") as f:
    for r in csv.DictReader(f):
        lib_ids.append(r["id"]); lib_smiles.append(r["smiles"].strip())
panel = json.load(open(MTB / "target_panel.json"))
N = len(lib_smiles)
log(f"候选库 {N} 药; 面板 {len(panel)} 靶点")

# ── 载入 L1/L3/MOLF 缓存 ──
CACHE = MTB / "hts_l1l3_cache.pkl"
_loaded = pickle.load(open(CACHE, "rb"))
if len(_loaded) == 3:
    L1c, L3c, MOLF = _loaded
else:
    L1c, L3c = _loaded; MOLF = {}
log(f"载入 L1/L3/MOLF 缓存 {len(L1c)} 分子")

feat = BindingDBFeature()
L1_arr = np.array([L1c[s] for s in lib_smiles], dtype=np.float32)
L3_arr = np.array([L3c[s] for s in lib_smiles], dtype=np.float32)
_z = np.zeros(feat.DIM - feat.N_TARGET, dtype=np.float32)
molfeat = np.array([MOLF.get(s, _z) for s in lib_smiles], dtype=np.float32)  # (N,520)

from rdkit import Chem
def canon(s):
    try: return Chem.MolToSmiles(Chem.MolFromSmiles(s))
    except Exception: return None
lib_canon = {}
for i, s in enumerate(lib_smiles):
    c = canon(s)
    if c: lib_canon[c] = i
log(f"库内 canon 索引 {len(lib_canon)}")

# ── 指标 ──
def bedroc(pos_scores, neg_scores, alpha=20.0):
    all_s = np.concatenate([pos_scores, neg_scores])
    n = len(all_s); npos = len(pos_scores)
    if n == 0 or npos == 0: return 0.0
    order = np.argsort(-all_s)
    pos_set = set(range(npos)); ra = 0.0
    for rank, idx in enumerate(order):
        if idx in pos_set: ra += math.exp(-alpha * (rank + 1) / n)
    ra_max = (1 - math.exp(-alpha * npos / n)) / (1 - math.exp(-alpha))
    return min(1.0, ra / ra_max) if ra_max > 0 else 0.0

def recall_at(pos_idx, order, fracs):
    out = {}; npos = len(pos_idx)
    if npos == 0: return {f"recall@{int(f*100)}%": 0.0 for f in fracs}
    nlib = len(order)
    for f in fracs:
        K = max(1, int(round(nlib * f)))
        top = set(order[:K].tolist())
        out[f"recall@{int(f*100)}%"] = sum(1 for i in pos_idx if i in top) / npos
    return out

def ef_at(pos_idx, score_vec, frac=0.05):
    nlib = len(score_vec); npos = len(pos_idx)
    if npos == 0: return 0.0
    order = np.argsort(-score_vec)
    K = max(1, int(round(nlib * frac)))
    pos_set = set(pos_idx)
    hits = sum(1 for i in order[:K] if i in pos_set)
    return (hits / K) / (npos / nlib)

# ── 逐靶点闭环 (序列 L2) ──
l2 = Layer2BindingDBSeq()
log("L2 模型: BindingDB-L2-SeqCNN (128维序列嵌入)")
fracs = [0.01, 0.05, 0.10]
pairs_total = 0
records = []

for cid, v in panel.items():
    pidx = [lib_canon[c] for c in (canon(p) for p in v.get("pos", [])) if c and c in lib_canon]
    nidx = [lib_canon[c] for c in (canon(n) for n in v.get("neg", [])) if c and c in lib_canon]
    if len(pidx) < 3:
        continue
    l2p = l2.predict_matrix(molfeat, chembl_id=cid)      # (N,)
    final = 0.20 * L1_arr + 0.50 * l2p + 0.20 * L3_arr

    order_final = np.argsort(-final)
    order_l1 = np.argsort(-L1_arr)
    order_l2 = np.argsort(-l2p)
    rf = recall_at(pidx, order_final, fracs)
    rl1 = recall_at(pidx, order_l1, fracs)
    rl2 = recall_at(pidx, order_l2, fracs)
    ef = ef_at(pidx, final, 0.05)
    neg_mask = np.ones(N, dtype=bool); neg_mask[pidx] = False
    bed = bedroc(final[pidx], final[neg_mask], 20.0)
    y = np.zeros(N, dtype=int); y[pidx] = 1
    try: auc = float(roc_auc_score(y, final))
    except Exception: auc = float("nan")
    # L2-only AUC (直接看靶点条件化质量)
    try: auc_l2 = float(roc_auc_score(y, l2p))
    except Exception: auc_l2 = float("nan")

    s1 = np.where(L1_arr >= 0.45)[0]
    s1_pos = sum(1 for i in pidx if i in set(s1.tolist()))
    if len(s1) > 0:
        s1order = s1[np.argsort(-l2p[s1])]
        K2 = max(1, int(round(len(s1order) * 0.10)))
        s2 = set(s1order[:K2].tolist())
        s2_pos = sum(1 for i in pidx if i in s2)
    else:
        s2_pos = 0
    s3 = set(np.where(L3_arr >= 0.5)[0].tolist())
    s3_pos = sum(1 for i in pidx if i in s3)

    pairs_total += N
    records.append({
        "target": cid, "name": v["name"],
        "n_pos_in_lib": len(pidx), "n_neg_in_lib": len(nidx),
        "final_recall": rf, "l1_recall": rl1, "l2_recall": rl2,
        "ef5": ef, "bedroc": round(bed, 3),
        "auc_full": round(auc, 3), "auc_l2_only": round(auc_l2, 3),
        "cascade": {
            "l1_keep_frac": round(len(s1) / N, 3), "l1_pos_kept": round(s1_pos / len(pidx), 3),
            "l2_top10_pos_kept": round(s2_pos / len(pidx), 3),
            "l3_keep_frac": round(len(s3) / N, 3), "l3_pos_kept": round(s3_pos / len(pidx), 3),
        },
    })

wall = time.time() - t0
log(f"闭环评估完成 {len(records)} 靶点; 配对 {pairs_total:,}; 墙钟 {wall/60:.1f}min; 吞吐 {pairs_total/wall:,.0f} pair/s")

# ── 聚合 ──
def agg(fn):
    vals = [fn(r) for r in records]
    vals = [x for x in vals if x is not None and not (isinstance(x, float) and math.isnan(x))]
    if not vals: return None
    vals.sort(); n = len(vals)
    return {"n": n, "median": round(vals[n//2], 3), "mean": round(sum(vals)/n, 3),
            "min": round(vals[0], 3), "max": round(vals[-1], 3)}
def agg_recall(field):
    return {f: agg(lambda r: r[field][f]) for f in ["recall@1%", "recall@5%", "recall@10%"]}

summary = {
    "setup": {"library_size": N, "n_targets_evaluated": len(records),
              "n_pos_total": sum(r["n_pos_in_lib"] for r in records),
              "l2_model": "Layer2BindingDBSeq (128维序列嵌入 CNN, 端到端 MLP)"},
    "throughput": {"pairs_total": pairs_total, "wall_seconds": round(wall, 1),
                   "pairs_per_second": round(pairs_total / wall, 1)},
    "recall_final": agg_recall("final_recall"),
    "recall_l2_only": agg_recall("l2_recall"),
    "recall_l1_baseline": agg_recall("l1_recall"),
    "ef5_final": agg(lambda r: r["ef5"]),
    "bedroc_final": agg(lambda r: r["bedroc"]),
    "auc_full_final": agg(lambda r: r["auc_full"]),
    "auc_l2_only": agg(lambda r: r["auc_l2_only"]),
    "frac_targets_recall1pct_ge_0.5": round(
        sum(1 for r in records if r["final_recall"]["recall@1%"] >= 0.5) / len(records), 3),
}
out = {"summary": summary, "per_target": records,
       "meta": {"script": "hts_closed_loop_seq.py", "l1l3_cache": str(CACHE)}}
json.dump(out, open(MTB / "hts_closed_loop_seq_results.json", "w"), indent=1)

# ── 与旧哈希版对比 ──
cmp_lines = []
old_path = MTB / "hts_closed_loop_results.json"
if old_path.exists():
    old = json.load(open(old_path))
    old_auc = {r["target"]: r.get("auc_full") for r in old["per_target"]}
    reversed_fixed = 0; reversed_total = 0; improved = 0
    for r in records:
        oa = old_auc.get(r["target"])
        na = r["auc_full"]
        if oa is None or na is None: continue
        if oa < 0.5:
            reversed_total += 1
            if na >= 0.5: reversed_fixed += 1
        if na > oa: improved += 1
    old_med = old["summary"]["auc_full_final"]["median"]
    new_med = summary["auc_full_final"]["median"]
    cmp_lines = [
        f"旧哈希AUC中位={old_med} -> 新序列AUC中位={new_med}",
        f"AUC提升的靶点: {improved}/{len(records)}",
        f"旧反向(AUC<0.5)靶点: {reversed_total} 个, 其中修复(新>=0.5): {reversed_fixed} 个",
    ]
    summary["compare_vs_hash"] = {
        "old_auc_median": old_med, "new_auc_median": new_med,
        "n_improved": improved, "n_targets": len(records),
        "reversed_old": reversed_total, "reversed_fixed": reversed_fixed,
    }
    json.dump(out, open(MTB / "hts_closed_loop_seq_results.json", "w"), indent=1)

print("\n" + "=" * 70)
print("高通量测试闭环结果 (序列嵌入 L2 升级版)")
print("=" * 70)
print(f"库 {N} × {len(records)} 靶点 = {pairs_total:,} 配对 | 吞吐 {pairs_total/wall:,.0f} pair/s")
print("\n[全库回收率 recall@K]  median")
for label, key in [("final(四级漏斗)", "recall_final"),
                   ("L2-only(靶点感知)", "recall_l2_only"),
                   ("L1-only(质量基线)", "recall_l1_baseline")]:
    a = summary[key]
    print(f"  {label:22s} @1%={a['recall@1%']['median']:.3f}  "
          f"@5%={a['recall@5%']['median']:.3f}  @10%={a['recall@10%']['median']:.3f}")
print(f"\nEF@5% median={summary['ef5_final']['median']}  "
      f"BEDROC median={summary['bedroc_final']['median']}")
print(f"全库AUC median={summary['auc_full_final']['median']}  "
      f"L2-only AUC median={summary['auc_l2_only']['median']}")
print(f"达到 recall@1%>=0.5 的靶点比例: {summary['frac_targets_recall1pct_ge_0.5']}")
if cmp_lines:
    print("\n[对比旧哈希版]")
    for l in cmp_lines: print("  " + l)
print("done.")

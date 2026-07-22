# -*- coding: utf-8 -*-
"""
用扩大训练面后的 L2 (BindingDB+ChEMBL+OpenTargets) 在 (24k候选 × 164靶点) 重跑,
验证 AUC 是否提升, 并单独看此前 48 个未映射靶点的改善。

改进点 (相对上一轮 rerun_panel_real_l2.py):
  - 靶点解析改用 target_resolver.TargetResolver:
      面板 target_chembl_id -> ChEMBL target_text 精确命中 (覆盖那 48 个)
  - 报告分组: 全量 / 此前 mapped / 此前 unmapped(现应 mapped)
"""
import json, csv, random, sys, time
from pathlib import Path
from collections import Counter, defaultdict

import numpy as np
from sklearn.metrics import roc_auc_score

REPO = Path("<validated-workspace>")
SC_DIR = REPO / "评分_work_package/评分"
sys.path.insert(0, str(SC_DIR))
from l2_bindingdb import BindingDBFeature, Layer2BindingDB
from scoring import MoleculeScorer, Layer1Scorer, Layer3Scorer
from target_resolver import TargetResolver

random.seed(0)
t0 = time.time()

panel = json.load(open(REPO / "scientific_validation/multitarget_benchmark/target_panel.json"))
lib_rows = []
with open(REPO / "scientific_validation/multitarget_benchmark/candidates_10k.csv") as f:
    for row in csv.DictReader(f):
        lib_rows.append((row["id"], row["smiles"].strip()))
lib_smiles = [s for _, s in lib_rows]
print(f"[{time.time()-t0:.0f}s] 候选库 {len(lib_smiles)} 药; 面板 {len(panel)} 靶点")

# 解析器: 面板 chembl_id -> ChEMBL target_text
resolver = TargetResolver()
mapping = {}
methods = {}
for cid, v in panel.items():
    text, method = resolver.resolve(v["name"], chembl_id=cid)
    mapping[cid] = text
    methods[cid] = method
n_mapped = sum(1 for t in mapping.values() if t)
print(f"[{time.time()-t0:.0f}s] 解析到 target_text: {n_mapped}/{len(panel)}  (此前 116)")

# 共享 decoy 池
panel_smiles = set()
for v in panel.values():
    panel_smiles.update(v.get("pos", [])); panel_smiles.update(v.get("neg", []))
pool_candidates = [s for s in lib_smiles if s not in panel_smiles]
decoy_pool = random.sample(pool_candidates, 500)
print(f"[{time.time()-t0:.0f}s] decoy 池 {len(decoy_pool)}")

feat = BindingDBFeature()
l2 = Layer2BindingDB(prefer="mlp")
sc = MoleculeScorer(use_unimol=False, l2_method="bindingdb", default_target_text="")
l1sc = sc.l1
l3sc = sc.l3

need = set(decoy_pool)
for v in panel.values():
    need.update(v.get("pos", [])); need.update(v.get("neg", []))
need = list(need)
print(f"[{time.time()-t0:.0f}s] 需 L1/L3 唯一分子 {len(need)}")
mol_feat = {}; L1c, L3c = {}, {}
for i, smi in enumerate(need):
    mol_feat[smi] = feat.mol_features(smi)
    try: L1c[smi] = float(l1sc.score(smi).get("layer1_score") or 0.0)
    except Exception: L1c[smi] = 0.0
    try:
        a = l3sc.score(smi).get("admet_score")
        L3c[smi] = float(a) if a is not None else 0.5
    except Exception: L3c[smi] = 0.5
    if (i+1) % 1000 == 0:
        print(f"[{time.time()-t0:.0f}s] L1/L3 {i+1}/{len(need)}")
print(f"[{time.time()-t0:.0f}s] L1/L3 完成; 无效 {sum(1 for v in mol_feat.values() if v is None)}")

tt_cache = {tt: feat.target_features(tt) for tt in set(mapping.values()) if tt}

def auc_safe(y, s):
    y = np.asarray(y); s = np.asarray(s)
    if len(np.unique(y)) < 2:
        return None
    return float(roc_auc_score(y, s))

def eval_target(cid):
    tt = mapping[cid]
    if not tt:
        return None
    v = panel[cid]
    pos = [s for s in v.get("pos", []) if mol_feat.get(s) is not None]
    neg = [s for s in v.get("neg", []) if mol_feat.get(s) is not None]
    if len(pos) < 3 or len(neg) < 3:
        return None
    cands = pos + neg + decoy_pool
    X = np.array([np.concatenate([mol_feat[s], tt_cache[tt]]) for s in cands], dtype=np.float32)
    l2p = l2.predict_proba(X)
    l1 = np.array([L1c[s] for s in cands])
    l3 = np.array([L3c[s] for s in cands])
    final = 0.20 * l1 + 0.50 * l2p + 0.20 * l3
    y = [1]*len(pos) + [0]*(len(neg)+len(decoy_pool))
    return {"target": cid, "name": v["name"], "method": methods[cid],
            "n_pos": len(pos), "n_neg": len(neg),
            "auc_l2": auc_safe(y, l2p), "auc_final": auc_safe(y, final),
            "auc_l1_baseline": auc_safe(y, l1)}

results = []
for cid in panel:
    r = eval_target(cid)
    if r:
        results.append(r)
print(f"[{time.time()-t0:.0f}s] 评估完成 {len(results)} 靶点")

def summarize(lst, key):
    vals = [r[key] for r in lst if r[key] is not None]
    if not vals: return None
    vals.sort(); n = len(vals)
    return {"n": n, "median": round(vals[n//2],3), "mean": round(sum(vals)/n,3),
            "min": round(vals[0],3), "max": round(vals[-1],3),
            "frac_gt_0.5": round(sum(1 for x in vals if x>0.5)/n,3),
            "frac_gt_0.7": round(sum(1 for x in vals if x>0.7)/n,3)}

# 区分"此前未映射组": method 为 chembl_id_exact 且不在上一轮 bindingdb 映射中的,
# 简化: 报告 by method
by_method = defaultdict(list)
for r in results:
    by_method[r["method"]].append(r)

summary = {"all_l2": summarize(results,"auc_l2"),
           "all_final": summarize(results,"auc_final"),
           "all_l1_baseline": summarize(results,"auc_l1_baseline"),
           "by_method": {m: summarize(lst,"auc_l2") for m,lst in by_method.items()}}
out = {"summary": summary, "per_target": results,
       "meta": {"model": "combined BindingDB+ChEMBL+OpenTargets MLP(256,128)",
                "n_mapped": n_mapped, "n_panel": len(panel),
                "decoy_pool": len(decoy_pool)}}
json.dump(out, open(REPO/"scientific_validation/multitarget_benchmark/real_l2_combined_rerun.json","w"), indent=1)

print(f"\n[{time.time()-t0:.0f}s] === 结果 (扩大训练面后) ===")
print("全量 L2 AUC      :", summary["all_l2"])
print("全量 集成final AUC:", summary["all_final"])
print("全量 L1基线 AUC  :", summary["all_l1_baseline"])
for m, s in summary["by_method"].items():
    print(f"  method={m:20s} L2 AUC: {s}")
print("done.")

# -*- coding: utf-8 -*-
"""
用升级后的真实 L2 (BindingDB MLP, 已接入 MoleculeScorer L2槽位权重0.50) 在
同一组 (24k候选药库 × 164靶点面板) 上重跑, 验证 AUC 回到 0.5 以上。

效率策略:
  - 共享 decoy 池 (500个随机库药) 作为各靶点硬阴性, 避免每靶点重采样爆炸。
  - 每个唯一 SMILES 只算一次 L1/L3 并缓存; L2 用批量推理。
  - 靶点文本: 用面板 name 对齐到 BindingDB target_text (全串+尾部规范名两遍匹配)。
    mapped 靶点(模型训练时见过) 为主结果; unmapped 用 name 作文本(弱信号)单列。

产出: scientific_validation/multitarget_benchmark/real_l2_rerun_results.json
"""
import json, csv, random, sys, time, re
from pathlib import Path
from collections import Counter

import numpy as np
from sklearn.metrics import roc_auc_score

REPO = Path("<validated-workspace>")
SC_DIR = REPO / "评分_work_package/评分"
sys.path.insert(0, str(SC_DIR))
from l2_bindingdb import BindingDBFeature, Layer2BindingDB
from scoring import MoleculeScorer, Layer1Scorer, Layer3Scorer

random.seed(0)
t0 = time.time()

# ---------- 数据 ----------
panel = json.load(open(REPO / "scientific_validation/multitarget_benchmark/target_panel.json"))
lib_rows = []
with open(REPO / "scientific_validation/multitarget_benchmark/candidates_10k.csv") as f:
    for row in csv.DictReader(f):
        lib_rows.append((row["id"], row["smiles"].strip()))
lib_smiles = [s for _, s in lib_rows]
print(f"[{time.time()-t0:.0f}s] 候选库 {len(lib_smiles)} 药; 面板 {len(panel)} 靶点")

# ---------- 靶点文本映射 (面板name -> BindingDB target_text) ----------
ex_csv = REPO / "data_lake/bindingdb/aligned_model_input/bindingdb_202606_target_match_examples.csv"
tc = Counter()
with open(ex_csv) as f:
    for row in csv.DictReader(f):
        tc[row["target_text"]] += 1
print(f"[{time.time()-t0:.0f}s] BindingDB 不同靶点文本 {len(tc)}")

UNIPROT = re.compile(r"\b[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}\b")
PDB = re.compile(r"\b[0-9][A-Z0-9]{3}\b")
def trailing_name(text):
    # 取 UniProt 之后、PDB 码之前的片段作为规范名
    m = UNIPROT.search(text)
    seg = text[m.end():] if m else text
    # 去掉尾随 PDB 码
    seg = PDB.sub("", seg)
    return seg.strip(" ,")

trailing = {}
for t in tc:
    tr = trailing_name(t)
    if tr and (tr not in trailing or tc[t] > tc[trailing[tr]]):
        trailing[tr] = t

def map_target(name):
    ln = name.lower()
    cands = [(t, tc[t]) for t in tc if ln in t.lower()]
    if not cands:
        cands = [(trailing[tr], tc[trailing[tr]]) for tr in trailing
                 if ln in tr.lower() and tr]
    if not cands:
        return None
    cands.sort(key=lambda x: -x[1])
    return cands[0][0]

mapping = {tid: map_target(v["name"]) for tid, v in panel.items()}
mapped = {tid: t for tid, t in mapping.items() if t}
print(f"[{time.time()-t0:.0f}s] 映射到 BindingDB 靶点文本: {len(mapped)}/{len(panel)}")

# ---------- 共享 decoy 池 ----------
panel_smiles = set()
for v in panel.values():
    panel_smiles.update(v.get("pos", [])); panel_smiles.update(v.get("neg", []))
pool_candidates = [s for s in lib_smiles if s not in panel_smiles]
decoy_pool = random.sample(pool_candidates, 500)
print(f"[{time.time()-t0:.0f}s] decoy 池 {len(decoy_pool)} (排除面板已知药避免泄漏)")

# ---------- 特征/打分器 ----------
feat = BindingDBFeature()
l2 = Layer2BindingDB(prefer="mlp")
sc = MoleculeScorer(use_unimol=False, l2_method="bindingdb", default_target_text="")
l1sc = sc.l1
l3sc = sc.l3

# 收集所有需要的唯一 SMILES
need = set(decoy_pool)
for v in panel.values():
    need.update(v.get("pos", [])); need.update(v.get("neg", []))
need = list(need)
print(f"[{time.time()-t0:.0f}s] 需计算 L1/L3 的唯一分子 {len(need)}")

mol_feat = {}
L1c, L3c = {}, {}
for i, smi in enumerate(need):
    mf = feat.mol_features(smi)
    mol_feat[smi] = mf  # None 表示无效
    try:
        L1c[smi] = float(l1sc.score(smi).get("layer1_score") or 0.0)
    except Exception:
        L1c[smi] = 0.0
    try:
        a = l3sc.score(smi).get("admet_score")
        L3c[smi] = float(a) if a is not None else 0.5
    except Exception:
        L3c[smi] = 0.5
    if (i + 1) % 1000 == 0:
        print(f"[{time.time()-t0:.0f}s] L1/L3 进度 {i+1}/{len(need)}")
print(f"[{time.time()-t0:.0f}s] L1/L3 完成; 无效分子 {sum(1 for v in mol_feat.values() if v is None)}")

# 靶点文本特征缓存
tt_cache = {tt: feat.target_features(tt) for tt in set(mapping.values())}

def auc_safe(y, s):
    y = np.asarray(y); s = np.asarray(s)
    if len(np.unique(y)) < 2:
        return None
    return float(roc_auc_score(y, s))

# ---------- 逐靶点评估 ----------
def eval_target(tid, tt):
    v = panel[tid]
    pos = [s for s in v.get("pos", []) if mol_feat.get(s) is not None]
    neg = [s for s in v.get("neg", []) if mol_feat.get(s) is not None]
    if len(pos) < 3 or len(neg) < 3:
        return None
    cands = pos + neg + decoy_pool
    X = np.array([np.concatenate([mol_feat[s], tt_cache[tt]]) for s in cands], dtype=np.float32)
    l2p = l2.predict_proba(X)
    l1 = np.array([L1c[s] for s in cands])
    l3 = np.array([L3c[s] for s in cands])
    final = 0.20 * l1 + 0.50 * l2p + 0.20 * l3  # MoleculeScorer 集成公式
    y = [1] * len(pos) + [0] * (len(neg) + len(decoy_pool))
    return {
        "target": tid, "name": v["name"],
        "n_pos": len(pos), "n_neg": len(neg), "n_decoy": len(decoy_pool),
        "auc_l2": auc_safe(y, l2p),
        "auc_final": auc_safe(y, final),
        "auc_l1_baseline": auc_safe(y, l1),  # 无L2的"瞎选"基线
    }

results = {"mapped": [], "unmapped": []}
for tid, tt in mapping.items():
    r = eval_target(tid, tt)
    if r is None:
        continue
    (results["mapped"] if tt else results["unmapped"]).append(r)

def summarize(lst, key):
    vals = [r[key] for r in lst if r[key] is not None]
    if not vals:
        return None
    vals.sort()
    n = len(vals)
    med = vals[n // 2]
    gt5 = sum(1 for x in vals if x > 0.5) / n
    gt7 = sum(1 for x in vals if x > 0.7) / n
    return {"n": n, "median": round(med, 3), "mean": round(sum(vals) / n, 3),
            "min": round(vals[0], 3), "max": round(vals[-1], 3),
            "frac_gt_0.5": round(gt5, 3), "frac_gt_0.7": round(gt7, 3)}

summary = {
    "mapped_l2": summarize(results["mapped"], "auc_l2"),
    "mapped_final": summarize(results["mapped"], "auc_final"),
    "mapped_l1_baseline": summarize(results["mapped"], "auc_l1_baseline"),
    "unmapped_l2": summarize(results["unmapped"], "auc_l2"),
    "unmapped_final": summarize(results["unmapped"], "auc_final"),
}
out = {"summary": summary, "per_target": results,
       "meta": {"decoy_pool": len(decoy_pool), "model": "BindingDB MLP(256,128)",
                "n_mapped": len(results["mapped"]), "n_unmapped": len(results["unmapped"]),
                "trained_auc_mlp": 0.9354}}
json.dump(out, open(REPO / "scientific_validation/multitarget_benchmark/real_l2_rerun_results.json", "w"), indent=1)

print(f"\n[{time.time()-t0:.0f}s] === 结果 ===")
print("映射靶点(模型见过) L2 AUC :", summary["mapped_l2"])
print("映射靶点(模型见过) 集成final AUC :", summary["mapped_final"])
print("映射靶点 L1基线(无L2) AUC :", summary["mapped_l1_baseline"])
if summary["unmapped_l2"]:
    print("未映射靶点 L2 AUC :", summary["unmapped_l2"])
print("done.")

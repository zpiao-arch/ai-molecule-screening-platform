# -*- coding: utf-8 -*-
"""
② 柔性/刚性对接精排 (Docking Rerank) — 靶点 CHEMBL2051 神经氨酸酶 (Neuraminidase)
================================================================================
动机: 面板闭环显示 L2 (无论哈希 or 序列嵌入) 对 NA 的 top-K 富集很弱 (recall@10%≈0)。
本实验用真实物理对接 (smina + AutoDock Vina 打分函数, 3TI6 受体) 对
  21 个库内已知 NA 活性分子 + 100 个随机诱饵 (decoy)
逐一 SMILES->3D(obabel gen3d, pH7.4 质子化)->pdbqt->刚性对接, 取最优 pose affinity。
比较三种排序对活性分子的早期富集:
  (a) L2 靶点感知分数 (哈希版)
  (b) 对接 affinity (越负越好 -> 用 -affinity 排序)
  (c) L2 与 docking 的 rank 融合
指标: ROC-AUC, EF@10%, top-K recall。证明 "对接精排把弱 L2 信号转成真实 top-K 命中"。
"""
import json, csv, sys, os, subprocess, time, random, tempfile
from pathlib import Path
import numpy as np

REPO = Path("<validated-workspace>")
MTB = REPO / "scientific_validation/multitarget_benchmark"
SC = REPO / "评分_work_package/评分"
SMINA_BIN = REPO / "ai_drug_eval_tools/micromamba-root/envs/smina-local/bin/smina"
OBABEL_BIN = REPO / "ai_drug_eval_tools/micromamba-root/envs/smina-local/bin/obabel"
RECEPTOR = REPO / "ai_mol_loop/influenza_na_2000_project_20260625/stage4/receptors/3TI6_protein_only_obabel.pdbqt"
BOX = dict(cx=-28.914, cy=14.334, cz=20.794, sx=23.585, sy=20.45, sz=24.18)
WORK = MTB / "dock_rerank_work"
WORK.mkdir(exist_ok=True)
OUT = MTB / "dock_rerank_na_results.json"
CID = "CHEMBL2051"
N_DECOY = 100
SEED = 42

t0 = time.time()
def log(*a): print(f"[{time.time()-t0:6.0f}s]", *a, flush=True)

sys.path.insert(0, str(SC))
from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

def canon(s):
    try: return Chem.MolToSmiles(Chem.MolFromSmiles(s))
    except Exception: return None

# ── 1. 构建配体集 (21 活性 + 100 诱饵) ──
panel = json.load(open(MTB / "target_panel.json"))
v = panel[CID]
lib = []
with open(MTB / "candidates_10k.csv") as f:
    for r in csv.DictReader(f):
        lib.append((r["id"], r["smiles"].strip()))
libc = {}
for i, (lid, smi) in enumerate(lib):
    c = canon(smi)
    if c: libc[c] = i
active_idx = set()
actives = []
for p in v.get("pos", []):
    c = canon(p)
    if c and c in libc:
        i = libc[c]; active_idx.add(i)
        actives.append((lib[i][0], lib[i][1]))
random.seed(SEED)
pool = [i for i in range(len(lib)) if i not in active_idx]
decoy_i = random.sample(pool, N_DECOY)
decoys = [(lib[i][0], lib[i][1]) for i in decoy_i]
ligands = [(lid, smi, 1) for lid, smi in actives] + [(lid, smi, 0) for lid, smi in decoys]
log(f"配体集: 活性 {len(actives)} + 诱饵 {len(decoys)} = {len(ligands)}")

# ── 2. L2 分数 (哈希版, 面板上更优) ──
from l2_bindingdb import Layer2BindingDB, BindingDBFeature
from target_resolver import TargetResolver
res = TargetResolver()
tt, _ = res.resolve(v["name"], chembl_id=CID)
l2 = Layer2BindingDB(prefer="mlp")
feat = BindingDBFeature()
tf = feat.target_features(tt)
l2_scores = {}
for lid, smi, lab in ligands:
    mf = feat.mol_features(smi)
    if mf is None:
        l2_scores[lid] = 0.0; continue
    x = np.hstack([mf, tf]).reshape(1, -1)
    l2_scores[lid] = float(l2.predict_proba(x)[0])
log("L2 (哈希) 打分完成")

# ── 3. 逐配体 prep + dock ──
def prep(smi, out_pdbqt):
    r = subprocess.run([str(OBABEL_BIN), f"-:{smi}", "-O", str(out_pdbqt),
                        "--gen3d", "-p", "7.4"],
                       capture_output=True, text=True, timeout=90)
    return out_pdbqt.exists() and out_pdbqt.stat().st_size > 0

def dock(lig_pdbqt, out_pose):
    cmd = [str(SMINA_BIN), "--receptor", str(RECEPTOR), "--ligand", str(lig_pdbqt),
           "--center_x", str(BOX["cx"]), "--center_y", str(BOX["cy"]), "--center_z", str(BOX["cz"]),
           "--size_x", str(BOX["sx"]), "--size_y", str(BOX["sy"]), "--size_z", str(BOX["sz"]),
           "--exhaustiveness", "8", "--cpu", "4", "--num_modes", "3", "--seed", str(SEED),
           "--out", str(out_pose)]
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    best = None
    for line in p.stdout.splitlines():
        s = line.split()
        if len(s) >= 2 and s[0] == "1":
            try: best = float(s[1]); break
            except ValueError: pass
    return best

records = []
for k, (lid, smi, lab) in enumerate(ligands, 1):
    lig = WORK / f"{lid}.pdbqt"; pose = WORK / f"{lid}_pose.pdbqt"
    aff = None; status = "ok"
    try:
        if prep(smi, lig):
            aff = dock(lig, pose)
            if aff is None: status = "dock_no_score"
        else:
            status = "prep_failed"
    except subprocess.TimeoutExpired:
        status = "timeout"
    except Exception as e:
        status = f"err:{type(e).__name__}"
    records.append({"id": lid, "label": lab, "smiles": smi,
                    "l2": round(l2_scores.get(lid, 0.0), 4),
                    "affinity": aff, "status": status})
    if k % 10 == 0 or k == len(ligands):
        ok = sum(1 for r in records if r["affinity"] is not None)
        log(f"进度 {k}/{len(ligands)}  成功对接 {ok}")

# ── 4. 指标 ──
from sklearn.metrics import roc_auc_score
docked = [r for r in records if r["affinity"] is not None]
y = np.array([r["label"] for r in docked])
l2v = np.array([r["l2"] for r in docked])
dockv = np.array([-r["affinity"] for r in docked])  # 越负affinity -> 越大分数
# rank 融合 (两者 rank 平均, 越小越好 -> 取负)
def ranks(v):  # 高分->低rank(好)
    order = np.argsort(-v); rk = np.empty_like(order); rk[order] = np.arange(len(v)); return rk
comb = -(ranks(l2v) + ranks(dockv)).astype(float)

def ef_at(score, frac=0.10):
    n = len(score); npos = int(y.sum())
    if npos == 0: return 0.0
    order = np.argsort(-score); K = max(1, int(round(n*frac)))
    hits = int(y[order[:K]].sum())
    return round((hits/K)/(npos/n), 3)

def recall_topk(score, K):
    order = np.argsort(-score)
    return round(int(y[order[:K]].sum())/int(y.sum()), 3)

metrics = {}
for name, sc in [("L2_hash", l2v), ("docking", dockv), ("L2+docking", comb)]:
    try: auc = round(float(roc_auc_score(y, sc)), 3)
    except Exception: auc = None
    metrics[name] = {"auc": auc, "ef@10%": ef_at(sc, 0.10),
                     "recall@10": recall_topk(sc, 10), "recall@20": recall_topk(sc, 20),
                     "recall@30": recall_topk(sc, 30)}

summary = {
    "target": CID, "target_name": v["name"], "receptor": "3TI6",
    "n_actives_docked": int(y.sum()), "n_decoys_docked": int((1-y).sum()),
    "n_prep_failed": sum(1 for r in records if r["status"] != "ok"),
    "docking_engine": "smina 2020.12.10 (Vina scoring)", "box": BOX,
    "metrics": metrics,
    "active_affinity_median": round(float(np.median([-x for x, l in zip(dockv, y) if l == 1])), 2),
    "decoy_affinity_median": round(float(np.median([-x for x, l in zip(dockv, y) if l == 0])), 2),
}
json.dump({"summary": summary, "records": records}, open(OUT, "w"), indent=1)

log("对接精排完成")
print("\n" + "="*68)
print(f"② 对接精排结果 — {v['name']} ({CID}) / 受体 3TI6")
print("="*68)
print(f"成功对接: 活性 {summary['n_actives_docked']}/{len(actives)}  诱饵 {summary['n_decoys_docked']}/{len(decoys)}")
print(f"活性 affinity 中位 {summary['active_affinity_median']} kcal/mol | 诱饵 {summary['decoy_affinity_median']} kcal/mol")
print(f"\n{'排序方式':14s} {'AUC':>6s} {'EF@10%':>7s} {'R@10':>6s} {'R@20':>6s} {'R@30':>6s}")
for name in ["L2_hash", "docking", "L2+docking"]:
    m = metrics[name]
    print(f"{name:14s} {str(m['auc']):>6s} {str(m['ef@10%']):>7s} "
          f"{str(m['recall@10']):>6s} {str(m['recall@20']):>6s} {str(m['recall@30']):>6s}")
print("done.")

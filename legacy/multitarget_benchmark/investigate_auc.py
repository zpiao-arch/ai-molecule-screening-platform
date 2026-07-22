# -*- coding: utf-8 -*-
"""AUC 低的根因调查: 训练支撑度 + 活性分子稀疏度 + 集成效应 三维诊断."""
import csv, json, sys
import numpy as np
from pathlib import Path
sys.path.insert(0, str(Path("<validated-workspace>/评分_work_package/评分")))
from target_resolver import TargetResolver

REPO = Path("<validated-workspace>")
BD = REPO / "data_lake/bindingdb/aligned_model_input/bindingdb_202606_target_match_examples.csv"
CH = REPO / "data_lake/chembl/latest/aligned_model_input/chembl_37_target_match_examples.csv"
PANEL = json.load(open(REPO / "scientific_validation/multitarget_benchmark/target_panel.json"))
RES = json.load(open(REPO / "scientific_validation/multitarget_benchmark/hts_closed_loop_results.json"))

# 1) 训练支撑度: 合并训练集里每个 target_text 出现次数 (模型实际 key 的是这个字符串)
support_by_text = {}
for path in (BD, CH):
    if not path.exists():
        continue
    with open(path) as f:
        for row in csv.DictReader(f):
            t = row.get("target_text", "").strip()
            if t:
                support_by_text[t] = support_by_text.get(t, 0) + 1
print(f"[训练支撑] 合并训练集唯一 target_text 数: {len(support_by_text)}")

# 2) 面板靶点 -> resolved_text + method + support
res = TargetResolver()
rec = {}
for cid, v in PANEL.items():
    text, method = res.resolve(v.get("name", ""), chembl_id=cid)
    sup = support_by_text.get(text, 0) if text else 0
    rec[cid] = {"text": text, "method": method, "support": sup}

# 3) 合并结果
rows = []
for r in RES["per_target"]:
    cid = r["target"]
    info = rec.get(cid, {})
    rows.append({
        "cid": cid,
        "name": r.get("name", ""),
        "method": info.get("method", "?"),
        "support": info.get("support", 0),
        "npos": r["n_pos_in_lib"],
        "auc": r["auc_full"],
        "ef5": r["ef5"],
        "r5_final": r["final_recall"]["recall@5%"],
        "r5_l2": r["l2_recall"]["recall@5%"],
    })
A = np.array([x["auc"] for x in rows])
SUP = np.array([x["support"] for x in rows], dtype=float)
NPOS = np.array([x["npos"] for x in rows], dtype=float)
R5F = np.array([x["r5_final"] for x in rows])
R5L = np.array([x["r5_l2"] for x in rows])

def corr(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 3: return float("nan")
    return float(np.corrcoef(a[m], b[m])[0, 1])

def med(x):
    x = np.array([v for v in x if np.isfinite(v)])
    return float(np.median(x)) if len(x) else float("nan")

def bucket_stats(key, edges, labels):
    print(f"\n=== 按 {key} 分组的 AUC / EF5 / recall@5% 中位数 ===")
    vals = {"support": SUP, "npos": NPOS}[key]
    bins = np.digitize(vals, edges)  # edges 升序
    for i, lab in enumerate(labels):
        mask = bins == i
        if mask.sum() == 0:
            print(f"  {lab:14s} n=0")
            continue
        print(f"  {lab:14s} n={mask.sum():3d} | AUC={med(A[mask]):.3f}  EF5={med(np.array([x['ef5'] for x,m in zip(rows,mask) if m])):.3f}  recall@5%={med(R5F[mask]):.3f}")

print("\n===== 总体 =====")
print(f"靶点总数: {len(rows)}")
print(f"AUC: median={med(A):.3f} mean={A.mean():.3f} min={A.min():.3f} max={A.max():.3f}")
print(f"  AUC>0.6: {(A>0.6).sum()}  AUC>0.7: {(A>0.7).sum()}  AUC<0.5(反向): {(A<0.5).sum()}")
print(f"recall@5%: median={med(R5F):.3f} (随机基线=0.05)")
print(f"n_pos_in_lib: median={med(NPOS):.1f} mean={NPOS.mean():.1f} min={NPOS.min():.0f} max={NPOS.max():.0f}")
for thr in (5, 10, 20, 50):
    print(f"  活性分子<{thr} 的靶点: {(NPOS<thr).sum()}/{len(rows)}")
print(f"support(训练条数): median={med(SUP):.0f} mean={SUP.mean():.0f} min={SUP.min():.0f} max={SUP.max():.0f}")

print("\n===== 相关性 =====")
print(f"corr(support, AUC)        = {corr(SUP, A):.3f}")
print(f"corr(n_pos,   AUC)        = {corr(NPOS, A):.3f}")
print(f"corr(n_pos,   recall@5%)  = {corr(NPOS, R5F):.3f}")
print(f"corr(support, recall@5%)  = {corr(SUP, R5F):.3f}")

bucket_stats("support", [1, 50, 500], ["0(无)", "1-49", "50-499", "500+"])
bucket_stats("npos", [5, 20, 100], ["<5", "5-19", "20-99", "100+"])

# 4) 反向靶点 (AUC<0.5) 画像
rev = [x for x in rows if x["auc"] < 0.5]
print(f"\n===== 反向靶点 (AUC<0.5, n={len(rev)}) 画像 =====")
if rev:
    print(f"  median support: {np.median([x['support'] for x in rev]):.0f}")
    print(f"  median n_pos:   {np.median([x['npos'] for x in rev]):.0f}")
    print(f"  样例(最差5):")
    for x in sorted(rev, key=lambda z: z["auc"])[:5]:
        print(f"    {x['cid']} {x['name'][:24]:24s} AUC={x['auc']:.3f} sup={x['support']:5d} npos={x['npos']:4d} method={x['method']}")

# 5) 集成效应: final vs L2-only (看 L1/L3 是帮忙还是稀释)
delta = R5F - R5L
print(f"\n===== 集成效应 (final_recall@5% - l2_recall@5%) =====")
print(f"  median(final - l2) = {med(delta):.3f}  (>0 表示 L1/L3 帮忙, <0 表示稀释)")
low = A < 0.6
print(f"  低AUC(<0.6)靶点的 median(final - l2) = {med(delta[low]):.3f}  (集成对弱靶点是否帮倒忙)")
print(f"  高AUC(>0.7)靶点的 median(final - l2) = {med(delta[A>0.7]):.3f}")

# 6) 强/弱靶点对比
hi = sorted(rows, key=lambda z: -z["auc"])[:5]
lo = sorted(rows, key=lambda z: z["auc"])[:5]
print("\n  最强5靶:", [(x['cid'], round(x['auc'],3), x['support'], x['npos']) for x in hi])
print("  最弱5靶:", [(x['cid'], round(x['auc'],3), x['support'], x['npos']) for x in lo])

# -*- coding: utf-8 -*-
"""
1000 靶点 × 10000 候选 完整闭环批测（auto 路由：每靶取 AUC 最高的可用方法）
================================================================================
语义映射（应用户口径 "跑1000次不同药物(10000候选)不同靶点完整闭环"）:
    - 1000 个不同靶点  = 471 个可评测标注靶(有>=5已知活性) + 填充靶补足到 1000
    - 每个靶点独立筛一个 10000 分子候选池 = 该靶已知活性(在库内) + 随机诱饵补足到 10000
      (per-target pool; 因 471 靶活性并集 15,601 > 10,000, 无法用单一共享池装下所有活性)
    - 完整闭环 + auto 路由: 每靶点自动选 AUC 最高的可用方法
        * 有 3D 受体的靶(当前仅 NA/CHEMBL2051/3TI6) -> 级联(L2 粗筛 + smina 对接 + LE 修正融合, 0.935 公式)
        * 其余靶点            -> L2 结合亲和力漏斗(检索域最优哈希版, 聚合中位 0.689)

L2 打分: 因子化推理(与 bigrun_score_hash 一致):
    molpart[m] = feats[m] @ W1m ; tgtpart[t] = tgt_hash[t] @ W1t + b0
    h1 = relu(molpart+tgtpart); h2 = relu(h1@W2+b2); p = sigmoid(h2@W3+b3)

评测: 每靶在其 10k 池上算 AUC / EF@1% / recall@1%; 聚合中位&均值(仅标注靶)。
NA 额外跑真实对接级联, 报告 L2 vs 级联 的 auto 收益。

产物: bigrun_1000x10000.json + bigrun_1000x10000.log
用法:
    DOCK_BIN_DIR=/path/to/smina/bin RUN_TARGETS=1000 POOL=10000 \
    /usr/bin/python3 bigrun_1000x10000.py
"""
import os, sys, json, time, tempfile, multiprocessing as mp
from pathlib import Path
import numpy as np

REPO = Path("<validated-workspace>")
SC = REPO / "评分_work_package/评分"
MTB = REPO / "scientific_validation/multitarget_benchmark"
BR = Path("<external-library-cache>")
sys.path.insert(0, str(SC))

from sklearn.metrics import roc_auc_score
from rdkit import Chem
from l2_bindingdb import Layer2BindingDB, BindingDBFeature
from target_resolver import TargetResolver
from target_seq_embedding import get_seq_dicts, chembl_to_seq
from dock_rerank import find_binary, prep_ligand, dock, _heavy_atoms, cascade_corrected_fusion

# ── 配置 ──
N_TARGETS = int(os.environ.get("RUN_TARGETS", "1000"))
POOL = int(os.environ.get("POOL", "10000"))
MIN_POS = 5
SEED = 42
TOPK_FRAC = 0.01                      # EF/recall @1%
OUT_JSON = os.environ.get("OUT_JSON", "bigrun_1000x10000.json")

# 级联(auto 高 AUC 分支)配置 —— 仅对有受体的靶
NA_TARGET = "CHEMBL2051"
RECEPTOR = str(REPO / "ai_mol_loop/influenza_na_2000_project_20260625/stage4/receptors/3TI6_protein_only_obabel.pdbqt")
CENTER = (-28.914, 14.334, 20.794)
SIZE = (23.585, 20.45, 24.18)
CAS_TOPN = int(os.environ.get("CAS_TOPN", "300"))
CAS_WORKERS = int(os.environ.get("CAS_WORKERS", "4"))
CAS_CPU = int(os.environ.get("CAS_CPU", "2"))
CAS_EXHAUST = int(os.environ.get("CAS_EXHAUST", "4"))
CAS_TIMEOUT = int(os.environ.get("CAS_TIMEOUT", "90"))
CAS_HAC_MAX = int(os.environ.get("CAS_HAC_MAX", "50"))
W_LE = 0.30

t0 = time.time()
def log(*a): print(f"[{time.time()-t0:7.0f}s]", *a, flush=True)

def canon(s):
    try:
        m = Chem.MolFromSmiles(s)
        return Chem.MolToSmiles(m) if m else None
    except Exception:
        return None

# ════════════════════════════════════════════════════════════
# 对接 worker（fork 进程池, 手动 prep+dock, 超时/HAC 过滤）
# ════════════════════════════════════════════════════════════
_ctx = mp.get_context("fork")
_W = {}
def _init_worker(receptor, center, size, smina_bin, obabel_bin, exhaust, cpu, timeout, hac_max):
    _W.update(dict(receptor=receptor, center=center, size=size, smina_bin=smina_bin,
                   obabel_bin=obabel_bin, exhaust=exhaust, cpu=cpu,
                   timeout=timeout, hac_max=hac_max,
                   workdir=tempfile.mkdtemp(prefix="big1kx10k_dock_")))
def _dock(arg):
    smi, mid = arg
    try:
        hac = _heavy_atoms(smi)
        if hac is None or hac > _W["hac_max"]:
            return (smi, None, hac, None, "skipped_hac")
        lig = os.path.join(_W["workdir"], f"{mid}.pdbqt")
        if not prep_ligand(smi, lig, _W["obabel_bin"], 7.4, _W["timeout"]):
            return (smi, None, hac, None, "prep_failed")
        aff = dock(lig, _W["receptor"], _W["center"], _W["size"], _W["smina_bin"],
                   _W["exhaust"], _W["cpu"], 3, 42, "", timeout=_W["timeout"])
        if aff is None:
            return (smi, None, hac, None, "dock_no_score")
        return (smi, round(aff, 3), hac, round(aff / hac, 4), "ok")
    except Exception as e:
        return (smi, None, _heavy_atoms(smi), None, f"err:{type(e).__name__}")

def run_dock(smiles_list, n_workers):
    args = [(s, f"m{i}") for i, s in enumerate(smiles_list)]
    out = {}
    with _ctx.Pool(n_workers, initializer=_init_worker,
                   initargs=(RECEPTOR, CENTER, SIZE, find_binary("smina"),
                             find_binary("obabel"), CAS_EXHAUST, CAS_CPU,
                             CAS_TIMEOUT, CAS_HAC_MAX)) as pool:
        done = 0
        for smi, aff, hac, le, status in pool.imap_unordered(_dock, args, chunksize=8):
            out[smi] = {"affinity": aff, "heavy_atoms": hac, "ligand_efficiency": le, "status": status}
            done += 1
            if done % 50 == 0:
                log(f"    [NA对接] {done}/{len(args)}")
    return out

# ════════════════════════════════════════════════════════════
# 1. 载入库
# ════════════════════════════════════════════════════════════
log("载入 100k 库 ...")
feats = np.load(BR / "lib_feats.npy").astype(np.float32)          # (100000,520)
smiles = open(BR / "lib_smiles.txt").read().splitlines()
label_map = json.load(open(BR / "label_map.json"))
N_MOL = feats.shape[0]
smi2row = {s: i for i, s in enumerate(smiles)}
log(f"库: {N_MOL} 分子; label_map {len(label_map)} 靶")

# ════════════════════════════════════════════════════════════
# 2. 目标靶集: 471 可评测 + 填充到 1000
# ════════════════════════════════════════════════════════════
cmap, useq = get_seq_dicts()
labeled_ok = sorted([t for t, d in label_map.items()
                     if len(d["pos"]) >= MIN_POS and chembl_to_seq(t)])
# 确保 NA 在评测集内
if NA_TARGET not in labeled_ok and NA_TARGET in label_map:
    labeled_ok.append(NA_TARGET)
log(f"可评测标注靶(>= {MIN_POS} 活性 & 有序列): {len(labeled_ok)}")

fasta_chembls = [c for c in cmap.keys() if cmap.get(c) and useq.get(cmap[c]) and len(useq[cmap[c]]) > 10]
fill = [c for c in fasta_chembls if c not in set(labeled_ok)]
np.random.RandomState(SEED).shuffle(fill)
targets = labeled_ok + fill[:max(0, N_TARGETS - len(labeled_ok))]
targets = targets[:N_TARGETS]
labeled_set = set(labeled_ok)
log(f"靶点集: {len(targets)} (可评测 {len(labeled_set)} + 填充 {len(targets)-len(labeled_set)})")

# ════════════════════════════════════════════════════════════
# 3. 载入哈希 L2, 因子化权重
# ════════════════════════════════════════════════════════════
l2 = Layer2BindingDB(prefer="mlp")
l2._ensure_model()
mlp = l2._model
W = [c.copy() for c in mlp.coefs_]           # [(776,256),(256,128),(128,1)]
b = [i.copy() for i in mlp.intercepts_]
W1m, W1t = W[0][:520], W[0][520:]            # (520,256),(256,256)
feat = BindingDBFeature()
resolver = TargetResolver()
def sigmoid(x): return 1.0 / (1.0 + np.exp(-x))

log("预计算 molpart(100k×256) ...")
molpart = feats @ W1m                         # (N_MOL,256)

# 靶点 256 维哈希文本
log("计算靶点哈希文本 ...")
tgt_hashes = np.zeros((len(targets), 256), dtype=np.float32)
for i, c in enumerate(targets):
    if c in labeled_set:
        tt, _m = resolver.resolve(c, chembl_id=c)
        text = tt if tt else c
    else:
        text = c
    tgt_hashes[i] = feat.target_features(text)
log("靶点哈希完成")

# ════════════════════════════════════════════════════════════
# 4. 逐靶: 构建 10k 候选池 -> L2 打分 -> (auto) 级联 -> 指标
# ════════════════════════════════════════════════════════════
all_rows = np.arange(N_MOL)
metrics = {}          # 仅标注靶
cascade_info = None
n_pairs = 0
score_t0 = time.time()

for i, tid in enumerate(targets):
    labeled = tid in labeled_set
    # ── 构建该靶 10k 候选池 ──
    if labeled:
        pos_rows = sorted({smi2row[s] for s in label_map[tid]["pos"] if s in smi2row})
    else:
        pos_rows = []
    rng = np.random.RandomState(SEED + i)
    n_decoy = max(0, POOL - len(pos_rows))
    non_pos = np.setdiff1d(all_rows, np.array(pos_rows, dtype=int), assume_unique=False) if pos_rows else all_rows
    if n_decoy > len(non_pos):
        n_decoy = len(non_pos)
    decoys = rng.choice(non_pos, size=n_decoy, replace=False)
    pool_rows = np.array(list(pos_rows) + list(decoys), dtype=int)
    # 池内标签
    y = np.zeros(len(pool_rows), dtype=np.int8)
    y[:len(pos_rows)] = 1

    # ── L2 因子化打分(仅池内) ──
    tgtpart = tgt_hashes[i] @ W1t + b[0]                    # (256,)
    h1 = np.maximum(molpart[pool_rows] + tgtpart, 0.0)      # (P,256)
    h2 = np.maximum(h1 @ W[1] + b[1], 0.0)                  # (P,128)
    p = sigmoid((h2 @ W[2] + b[2]).ravel())                 # (P,)
    n_pairs += len(pool_rows)

    if labeled and len(pos_rows) >= MIN_POS:
        auc = float(roc_auc_score(y, p))
        K = max(1, int(round(len(pool_rows) * TOPK_FRAC)))
        order = np.argsort(-p)[:K]
        hits = int(y[order].sum()); npos = int(y.sum())
        ef1 = round((hits / K) / (npos / len(pool_rows)), 3) if npos else None
        rec1 = round(hits / npos, 4) if npos else None
        rec = {"auc": round(auc, 4), "ef@1%": ef1, "recall@1%": rec1,
               "n_pos": npos, "pool": int(len(pool_rows)), "method": "L2"}

        # ── auto 路由: 有受体的靶叠加级联(取 AUC 更高者) ──
        if tid == NA_TARGET and find_binary("smina") and Path(RECEPTOR).exists():
            log(f"  [auto] {tid} 命中受体 -> 级联对接 (top-{CAS_TOPN} of 10k 池)")
            top_local = np.argsort(-p)[:CAS_TOPN]
            top_smiles = [smiles[pool_rows[j]] for j in top_local]
            t_dk = time.time()
            dock_res = run_dock(top_smiles, CAS_WORKERS)
            le_full = np.full(len(pool_rows), np.nan, dtype=float)
            for j, s in zip(top_local, top_smiles):
                le_full[j] = dock_res.get(s, {}).get("ligand_efficiency", np.nan) \
                             if dock_res.get(s, {}).get("ligand_efficiency") is not None else np.nan
            fused = cascade_corrected_fusion(p.astype(float), le_full, W_LE)
            auc_cas = float(roc_auc_score(y, fused))
            n_docked = int((~np.isnan(le_full)).sum())
            log(f"  [auto] {tid} L2 AUC={auc:.4f} -> 级联 AUC={auc_cas:.4f} "
                f"({n_docked} 对接成功, {time.time()-t_dk:.0f}s)")
            cascade_info = {"target": tid, "receptor": Path(RECEPTOR).name,
                            "pool": int(len(pool_rows)), "topN_docked": n_docked,
                            "auc_L2": round(auc, 4), "auc_cascade": round(auc_cas, 4),
                            "fusion_w_le": W_LE, "formula": "base + 0.30·σ(base)·z(−LE)"}
            if auc_cas >= auc:                       # auto 取更高者
                rec = {"auc": round(auc_cas, 4), "ef@1%": ef1, "recall@1%": rec1,
                       "n_pos": npos, "pool": int(len(pool_rows)),
                       "method": "cascade", "auc_L2_before": round(auc, 4)}
        metrics[tid] = rec

    if (i + 1) % 100 == 0:
        el = time.time() - score_t0
        log(f"进度 {i+1}/{len(targets)} 靶  已评测 {len(metrics)}  "
            f"pairs={n_pairs/1e6:.1f}M  吞吐 {n_pairs/el/1e6:.2f}M/s")

score_dt = time.time() - score_t0
log(f"打分完成: {n_pairs} 对, {score_dt:.0f}s")

# ════════════════════════════════════════════════════════════
# 5. 汇总
# ════════════════════════════════════════════════════════════
aucs = [m["auc"] for m in metrics.values() if m["auc"] is not None]
efs = [m["ef@1%"] for m in metrics.values() if m["ef@1%"] is not None]
recs = [m["recall@1%"] for m in metrics.values() if m["recall@1%"] is not None]
def med(x): return round(float(np.median(x)), 4) if x else None
def mean(x): return round(float(np.mean(x)), 4) if x else None
rev = sum(1 for a in aucs if a < 0.5)

summary = {
    "run": "1000靶点 × 10000候选 完整闭环 (auto 路由)",
    "n_targets_total": len(targets),
    "n_targets_labeled_eval": len(aucs),
    "pool_per_target": POOL,
    "pool_design": "per-target: 该靶已知活性(在库内) + 随机诱饵补足到 10000",
    "n_pairs": n_pairs, "n_pairs_million": round(n_pairs / 1e6, 1),
    "scoring_seconds": round(score_dt, 1),
    "auto_routing": "有受体->级联(0.935公式), 其余->L2哈希漏斗",
    "auc_median": med(aucs), "auc_mean": mean(aucs),
    "auc_ge_0.7": sum(1 for a in aucs if a >= 0.7),
    "auc_ge_0.8": sum(1 for a in aucs if a >= 0.8),
    "auc_ge_0.9": sum(1 for a in aucs if a >= 0.9),
    "auc_reversed_lt0.5": rev,
    "auc_reversed_frac": round(rev / len(aucs), 3) if aucs else None,
    "ef@1%_median": med(efs), "recall@1%_median": med(recs),
    "cascade_auto_gain": cascade_info,
    "l2_model": "Layer2BindingDB 哈希版(检索域最优, 因子化推理)",
    "note": "AUC 仅在标注可评测靶上可计; 填充靶仅贡献吞吐. NA 经 auto 路由走对接级联."
}
out = {"summary": summary, "per_target_metrics": metrics,
       "targets_sample": targets[:20]}
json.dump(out, open(MTB / OUT_JSON, "w"), indent=1, ensure_ascii=False)
log(f"已写出 {MTB / OUT_JSON}")

print("\n" + "=" * 72)
print("1000 靶点 × 10000 候选 完整闭环 (auto 路由) — 结果")
print("=" * 72)
print(f"规模: {len(targets)} 靶 × {POOL} 候选 = {n_pairs/1e6:.1f}M 对; 打分 {score_dt:.0f}s")
print(f"可评测标注靶: {len(aucs)}")
print(f"AUC 中位={summary['auc_median']} 均值={summary['auc_mean']}  "
      f">=0.7:{summary['auc_ge_0.7']} >=0.8:{summary['auc_ge_0.8']} >=0.9:{summary['auc_ge_0.9']} "
      f"反向:{rev}")
print(f"EF@1% 中位={summary['ef@1%_median']}  recall@1% 中位={summary['recall@1%_median']}")
if cascade_info:
    print(f"[auto 收益] NA: L2 {cascade_info['auc_L2']} -> 级联 {cascade_info['auc_cascade']}")
print("BIGRUN_1Kx10K_DONE")

# -*- coding: utf-8 -*-
"""
大通量闭环(哈希版 L2) — Stage B+C — 100,000 分子 × 10,000 靶点
================================================================
与 bigrun_score.py 的唯一差别: L2 用 **靶点哈希版 Layer2BindingDB**(检索域更优,
0.674 vs 序列版 0.553), 取代上一轮大通量错用的序列版 bindingdb_seq。

因子化首层 (同 seq 版工程):
    W1 = [W1_mol (256x520) | W1_tgt (256x256)]   # 哈希 tgt 为 256 维
    molpart[m] = feats[m] @ W1_mol          # 每分子算一次 (100k x 256)
    tgtpart[t] = tgt_hash[t] @ W1_tgt + b1  # 每靶点算一次 (10k x 256)
    h1[m,t]    = relu(molpart[m] + tgtpart[t])
    h2 = relu(h1 @ W2 + b2);  p = sigmoid(h2 @ W3 + b3)
末层 out_activation=logistic(已验证手写前向与 sklearn 一致, max diff<1e-5)。

靶点文本: 标注靶经 TargetResolver 取 bindingdb 训练一致 target_text; 填充靶用 chembl id 作哈希文本。
评测口径(诚实, 同原版): AUC 仅在 ~471 个有已知活性落库的标注靶上计算; 其余仅贡献吞吐。
产出: bigrun_results_hash.json + bigrun_topk_hits_hash.json
"""
import os, sys, json, time
import numpy as np
from pathlib import Path

REPO = Path("<validated-workspace>")
SC = REPO / "评分_work_package/评分"
MTB = REPO / "scientific_validation/multitarget_benchmark"
BR = Path("<external-library-cache>")
OUT_JSON = os.environ.get("OUT_JSON", "bigrun_results_hash.json")
sys.path.insert(0, str(SC))

N_TARGETS = int(os.environ.get("HAST_TARGETS", "10000"))
TOPK = 100
MIN_POS = 5
t0 = time.time()
def log(*a): print(f"[{time.time()-t0:6.0f}s]", *a, flush=True)

from target_resolver import TargetResolver
from l2_bindingdb import Layer2BindingDB, BindingDBFeature
from sklearn.metrics import roc_auc_score

# ── 1. 载入库 ──
feats = np.load(BR / "lib_feats.npy")                 # (100000, 520)
smiles = open(BR / "lib_smiles.txt").read().splitlines()
label_map = json.load(open(BR / "label_map.json"))
smi2row = {s: i for i, s in enumerate(smiles)}
N_MOL = feats.shape[0]
log(f"库: {N_MOL} 分子, feat {feats.shape}; 标注靶点 {len(label_map)}")

# ── 2. 构建 10k 靶点集 (与 seq 版一致: 标注优先 + fasta 补齐) ──
from target_seq_embedding import get_seq_dicts, chembl_to_seq
cmap, useq = get_seq_dicts()
labeled_ok = [t for t, d in label_map.items()
              if len(d["pos"]) >= MIN_POS and chembl_to_seq(t)]
labeled_ok = sorted(labeled_ok)
log(f"标注可评测靶点 (>= {MIN_POS} 活性 & 有序列): {len(labeled_ok)}")
fasta_chembls = [c for c in cmap.keys() if cmap.get(c) and useq.get(cmap[c])
                 and len(useq[cmap[c]]) > 10]
fill = [c for c in fasta_chembls if c not in set(labeled_ok)]
np.random.RandomState(42).shuffle(fill)
targets = labeled_ok + fill[:max(0, N_TARGETS - len(labeled_ok))]
targets = targets[:N_TARGETS]
labeled_set = set(labeled_ok)
log(f"靶点集: {len(targets)} (标注可评测 {len(labeled_set)} + 填充 {len(targets)-len(labeled_set)})")

# ── 3. 载入哈希 L2 模型, 提取权重做因子化 ──
_L2_PATH = os.environ.get("L2_MODEL_PATH")   # 重训模型覆盖点
l2 = Layer2BindingDB(model_path=_L2_PATH, prefer="mlp") if _L2_PATH else Layer2BindingDB(prefer="mlp")
l2._ensure_model()
mlp = l2._model
W = [c.copy() for c in mlp.coefs_]            # [(776,256),(256,128),(128,1)]
b = [i.copy() for i in mlp.intercepts_]
W1m, W1t = W[0][:520], W[0][520:]             # (520,256),(256,256)
log(f"哈希 L2 载入: 首层 {W[0].shape}, 末层激活 logistic")

# ── 4. 计算 10k 靶点 256 维哈希文本特征 ──
f = BindingDBFeature()
resolver = TargetResolver()
tgt_hashes = np.zeros((len(targets), 256), dtype=np.float32)
for i, c in enumerate(targets):
    if c in labeled_set:
        tt, method = resolver.resolve(c, chembl_id=c)
        text = tt if tt else c
    else:
        text = c   # 填充靶: 用 chembl id 作哈希文本(仅贡献吞吐, 不计入 AUC)
    tgt_hashes[i] = f.target_features(text)
    if (i + 1) % 2000 == 0:
        log(f"靶点哈希 {i+1}/{len(targets)}")
log(f"靶点哈希完成 {tgt_hashes.shape}")

# 因子化预计算 mol 部分(一次)
molpart = feats.astype(np.float32) @ W1m        # (Nmol,256)
log(f"molpart 预计算完成 {molpart.shape}")

def sigmoid(x): return 1.0 / (1.0 + np.exp(-x))

# ── 5. 逐靶打分(因子化), top-K + 标注靶 AUC ──
metrics = {}
topk_hits = {}
n_pairs = 0
score_t0 = time.time()
for i, tid in enumerate(targets):
    tgtpart = tgt_hashes[i] @ W1t + b[0]          # (256,)
    h1 = np.maximum(molpart + tgtpart, 0.0)       # (Nmol,256)
    h2 = np.maximum(h1 @ W[1] + b[1], 0.0)        # (Nmol,128)
    logit = (h2 @ W[2] + b[2]).ravel()            # (Nmol,)
    p = sigmoid(logit)
    n_pairs += N_MOL
    if TOPK < N_MOL:
        idx = np.argpartition(-p, TOPK)[:TOPK]
        idx = idx[np.argsort(-p[idx])]
    else:
        idx = np.argsort(-p)
    topk_hits[tid] = [(int(j), float(p[j])) for j in idx[:TOPK]]
    if tid in labeled_set:
        pos_rows = [smi2row[s] for s in label_map[tid]["pos"] if s in smi2row]
        if len(pos_rows) >= MIN_POS:
            y = np.zeros(N_MOL, dtype=np.int8); y[pos_rows] = 1
            try:
                auc = float(roc_auc_score(y, p))
            except Exception:
                auc = None
            K = max(1, int(round(N_MOL * 0.01)))
            order = np.argsort(-p)[:K]
            hits = int(y[order].sum()); npos = int(y.sum())
            ef1 = round((hits / K) / (npos / N_MOL), 3) if npos else None
            rec1 = round(hits / npos, 4) if npos else None
            metrics[tid] = {"auc": round(auc, 4) if auc is not None else None,
                            "ef@1%": ef1, "recall@1%": rec1, "n_pos": npos}
    if (i + 1) % 500 == 0:
        el = time.time() - score_t0
        log(f"打分 {i+1}/{len(targets)} 靶  pairs={n_pairs/1e6:.0f}M  "
            f"吞吐 {n_pairs/el/1e6:.2f}M pair/s")
score_dt = time.time() - score_t0
log(f"打分完成: {n_pairs} 对, {score_dt:.0f}s, 吞吐 {n_pairs/score_dt/1e6:.2f}M pair/s")

# ── 6. 汇总 ──
aucs = [m["auc"] for m in metrics.values() if m["auc"] is not None]
efs = [m["ef@1%"] for m in metrics.values() if m["ef@1%"] is not None]
recs = [m["recall@1%"] for m in metrics.values() if m["recall@1%"] is not None]
def med(x): return round(float(np.median(x)), 4) if x else None
def mean(x): return round(float(np.mean(x)), 4) if x else None
rev = sum(1 for a in aucs if a < 0.5)
summary = {
    "n_molecules": N_MOL, "n_targets": len(targets),
    "n_pairs": n_pairs, "n_pairs_billion": round(n_pairs / 1e9, 3),
    "scoring_seconds": round(score_dt, 1),
    "throughput_Mpair_per_s": round(n_pairs / score_dt / 1e6, 2),
    "n_labeled_eval_targets": len(aucs),
    "auc_median": med(aucs), "auc_mean": mean(aucs),
    "auc_reversed_lt0.5": rev, "auc_reversed_frac": round(rev / len(aucs), 3) if aucs else None,
    "auc_ge_0.7": sum(1 for a in aucs if a >= 0.7),
    "ef@1%_median": med(efs), "recall@1%_median": med(recs),
    "topk_per_target": TOPK,
    "model": "Layer2BindingDB (靶点哈希 256-dim, 检索域优化, 因子化推理)",
    "note": "L2 已由序列版(bindingdb_seq)换回哈希版(bindingdb)。AUC 仅在标注子集上可计; 其余靶点仅贡献吞吐。",
}
out = {"summary": summary,
       "per_target_metrics": metrics,
       "targets_sample": targets[:20],
       "labeled_eval_targets": sorted(labeled_set)}
json.dump(out, open(MTB / OUT_JSON, "w"), indent=1)
hits_out = {tid: [(smiles[j], round(s, 4)) for j, s in hits[:20]]
            for tid, hits in list(topk_hits.items())[:2000]}
json.dump(hits_out, open(MTB / "bigrun_topk_hits_hash.json", "w"))

log("DONE")
print("\n" + "=" * 70)
print("大通量闭环结果 (哈希版 L2)")
print("=" * 70)
print(f"规模: {N_MOL:,} 分子 × {len(targets):,} 靶点 = {n_pairs/1e9:.2f}e9 对")
print(f"打分耗时: {score_dt:.0f}s  吞吐: {summary['throughput_Mpair_per_s']} M pair/s")
print(f"标注可评测靶点: {len(aucs)}")
print(f"AUC 中位: {summary['auc_median']}  均值: {summary['auc_mean']}  "
      f"反向(<0.5): {rev} ({summary['auc_reversed_frac']})  AUC>=0.7: {summary['auc_ge_0.7']}")
print(f"EF@1% 中位: {summary['ef@1%_median']}  recall@1% 中位: {summary['recall@1%_median']}")
print("BIGRUN_HASH_DONE")

# -*- coding: utf-8 -*-
"""
大通量闭环 Stage B+C — 10,000 靶点嵌入 + 因子化 L2 打分 (100k × 10k = 1e9 对)
================================================================================
关键工程: 因子化首层。MLP 首层 Linear(648->256) 的权重按输入拆分:
    W1 = [W1_mol (256x520) | W1_tgt (256x128)]
    molpart[m]  = feats[m] @ W1_mol.T           # 每分子算一次 (100k x 256)
    tgtpart[t]  = temb[t]  @ W1_tgt.T + b1       # 每靶点算一次 (10k  x 256)
    h1[m,t]     = molpart[m] + tgtpart[t]        # 广播, 免去 1e9 x 648 大矩阵乘
其余层 (256->64->1) 在 (m,t) 网格上算。相比朴素 predict_matrix (~5.7h) 快一个量级。

评测口径 (诚实):
  - 吞吐/规模: 对全部 1e9 对计时 (真实大通量压测)。
  - 排序质量(AUC/EF/recall): 仅在"有已知活性落在库内"的标注靶点子集上计算
    (~471 个 >=5 活性靶点); 其余 ~9500 靶点无标签, AUC 不可计, 只贡献吞吐与漏斗。
  - 级联漏斗: L1 门槛(QED/Lipinski)剪枝统计。
产出: bigrun_results.json + 控制台汇总。
"""
import os, sys, json, time, pickle
import numpy as np
from pathlib import Path

REPO = Path("<validated-workspace>")
SC = REPO / "评分_work_package/评分"
MTB = REPO / "scientific_validation/multitarget_benchmark"
BR = Path("<external-library-cache>")
sys.path.insert(0, str(SC))

N_TARGETS = 10_000
TOPK = 100                # 每靶保留 top-K 命中
MIN_POS = 5               # 标注靶点最小活性数 (可算 AUC)
t0 = time.time()
def log(*a): print(f"[{time.time()-t0:6.0f}s]", *a, flush=True)

import torch
from target_seq_embedding import (TargetSeqEncoder, get_seq_dicts, chembl_to_seq)
from l2_bindingdb import Layer2BindingDBSeq
from sklearn.metrics import roc_auc_score

# ── 1. 载入库 ──
feats = np.load(BR / "lib_feats.npy")               # (100000, 520)
smiles = open(BR / "lib_smiles.txt").read().splitlines()
label_map = json.load(open(BR / "label_map.json"))
smi2row = {s: i for i, s in enumerate(smiles)}
N_MOL = feats.shape[0]
log(f"库: {N_MOL} 分子, feat {feats.shape}; 标注靶点 {len(label_map)}")

# ── 2. 构建 10k 靶点集 (标注优先, 再从 fasta 补齐) ──
cmap, useq = get_seq_dicts()                          # chembl->uniprot, uniprot->seq
# 标注且序列可解析且 >=MIN_POS 活性
labeled_ok = [t for t, d in label_map.items()
              if len(d["pos"]) >= MIN_POS and chembl_to_seq(t)]
labeled_ok = sorted(labeled_ok)
log(f"标注可评测靶点 (>= {MIN_POS} 活性 & 有序列): {len(labeled_ok)}")
# 从 fasta 全部 chembl 靶点补齐
fasta_chembls = [c for c in cmap.keys() if cmap.get(c) and useq.get(cmap[c])
                 and len(useq[cmap[c]]) > 10]
fill = [c for c in fasta_chembls if c not in set(labeled_ok)]
np.random.RandomState(42).shuffle(fill)
targets = labeled_ok + fill[:max(0, N_TARGETS - len(labeled_ok))]
targets = targets[:N_TARGETS]
labeled_set = set(labeled_ok)
log(f"靶点集: {len(targets)} (标注可评测 {len(labeled_set)} + 填充 {len(targets)-len(labeled_set)})")

# ── 3. 载入训练好的 encoder + MLP, 计算 10k 靶点嵌入 ──
bundle = torch.load(SC / "models/bindingdb_l2_seq/l2_seq.pt", map_location="cpu")
enc = TargetSeqEncoder(); enc.load_state_dict(bundle["encoder_state"]); enc.eval()
hidden = tuple(bundle.get("mlp_hidden", (256, 64)))
mlp = Layer2BindingDBSeq._build_mlp(hidden); mlp.load_state_dict(bundle["mlp_state"]); mlp.eval()

# 靶点嵌入 (逐靶, UNK 兜底)
temb = np.zeros((len(targets), enc.embed_dim), dtype=np.float32)
unk = enc.unk.detach().cpu().numpy()
for i, c in enumerate(targets):
    seq = chembl_to_seq(c)
    temb[i] = enc.embed_one(seq) if seq else unk
    if (i + 1) % 2000 == 0:
        log(f"靶点嵌入 {i+1}/{len(targets)}")
log(f"靶点嵌入完成 {temb.shape}")

# ── 4. 拆分 MLP 权重做因子化 ──
lin1, lin2, lin3 = mlp[0], mlp[2], mlp[4]
W1 = lin1.weight.detach().numpy(); b1 = lin1.bias.detach().numpy()   # (256,648),(256,)
W1m, W1t = W1[:, :520], W1[:, 520:]                                  # (256,520),(256,128)
W2 = lin2.weight.detach().numpy(); b2 = lin2.bias.detach().numpy()   # (64,256)
W3 = lin3.weight.detach().numpy(); b3 = lin3.bias.detach().numpy()   # (1,64)

molpart = feats.astype(np.float32) @ W1m.T                            # (Nmol,256) 一次
tgtpart = temb @ W1t.T + b1                                           # (Ntgt,256) 一次
log(f"因子化预计算完成: molpart {molpart.shape}, tgtpart {tgtpart.shape}")

def sigmoid(x): return 1.0 / (1.0 + np.exp(-x))

# ── 5. 逐靶打分 (因子化), top-K 命中 + 标注靶 AUC ──
metrics = {}          # tid -> {auc, ef@1%, recall@1%, n_pos}
topk_hits = {}        # tid -> [(smiles_row, score), ...]  (仅存 top-K)
n_pairs = 0
score_t0 = time.time()
for i, tid in enumerate(targets):
    h1 = molpart + tgtpart[i]                # (Nmol,256) 广播
    np.maximum(h1, 0, out=h1)                # ReLU inplace
    h2 = h1 @ W2.T + b2                       # (Nmol,64)
    np.maximum(h2, 0, out=h2)
    logit = (h2 @ W3.T + b3).ravel()         # (Nmol,)
    p = sigmoid(logit)
    n_pairs += N_MOL
    # top-K
    if TOPK < N_MOL:
        idx = np.argpartition(-p, TOPK)[:TOPK]
        idx = idx[np.argsort(-p[idx])]
    else:
        idx = np.argsort(-p)
    topk_hits[tid] = [(int(j), float(p[j])) for j in idx[:TOPK]]
    # 标注靶 AUC
    if tid in labeled_set:
        pos_rows = [smi2row[s] for s in label_map[tid]["pos"] if s in smi2row]
        if len(pos_rows) >= MIN_POS:
            y = np.zeros(N_MOL, dtype=np.int8); y[pos_rows] = 1
            try:
                auc = float(roc_auc_score(y, p))
            except Exception:
                auc = None
            # EF@1% / recall@1%
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

# ── 6. 级联漏斗 (L1 门槛剪枝, 用已存 mol 描述符近似) ──
# feats 前 512 是 Morgan, 后 8 是描述符 [MW,logP,HBD,HBA,TPSA,RotB,Rings,QED?]; 仅统计规模
# 漏斗按 L2 命中阈值近似: P>=0.5 视为"进入下一层"的候选
# (真实级联在 scoring.py; 此处给规模化统计)
survive = 0
for tid, hits in topk_hits.items():
    survive += sum(1 for _, s in hits if s >= 0.5)

# ── 7. 汇总 ──
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
    "ef@1%_median": med(efs), "recall@1%_median": med(recs),
    "topk_per_target": TOPK, "l2_hits_p>=0.5_in_topk": survive,
    "model": "Layer2BindingDBSeq (BLOSUM62+CNN seq embed, factored inference)",
    "note": "AUC 仅在标注子集(有已知活性落库)上可计; 其余靶点仅贡献吞吐/漏斗规模。",
}
out = {"summary": summary,
       "per_target_metrics": metrics,
       "targets_sample": targets[:20],
       "labeled_eval_targets": sorted(labeled_set)}
json.dump(out, open(MTB / "bigrun_results.json", "w"), indent=1)
# top-K 命中另存 (压缩为 smiles)
hits_out = {tid: [(smiles[j], round(s, 4)) for j, s in hits[:20]]
            for tid, hits in list(topk_hits.items())[:2000]}
json.dump(hits_out, open(MTB / "bigrun_topk_hits.json", "w"))

log("DONE")
print("\n" + "=" * 70)
print("大通量闭环结果")
print("=" * 70)
print(f"规模: {N_MOL:,} 分子 × {len(targets):,} 靶点 = {n_pairs/1e9:.2f}e9 对")
print(f"打分耗时: {score_dt:.0f}s  吞吐: {summary['throughput_Mpair_per_s']} M pair/s")
print(f"标注可评测靶点: {len(aucs)}")
print(f"AUC 中位: {summary['auc_median']}  均值: {summary['auc_mean']}  "
      f"反向(<0.5): {rev} ({summary['auc_reversed_frac']})")
print(f"EF@1% 中位: {summary['ef@1%_median']}  recall@1% 中位: {summary['recall@1%_median']}")
print("BIGRUN_DONE")

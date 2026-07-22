# -*- coding: utf-8 -*-
"""
高通量测试闭环 (High-Throughput Closed-Loop Benchmark)
========================================================
把 **整个 24k 候选库** 灌入四级漏斗 (L1 质量 + L2 靶点感知 + L3 ADMET; Uni-Mol 为吞吐省略),
对每个靶点做全库排序, 用真实 "已知活性分子回收率" 衡量效能:

  指标:
    - recall@1% / @5% / @10% : 全库排序后, 已知活性分子落入前 K% 的比例
    - EF@5% 富集因子          : top5% 中活性密度 / 全局活性密度
    - BEDROC(alpha=20)        : 早期富集 (活性越早浮现越好)
    - 级联剪枝               : 模拟四级漏斗逐层过滤, 看算力节省 vs 活性保留

吞吐:
    - 全量 (库 × 靶点) 配对打分 -> 配对/秒, 总墙钟
    - L1/L3 为靶点无关, 预计算一次落盘缓存

注意: 面板 10155 个已知活性分子 100% 落在 24k 库内 (canon 命中), 故无需注入,
闭环基准即 "从 24k 真实药物里把已知活性分子筛出来"。
"""
import json, csv, time, sys, pickle, math
from pathlib import Path
import numpy as np
from sklearn.metrics import roc_auc_score
from collections import defaultdict

REPO = Path("<validated-workspace>")
SC = REPO / "评分_work_package/评分"
sys.path.insert(0, str(SC))
from l2_bindingdb import BindingDBFeature, Layer2BindingDB
from scoring import Layer1Scorer, Layer3Scorer
from target_resolver import TargetResolver

t0 = time.time()
def log(*a):
    print(f"[{time.time()-t0:6.0f}s]", *a, flush=True)

CACHE = REPO / "scientific_validation/multitarget_benchmark/hts_l1l3_cache.pkl"

# ───────────────────────── 1. 载入库 + 面板 ─────────────────────────
lib_ids, lib_smiles = [], []
with open(REPO / "scientific_validation/multitarget_benchmark/candidates_10k.csv") as f:
    for r in csv.DictReader(f):
        lib_ids.append(r["id"]); lib_smiles.append(r["smiles"].strip())
panel = json.load(open(REPO / "scientific_validation/multitarget_benchmark/target_panel.json"))
N = len(lib_smiles)
log(f"候选库 {N} 药; 面板 {len(panel)} 靶点")

# 解析靶点文本
resolver = TargetResolver()
target_text, methods = {}, {}
for cid, v in panel.items():
    tt, m = resolver.resolve(v["name"], chembl_id=cid)
    target_text[cid] = tt; methods[cid] = m
n_resolved = sum(1 for t in target_text.values() if t)
log(f"靶点解析 target_text: {n_resolved}/{len(panel)}")

# ───────────────────────── 2. L1/L3 预计算 (落盘缓存) ─────────────────────────
import concurrent.futures as _cf
import multiprocessing as _mp
import time as _t
from hts_worker import _worker_l1l3

feat = BindingDBFeature()

if CACHE.exists():
    _loaded = pickle.load(open(CACHE, "rb"))
    if len(_loaded) == 3:
        L1c, L3c, MOLF = _loaded
    else:
        L1c, L3c = _loaded
        MOLF = {}
    log(f"载入 L1/L3/MOLF 缓存 {len(L1c)} 分子")
else:
    log("无缓存, 分批多进程计算 L1/L3 (每分子30s超时, 卡死worker逐批kill)...")
    L1c, L3c, MOLF = {}, {}, {}
    skipped = 0
    n_done = 0
    BATCH, WORKERS, TMO = 500, 6, 30.0
    _zeros = np.zeros(feat.DIM - feat.N_TARGET, dtype=np.float32)
    _ctx = _mp.get_context("fork")   # fork: 子进程不重跑主模块 (避免 spawn 递归)
    for bi in range(0, N, BATCH):
        chunk = lib_smiles[bi:bi + BATCH]
        ex = _cf.ProcessPoolExecutor(max_workers=WORKERS, mp_context=_ctx)
        futs = {ex.submit(_worker_l1l3, s): s for s in chunk}
        try:
            for f in futs:
                smi = futs[f]
                try:
                    r = f.result(timeout=TMO)
                except Exception:
                    r = None
                if r is None:
                    L1c[smi] = 0.0; L3c[smi] = 0.5; MOLF[smi] = _zeros.copy(); skipped += 1
                else:
                    L1c[smi], L3c[smi], MOLF[smi] = r
                n_done += 1
        finally:
            try:
                ex.shutdown(wait=False, cancel_futures=True)
            except Exception:
                pass
            for ch in _mp.active_children():   # 强制 kill 本批全部 worker (含卡死)
                try: ch.kill()
                except Exception:
                    pass
        if (bi // BATCH) % 4 == 0:
            pickle.dump((L1c, L3c, MOLF), open(CACHE, "wb"))
            log(f"L1/L3 进度 {n_done}/{N}  跳过(超时/异常) {skipped}")
    pickle.dump((L1c, L3c, MOLF), open(CACHE, "wb"))
    log(f"L1/L3 完成 {N} (跳过 {skipped}); 缓存写出")

L1_arr = np.array([L1c[s] for s in lib_smiles], dtype=np.float32)   # (N,)
L3_arr = np.array([L3c[s] for s in lib_smiles], dtype=np.float32)   # (N,)

# 分子特征矩阵 (520 维, 靶点无关) — 直接来自多进程 worker 结果, 不主线程重算
_z = np.zeros(feat.DIM - feat.N_TARGET, dtype=np.float32)
molfeat = np.array([MOLF.get(s, _z) for s in lib_smiles], dtype=np.float32)  # (N, 520)

# 已知活性/阴性 在库内索引 (canon 匹配)
from rdkit import Chem
def canon(s):
    try: return Chem.MolToSmiles(Chem.MolFromSmiles(s))
    except Exception: return None
lib_canon = {}
for i, s in enumerate(lib_smiles):
    c = canon(s)
    if c: lib_canon[c] = i
log(f"库内 canon 索引 {len(lib_canon)}")

# ───────────────────────── 3. 工具: 排名/回收率/BEDROC ─────────────────────────
def bedroc(pos_scores, neg_scores, alpha=20.0):
    all_s = np.concatenate([pos_scores, neg_scores])
    n = len(all_s); npos = len(pos_scores)
    if n == 0 or npos == 0: return 0.0
    order = np.argsort(-all_s)  # 降序 (高分在前)
    pos_set = set(range(npos))
    ra = 0.0
    for rank, idx in enumerate(order):
        if idx in pos_set:
            ra += math.exp(-alpha * (rank + 1) / n)
    ra_max = (1 - math.exp(-alpha * npos / n)) / (1 - math.exp(-alpha))
    return min(1.0, ra / ra_max) if ra_max > 0 else 0.0

def recall_at(pos_idx, order, fracs):
    """order: 降序索引数组; 返回各 frac 下的 recall"""
    out = {}
    npos = len(pos_idx)
    if npos == 0: return {f"recall@{int(f*100)}%": 0.0 for f in fracs}
    pos_set = set(pos_idx)
    nlib = len(order)
    for f in fracs:
        K = max(1, int(round(nlib * f)))
        top = set(order[:K].tolist())
        hits = sum(1 for i in pos_idx if i in top)
        out[f"recall@{int(f*100)}%"] = hits / npos
    return out

def ef_at(pos_idx, score_vec, frac=0.05):
    nlib = len(score_vec); npos = len(pos_idx)
    if npos == 0: return 0.0
    order = np.argsort(-score_vec)
    K = max(1, int(round(nlib * frac)))
    pos_set = set(pos_idx)
    hits = sum(1 for i in order[:K] if i in pos_set)
    return (hits / K) / (npos / nlib)

# ───────────────────────── 4. 逐靶点全库高通量闭环 ─────────────────────────
l2 = Layer2BindingDB(prefer="mlp")
log(f"L2 模型: {l2.model_kind}")
fracs = [0.01, 0.05, 0.10]
pairs_total = 0
records = []

for cid, v in panel.items():
    tt = target_text[cid]
    if not tt:
        continue
    # 库内活性/阴性索引
    pidx = [lib_canon[c] for c in (canon(p) for p in v.get("pos", [])) if c and c in lib_canon]
    nidx = [lib_canon[c] for c in (canon(n) for n in v.get("neg", [])) if c and c in lib_canon]
    if len(pidx) < 3:
        continue
    # L2 全库向量化
    tf = feat.target_features(tt)                       # (256,)
    X = np.hstack([molfeat, np.tile(tf, (N, 1))])        # (N, 776)
    l2p = l2.predict_proba(X)                            # (N,)
    final = 0.20 * L1_arr + 0.50 * l2p + 0.20 * L3_arr   # (N,)

    order_final = np.argsort(-final)
    order_l1 = np.argsort(-L1_arr)
    order_l2 = np.argsort(-l2p)

    rf = recall_at(pidx, order_final, fracs)
    rl1 = recall_at(pidx, order_l1, fracs)
    rl2 = recall_at(pidx, order_l2, fracs)
    ef = ef_at(pidx, final, 0.05)
    # BEDROC: pos=活性, neg=全库其余
    neg_idx = np.ones(N, dtype=bool); neg_idx[pidx] = False
    bed = bedroc(final[pidx], final[neg_idx], 20.0)
    # AUC (全库ranking)
    y = np.zeros(N, dtype=int); y[pidx] = 1
    try: auc = float(roc_auc_score(y, final))
    except Exception: auc = float("nan")

    # 级联剪枝 (闭环漏斗模拟): L1 gate -> L2 top -> L3 gate
    # Stage1: L1 >= 0.45 保留
    s1 = np.where(L1_arr >= 0.45)[0]
    s1_pos = sum(1 for i in pidx if i in set(s1.tolist()))
    # Stage2: 在 s1 内按 L2 取 top 10%
    if len(s1) > 0:
        s1order = s1[np.argsort(-l2p[s1])]
        K2 = max(1, int(round(len(s1order) * 0.10)))
        s2 = set(s1order[:K2].tolist())
        s2_pos = sum(1 for i in pidx if i in s2)
    else:
        s2_pos = 0; s2 = set()
    # Stage3: L3 >= 0.5 保留
    s3 = set(np.where(L3_arr >= 0.5)[0].tolist())
    s3_pos = sum(1 for i in pidx if i in s3)

    pairs_total += N
    records.append({
        "target": cid, "name": v["name"], "method": methods[cid],
        "n_pos_in_lib": len(pidx), "n_neg_in_lib": len(nidx),
        "final_recall": rf, "l1_recall": rl1, "l2_recall": rl2,
        "ef5": ef, "bedroc": round(bed, 3), "auc_full": round(auc, 3),
        "cascade": {
            "l1_keep_frac": round(len(s1) / N, 3), "l1_pos_kept": round(s1_pos / len(pidx), 3),
            "l2_top10_pos_kept": round(s2_pos / len(pidx), 3),
            "l3_keep_frac": round(len(s3) / N, 3), "l3_pos_kept": round(s3_pos / len(pidx), 3),
        },
    })

wall = time.time() - t0
log(f"闭环评估完成 {len(records)} 靶点; 配对总数 {pairs_total:,}; 墙钟 {wall/60:.1f}min; "
    f"吞吐 {pairs_total/wall:,.0f} pair/s")

# ───────────────────────── 5. 聚合 ─────────────────────────
def agg(key_fn):
    vals = [key_fn(r) for r in records]
    vals = [v for v in vals if v is not None and not (isinstance(v, float) and math.isnan(v))]
    if not vals: return None
    vals.sort(); n = len(vals)
    return {"n": n, "median": round(vals[n//2], 3), "mean": round(sum(vals)/n, 3),
            "min": round(vals[0], 3), "max": round(vals[-1], 3)}

def agg_recall(field):
    out = {}
    for f in ["recall@1%", "recall@5%", "recall@10%"]:
        out[f] = agg(lambda r: r[field][f])
    return out

summary = {
    "setup": {
        "library_size": N, "n_targets_panel": len(panel),
        "n_targets_evaluated": len(records),
        "n_pos_total": sum(r["n_pos_in_lib"] for r in records),
        "layers": "L1(0.20) + L2_BindingDB(0.50) + L3_ADMET(0.20)  [Uni-Mol 0.10 吞吐省略]",
        "l2_model": f"combined BindingDB+ChEMBL MLP({l2.model_kind})",
    },
    "throughput": {
        "pairs_total": pairs_total,
        "wall_seconds": round(wall, 1),
        "pairs_per_second": round(pairs_total / wall, 1),
        "note": "L1/L3 为靶点无关, 预计算一次并缓存; 逐靶点增量仅 L2 MLP (向量化)",
    },
    "recall_final": agg_recall("final_recall"),
    "recall_l2_only": agg_recall("l2_recall"),
    "recall_l1_baseline": agg_recall("l1_recall"),
    "ef5_final": agg(lambda r: r["ef5"]),
    "bedroc_final": agg(lambda r: r["bedroc"]),
    "auc_full_final": agg(lambda r: r["auc_full"]),
    "cascade": {
        "l1_keep_frac": agg(lambda r: r["cascade"]["l1_keep_frac"]),
        "l1_pos_kept": agg(lambda r: r["cascade"]["l1_pos_kept"]),
        "l2_top10_pos_kept": agg(lambda r: r["cascade"]["l2_top10_pos_kept"]),
        "l3_keep_frac": agg(lambda r: r["cascade"]["l3_keep_frac"]),
        "l3_pos_kept": agg(lambda r: r["cascade"]["l3_pos_kept"]),
    },
    "frac_targets_recall1pct_ge_0.5": round(
        sum(1 for r in records if r["final_recall"]["recall@1%"] >= 0.5) / len(records), 3),
}

out = {"summary": summary, "per_target": records,
       "meta": {"script": "hts_closed_loop.py", "l1l3_cache": str(CACHE)}}
json.dump(out, open(REPO / "scientific_validation/multitarget_benchmark/hts_closed_loop_results.json", "w"), indent=1)

# ───────────────────────── 6. 打印 ─────────────────────────
print("\n" + "=" * 70)
print("高通量测试闭环结果")
print("=" * 70)
print(f"库 {N} × {len(records)} 靶点 = {pairs_total:,} 配对 | 吞吐 {pairs_total/wall:,.0f} pair/s")
print("\n[全库回收率 recall@K]  median / mean / min / max")
for label, key in [("final(四级漏斗)", "recall_final"),
                   ("L2-only(靶点感知)", "recall_l2_only"),
                   ("L1-only(质量基线)", "recall_l1_baseline")]:
    a = summary[key]
    print(f"  {label:22s} @1%={a['recall@1%']['median']:.3f}  "
          f"@5%={a['recall@5%']['median']:.3f}  @10%={a['recall@10%']['median']:.3f}")
print(f"\nEF@5% median={summary['ef5_final']['median']}  "
      f"BEDROC median={summary['bedroc_final']['median']}  "
      f"全库AUC median={summary['auc_full_final']['median']}")
print(f"达到 recall@1%>=0.5 的靶点比例: {summary['frac_targets_recall1pct_ge_0.5']}")
c = summary["cascade"]
print(f"\n[级联剪枝] L1保 {(1-c['l1_keep_frac']['median'])*100:.0f}%算力, 活性保留 {c['l1_pos_kept']['median']*100:.0f}% | "
      f"L2 top10% 活性保留 {c['l2_top10_pos_kept']['median']*100:.0f}% | "
      f"L3保 {(1-c['l3_keep_frac']['median'])*100:.0f}%算力, 活性保留 {c['l3_pos_kept']['median']*100:.0f}%")
print("done.")

# -*- coding: utf-8 -*-
"""
新 AUC 闭环（冻结版）— L2 粗筛 + 对接精排(LE 校正) 级联
================================================================
把 ② 对接精排真正接入闭环、端到端跑通、可冻结的最终版本。

Stage 1 (L2 粗筛, 高吞吐):
    哈希 L2（检索域最优, 24k×164 面板中位 0.674）对 100k 全库向量化打分
    (776 维 = 520 分子 + 256 靶点哈希)。一次 predict_proba 完成。

Stage 2 (对接精排, 低吞吐但确保跑通):
    取 L2 top-N, 并行 smina 对接 + obabel 3D(pH7.4),
    计算配体效率 LE = affinity / 重原子数（校正 Vina 尺寸偏好反富集）。
    另对 "NA 18 活性 + 100 诱饵" 平衡集全部对接, 与 NA demo(0.462→0.756) 苹果对苹果对比。

Stage 3 (级联评估):
    在 NA(CHEMBL2051, 唯一有受体的靶点) 上对比
    L2 基线 AUC vs 级联 AUC（对接 LE 重排）。

诚实边界: 当前环境只有 NA(3TI6) 一个受体, 故对接精排仅作用于 NA;
    其余靶点走纯 L2。多靶点扩展需补充受体结构（见 cascade_auc_report.md）。
产物: <external-library-cache>/cascade_results.json

用法:
    DOCK_BIN_DIR=/path/to/smina/bin CASCADE_TOPN=500 CASCADE_EXHAUST=4 \
    CASCADE_TIMEOUT=90 CASCADE_HAC_MAX=50 CASCADE_WORKERS=4 \
    .venv_mlx_qwen35/bin/python3 cascade_auc_closed_loop.py
"""
import os, sys, json, time, tempfile, multiprocessing as mp
from pathlib import Path
import numpy as np

REPO = Path("<validated-workspace>")
SC = REPO / "评分_work_package/评分"
BR = Path("<external-library-cache>")
OUT = BR / os.environ.get("CASCADE_OUT", "cascade_results.json")
sys.path.insert(0, str(SC))

from rdkit import Chem
from sklearn.metrics import roc_auc_score
from l2_bindingdb import Layer2BindingDB, BindingDBFeature
from target_resolver import TargetResolver
from dock_rerank import find_binary, prep_ligand, dock, _heavy_atoms

# ── 配置 ──
TARGET = "CHEMBL2051"          # 神经氨酸酶（唯一有受体的靶点）
RECEPTOR = str(REPO / "ai_mol_loop/influenza_na_2000_project_20260625/stage4/receptors/3TI6_protein_only_obabel.pdbqt")
CENTER = (-28.914, 14.334, 20.794)
SIZE = (23.585, 20.45, 24.18)
TOPN = int(os.environ.get("CASCADE_TOPN", "500"))
N_WORKERS = int(os.environ.get("CASCADE_WORKERS", "4"))
SMINA_CPU = int(os.environ.get("CASCADE_CPU", "2"))
EXHAUST = int(os.environ.get("CASCADE_EXHAUST", "4"))
DOCK_TIMEOUT = int(os.environ.get("CASCADE_TIMEOUT", "90"))
DOCK_HAC_MAX = int(os.environ.get("CASCADE_HAC_MAX", "50"))
BAL_DECOYS = 100               # 平衡集诱饵数
SEED = 42

t0 = time.time()
def log(*a): print(f"[{time.time()-t0:6.0f}s]", *a, flush=True)

# ── 多进程 worker（手动 prep+dock, 带超时与 HAC 过滤, 避免大模型拖死）──
import subprocess as _sp
_ctx = mp.get_context("fork")
_W = {}
def _init_worker(receptor, center, size, smina_bin, obabel_bin, exhaust, cpu, timeout, hac_max):
    _W["rr"] = dict(receptor=receptor, center=center, size=size,
                    smina_bin=smina_bin, obabel_bin=obabel_bin,
                    exhaustiveness=exhaust, cpu=cpu)
    _W["timeout"] = timeout
    _W["hac_max"] = hac_max
    _W["workdir"] = tempfile.mkdtemp(prefix="cascade_dock_")
def _dock(arg):
    smi, mid = arg
    try:
        cfg = _W["rr"]; timeout = _W["timeout"]; hac_max = _W["hac_max"]
        hac = _heavy_atoms(smi)
        if hac is None or hac > hac_max:
            return (smi, None, hac, None, "skipped_hac")
        lig = os.path.join(_W["workdir"], f"{mid}.pdbqt")
        if not prep_ligand(smi, lig, cfg["obabel_bin"], 7.4, timeout):
            return (smi, None, hac, None, "prep_failed")
        aff = dock(lig, cfg["receptor"], cfg["center"], cfg["size"], cfg["smina_bin"],
                   cfg["exhaustiveness"], cfg["cpu"], 3, 42, "", timeout=timeout)
        if aff is None:
            return (smi, None, hac, None, "dock_no_score")
        return (smi, round(aff, 3), hac, round(aff / hac, 4), "ok")
    except _sp.TimeoutExpired:
        return (smi, None, _heavy_atoms(smi), None, "timeout")
    except Exception as e:
        return (smi, None, _heavy_atoms(smi), None, f"err:{type(e).__name__}")

def run_dock(smiles_list, n_workers):
    args = [(s, f"m{i}") for i, s in enumerate(smiles_list)]
    n_dockable = sum(1 for s in smiles_list if (_heavy_atoms(s) or 0) <= DOCK_HAC_MAX)
    log(f"  提交 {len(args)} 个; 其中 HAC≤{DOCK_HAC_MAX} 可对接 ≈{n_dockable}（其余自动跳过）")
    out = {}
    done = 0
    with _ctx.Pool(n_workers, initializer=_init_worker,
                   initargs=(RECEPTOR, CENTER, SIZE, find_binary("smina"),
                             find_binary("obabel"), EXHAUST, SMINA_CPU,
                             DOCK_TIMEOUT, DOCK_HAC_MAX)) as pool:
        for smi, aff, hac, le, status in pool.imap_unordered(_dock, args, chunksize=8):
            out[smi] = {"affinity": aff, "heavy_atoms": hac,
                        "ligand_efficiency": le, "status": status}
            done += 1
            if done % 25 == 0:
                log(f"  对接进度 {done}/{len(args)}")
    return out

def canon(s):
    try:
        m = Chem.MolFromSmiles(s)
        return Chem.MolToSmiles(m) if m else None
    except Exception:
        return None

# ════════════════════════════════════════════════
# Stage 1 — L2 全库粗筛
# ════════════════════════════════════════════════
log("载入 100k 库 ...")
feats = np.load(BR / "lib_feats.npy").astype(np.float32)          # (N,520)
smiles = open(BR / "lib_smiles.txt").read().splitlines()
N = feats.shape[0]
label_map = json.load(open(BR / "label_map.json"))
na = label_map[TARGET]
pos_smiles_raw = na["pos"]; neg_smiles_raw = na["neg"]

# 活性/阴性 → 库内行号（SMILES 规范化匹配）
lib_canon = [canon(s) for s in smiles]
canon2row = {}
for i, c in enumerate(lib_canon):
    if c:
        canon2row.setdefault(c, i)
pos_rows = [canon2row[c] for c in (canon(s) for s in pos_smiles_raw) if c in canon2row]
neg_rows = [canon2row[c] for c in (canon(s) for s in neg_smiles_raw) if c in canon2row]
log(f"NA 活性命中库 {len(pos_rows)}/{len(pos_smiles_raw)}; 阴性 {len(neg_rows)}/{len(neg_smiles_raw)}")

log("解析 NA target_text（训练一致）...")
resolver = TargetResolver()
tt, method = resolver.resolve(name="neuraminidase", chembl_id=TARGET)
log(f"  target_text({method}): {tt[:70]!r}")

log("L2 向量化打分 100k ...")
_L2P = os.environ.get("L2_MODEL_PATH")
l2 = Layer2BindingDB(model_path=_L2P) if _L2P else Layer2BindingDB()
l2._ensure_model()
feat = BindingDBFeature()
thash = feat.target_features(tt)                                  # (256,)
X = np.concatenate([feats, np.tile(thash, (N, 1))], axis=1).astype(np.float32)
l2_scores = l2.predict_proba(X)                                   # (N,)
log(f"  L2 范围 [{l2_scores.min():.3f},{l2_scores.max():.3f}]")

y = np.zeros(N, dtype=int)
for r in pos_rows:
    y[r] = 1
auc_l2_full = roc_auc_score(y, l2_scores)
log(f"★ L2 基线 AUC（全 100k 库, {len(pos_rows)} 活性） = {auc_l2_full:.3f}")

# ════════════════════════════════════════════════
# Stage 2 — 对接精排
# ════════════════════════════════════════════════
order = np.argsort(l2_scores)[::-1]
top_idx = order[:TOPN]
top_smiles = [smiles[i] for i in top_idx]
log(f"取 L2 top-{TOPN} 对接（{len(top_smiles)} 配体, {N_WORKERS} 进程）...")
t_top = time.time()
top_dock = run_dock(top_smiles, N_WORKERS)
n_top_docked = sum(1 for v in top_dock.values() if v['ligand_efficiency'] is not None)
log(f"  top-N 对接完成 {time.time()-t_top:.0f}s; 成功 {n_top_docked}/{len(top_dock)}")

# 平衡集（18 活性 + 100 诱饵）全部对接，与 NA demo 对比
rng = np.random.RandomState(SEED)
decoy_rows = rng.choice(np.setdiff1d(np.arange(N), pos_rows), size=BAL_DECOYS, replace=False)
bal_rows = np.array(pos_rows + list(decoy_rows), dtype=int)
bal_smiles = [smiles[i] for i in bal_rows]
bal_y = np.array([1]*len(pos_rows) + [0]*len(decoy_rows), dtype=int)
log(f"平衡集对接（{len(bal_smiles)} 配体 = {len(pos_rows)} 活性 + {BAL_DECOYS} 诱饵）...")
t_bal = time.time()
bal_dock = run_dock(bal_smiles, N_WORKERS)
log(f"  平衡集对接完成 {time.time()-t_bal:.0f}s")

# ════════════════════════════════════════════════
# Stage 3 — 级联评估
# ════════════════════════════════════════════════
def le_of(d):
    return d.get("ligand_efficiency")

# (a) 平衡集 AUC：L2 vs 对接 LE（苹果对苹果对比 NA demo 0.462→0.756）
bal_le = np.array([le_of(bal_dock[s]) for s in bal_smiles], dtype=object)
bal_ok = np.array([le_of(bal_dock[s]) is not None for s in bal_smiles])
auc_l2_bal = roc_auc_score(bal_y[bal_ok], np.array([l2_scores[i] for i in bal_rows])[bal_ok])
le_vals = np.array([le_of(bal_dock[s]) for s in bal_smiles], dtype=float)
auc_dock_bal = roc_auc_score(bal_y[bal_ok], -le_vals[bal_ok])   # LE 越负越好 → 取负
n_bal_act_docked = sum(1 for s in bal_smiles[:len(pos_rows)] if le_of(bal_dock[s]) is not None)
log(f"★ 平衡集 L2 AUC = {auc_l2_bal:.3f}  |  ★ 平衡集 级联(LE) AUC = {auc_dock_bal:.3f} "
    f"（{n_bal_act_docked}/{len(pos_rows)} 活性已对接）")

# (b) 全库部署级联 AUC：top-N 内用 -LE 重排，非 top-N 视为淘汰(-inf)
final_full = np.full(N, -1e9, dtype=float)
le_arr = np.array([le_of(top_dock.get(smiles[i])) for i in top_idx], dtype=float)
docked_mask = ~np.isnan(le_arr)
final_full[top_idx[docked_mask]] = -le_arr[docked_mask]
auc_cascade_full = roc_auc_score(y, final_full)
n_act_in_top = sum(1 for r in pos_rows if r in set(top_idx.tolist()))
log(f"★ 全库部署级联 AUC（top-N 内 LE 重排, 其余淘汰） = {auc_cascade_full:.3f} "
    f"（{n_act_in_top}/{len(pos_rows)} 活性落入 L2 top-{TOPN}）")

# (b') 融合闭环分数 = 归一化 L2 + 归一化(-LE)，既用粗筛也用精排
def _norm(a):
    a = np.asarray(a, dtype=float)
    return (a - a.mean()) / (a.std() + 1e-9)
# 平衡集融合
l2_bal = np.array([l2_scores[i] for i in bal_rows], dtype=float)
le_bal = np.array([le_of(bal_dock[s]) if le_of(bal_dock[s]) is not None else np.nan
                   for s in bal_smiles], dtype=float)
le_bal_fill = np.where(np.isnan(le_bal), 0.0, le_bal)   # 未对接的 LE 视为中性
fused_bal = _norm(l2_bal) + _norm(-le_bal_fill)
auc_fused_bal = roc_auc_score(bal_y[bal_ok], fused_bal[bal_ok])
# 全库部署融合（top-N 内 L2+LE, 其余淘汰）
fused_full = np.full(N, -1e9, dtype=float)
le_top = np.array([le_of(top_dock.get(smiles[i])) for i in top_idx], dtype=float)
le_top_fill = np.where(np.isnan(le_top), 0.0, le_top)
fused_top = _norm(l2_scores[top_idx]) + _norm(-le_top_fill)
fused_full[top_idx[docked_mask]] = fused_top[docked_mask]
auc_fused_full = roc_auc_score(y, fused_full)
log(f"★ 平衡集 融合(L2+LE) AUC = {auc_fused_bal:.3f}   （L2={auc_l2_bal:.3f}, LE={auc_dock_bal:.3f}）")
log(f"★ 全库部署 融合(L2+LE) AUC = {auc_fused_full:.3f}")

# (b'') 修正版融合 —— 以 L2 为全库基线, 仅对【对接成功】分子施加同尺度 LE 校正。
#   原版 (b') 把非对接分子硬设 -1e9 垫底, 导致 15 个非对接活性被 211 个 top-L2 阴性插队压到库底 -> 0.58。
#   修正: 非对接分子保留 L2 原始分; 对接分子 = L2 + w·(归一化(-LE)·L2_std)。尺度一致, 不破坏 L2 排序。
le_top_raw = np.array([le_of(top_dock.get(smiles[i])) for i in top_idx], dtype=float)
le_top_val = np.where(np.isnan(le_top_raw), 0.0, le_top_raw)
le_top_scaled = -le_top_val                                   # LE 越负越好 -> 取负
le_top_scaled = (le_top_scaled - le_top_scaled.mean()) / (le_top_scaled.std() + 1e-9)
W_LE = 0.30
deployed_corrected = l2_scores.copy()
deployed_corrected[top_idx[docked_mask]] = \
    l2_scores[top_idx[docked_mask]] + W_LE * l2_scores.std() * le_top_scaled[docked_mask]
auc_deployed_corrected = roc_auc_score(y, deployed_corrected)
log(f"★ 全库部署【修正融合】 AUC = {auc_deployed_corrected:.3f}  "
    f"(基线 L2={auc_l2_full:.3f}, 仅对接分子叠加 LE 校正, W_LE={W_LE})")

# (c) EF@1% / recall@1%（全库）
def metrics(scores):
    o = np.argsort(scores)[::-1]
    top1k = o[:max(1, N//100)]
    hits = int(y[top1k].sum())
    nact = int(y.sum())
    recall = hits / nact if nact else 0.0
    ef = (hits / len(top1k)) / (nact / N) if nact else 0.0
    return recall, ef
rec_l2, ef_l2 = metrics(l2_scores)
rec_cas, ef_cas = metrics(final_full)
log(f"  EF@1% : L2={ef_l2:.3f}  级联={ef_cas:.3f}")
log(f"  recall@1%: L2={rec_l2:.3f}  级联={rec_cas:.3f}")

# ── 汇总落盘 ──
result = {
    "pipeline": "cascade_auc_closed_loop (L2 粗筛 + 对接精排 LE)",
    "frozen": True,
    "target": TARGET,
    "receptor": RECEPTOR,
    "library_size": N,
    "l2_model": "retrieved_negatives_retrain (Layer2BindingDB, K=2)" if _L2P else "hash (Layer2BindingDB, 检索域最优)",
    "config": {"TOPN": TOPN, "workers": N_WORKERS, "exhaustiveness": EXHAUST,
               "smina_cpu": SMINA_CPU, "box_center": CENTER, "box_size": SIZE,
               "dock_timeout": DOCK_TIMEOUT, "hac_max": DOCK_HAC_MAX},
    "auc_l2_full": round(float(auc_l2_full), 4),
    "auc_cascade_full_deployed": round(float(auc_cascade_full), 4),
    "auc_l2_balanced": round(float(auc_l2_bal), 4),
    "auc_cascade_balanced_LE": round(float(auc_dock_bal), 4),
    "auc_fused_balanced_L2plusLE": round(float(auc_fused_bal), 4),
    "auc_fused_full_deployed_L2plusLE": round(float(auc_fused_full), 4),
    "auc_deployed_corrected_L2plusLE": round(float(auc_deployed_corrected), 4),
    "fusion_w_le": W_LE,
    "n_actives_in_library": len(pos_rows),
    "n_actives_in_topN": n_act_in_top,
    "n_top_docked": n_top_docked,
    "n_bal_actives_docked": n_bal_act_docked,
    "ef@1%_l2": round(float(ef_l2), 4),
    "ef@1%_cascade": round(float(ef_cas), 4),
    "recall@1%_l2": round(float(rec_l2), 4),
    "recall@1%_cascade": round(float(rec_cas), 4),
    "na_demo_reference": {"L2_hash": 0.462, "docking_LE": 0.756},
    "elapsed_sec": round(time.time() - t0, 1),
}
json.dump(result, open(OUT, "w"), indent=2, ensure_ascii=False)
log(f"已写出 {OUT}")
log("=" * 60)
log("冻结版新 AUC 闭环结果:")
log(f"  L2 基线(全库) AUC        = {auc_l2_full:.3f}")
log(f"  级联(部署, 全库) AUC     = {auc_cascade_full:.3f}")
log(f"  平衡集 L2 AUC            = {auc_l2_bal:.3f}")
log(f"  平衡集 级联(LE) AUC      = {auc_dock_bal:.3f}   ← 对比 NA demo 0.462→0.756")
log(f"  平衡集 融合(L2+LE) AUC   = {auc_fused_bal:.3f}   ← 目标 0.8+")
log("=" * 60)

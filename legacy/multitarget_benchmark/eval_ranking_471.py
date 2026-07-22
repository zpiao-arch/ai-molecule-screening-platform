# -*- coding: utf-8 -*-
"""
评测 BPR 重训 L2: 与 bigrun_score_hash.py 同口径 (label_map + 100k 库)。
因子化前向 (W1 拆 mol/tgt), 已验证与 sklearn MLP 一致。
产出: bigrun_results_ranking.json (per_target_metrics: auc, n_pos, ...)
"""
import os, sys, time, json
from pathlib import Path
import numpy as np

REPO = Path("<validated-workspace>")
BR = Path("<external-library-cache>")
MODEL = REPO / os.environ.get("RANK_MODEL",
        os.environ.get("MODEL",
        "评分_work_package/评分/models/bindingdb_l2_ranking/l2_model.pt"))
OUT = REPO / "scientific_validation/multitarget_benchmark/bigrun_results_ranking_v2.json"
MIN_POS = int(os.environ.get("MIN_POS", "3"))

sys.path.insert(0, str(REPO / "评分_work_package/评分"))
import torch
from l2_bindingdb import BindingDBFeature

def sigmoid(x): return 1.0 / (1.0 + np.exp(-x))

def main():
    t0 = time.time()
    sd_bundle = torch.load(MODEL, map_location="cpu", weights_only=False)
    sd = {k: np.asarray(v, dtype=np.float32) for k, v in sd_bundle["state_dict"].items()}
    # fc1: (H1, D_IN) -> W1m (H1, D_MOL), W1t (H1, D_TGT)
    W1 = sd["fc1.weight"]; b1 = sd["fc1.bias"]
    D_MOL = sd_bundle["arch"]["d_mol"]; D_TGT = sd_bundle["arch"]["d_tgt"]
    W1m = W1[:, :D_MOL].T          # (D_MOL, H1)
    W1t = W1[:, D_MOL:].T          # (D_TGT, H1)
    W2 = sd["fc2.weight"].T        # (H1, H2)
    b2 = sd["fc2.bias"]
    W3 = sd["fc3.weight"].T        # (H2, 1)
    b3 = sd["fc3.bias"]
    H1 = W1m.shape[1]

    feats = np.load(BR / "lib_feats.npy").astype(np.float32)     # (100000, 520)
    smiles = [l.strip() for l in open(BR / "lib_smiles.txt")]
    smi2row = {s: i for i, s in enumerate(smiles)}
    N_MOL = len(smiles)
    label_map = json.load(open(BR / "label_map.json"))
    feat_obj = BindingDBFeature()

    def target_hash(t):
        # 与训练同口径: BindingDBFeature.target_features
        return feat_obj.target_features(t).astype(np.float32)

    # 因子化: molpart 一次算完全库
    molpart = feats @ W1m + b1      # (N_MOL, H1)
    metrics = {}
    aucs = []
    for tid, d in label_map.items():
        if "pos" not in d: continue
        pos_rows = [smi2row[s] for s in d["pos"] if s in smi2row]
        npos = len(pos_rows)
        if npos < MIN_POS: continue
        tt = target_hash(tid)
        tgtpart = tt @ W1t + b1      # (H1,)
        h1 = np.maximum(molpart + tgtpart[None, :], 0.0)
        h2 = np.maximum(h1 @ W2 + b2, 0.0)
        logit = h2 @ W3 + b3
        p = sigmoid(logit.ravel())   # (N_MOL,)
        y = np.zeros(N_MOL, dtype=np.int8); y[pos_rows] = 1
        try:
            from sklearn.metrics import roc_auc_score
            auc = float(roc_auc_score(y, p))
        except Exception:
            auc = None
        if auc is not None:
            aucs.append(auc)
            metrics[tid] = {"auc": round(auc, 4), "n_pos": npos}
    print(f"[{time.time()-t0:.0f}s] 评测靶点 {len(aucs)}, 中位 AUC = {np.median(aucs):.4f}")
    print(f"  均值 {np.mean(aucs):.4f}  反向 {sum(a<0.5 for a in aucs)}  AUC>=0.7 {sum(a>=0.7 for a in aucs)}")
    json.dump({"summary": {"n_labeled_eval_targets": len(aucs),
                           "auc_median": round(float(np.median(aucs)),4),
                           "auc_mean": round(float(np.mean(aucs)),4),
                           "auc_reversed": sum(a<0.5 for a in aucs),
                           "auc_ge_0.7": sum(a>=0.7 for a in aucs),
                           "model": "BPR ranking retrain (torch MLP 776->256->128->1)"},
                "per_target_metrics": metrics}, open(OUT, "w"), indent=1)
    print(f"[{time.time()-t0:.0f}s] 已保存 -> {OUT}")

if __name__ == "__main__":
    main()

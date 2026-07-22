# -*- coding: utf-8 -*-
"""对比 旧哈希模型 vs 重训模型 在 100k 药库上的打分分布 + NA 靶 AUC, 诊断重训为何 AUC 下降。"""
import os, sys, json, time
import numpy as np
from pathlib import Path
REPO = Path("<validated-workspace>")
SC = REPO / "评分_work_package/评分"
sys.path.insert(0, str(SC))
from l2_bindingdb import Layer2BindingDB, BindingDBFeature
from target_resolver import TargetResolver
from sklearn.metrics import roc_auc_score

BR = Path("<external-library-cache>")
feats = np.load(BR / "lib_feats.npy")          # (100000, 520)
smiles = open(BR / "lib_smiles.txt").read().splitlines()
label_map = json.load(open(BR / "label_map.json"))
smi2row = {s: i for i, s in enumerate(smiles)}
N = feats.shape[0]

f = BindingDBFeature(); resolver = TargetResolver()
# 选 NA (CHEMBL2051) 做具体靶诊断
TID = "CHEMBL2051"
tt, _ = resolver.resolve(TID, chembl_id=TID)
th = f.target_features(tt)
X = np.concatenate([feats, np.tile(th, (N, 1))], axis=1).astype(np.float32)

def evaluate(model_path, tag):
    l2 = Layer2BindingDB(model_path=model_path, prefer="mlp")
    p = l2.predict_proba(X)
    # 分布
    frac_hi = float((p > 0.9).mean())
    frac_lo = float((p < 0.1).mean())
    med = float(np.median(p))
    # NA 标注 AUC
    pos_rows = [smi2row[s] for s in label_map[TID]["pos"] if s in smi2row]
    y = np.zeros(N, dtype=int); y[pos_rows] = 1
    auc = float(roc_auc_score(y, p)) if len(pos_rows) >= 1 else float("nan")
    print(f"[{tag}] NA AUC={auc:.4f}  frac>0.9={frac_hi:.3f}  frac<0.1={frac_lo:.3f}  median={med:.3f}  "
          f"p-range=[{p.min():.3f},{p.max():.3f}]")
    return p

OLD = str(SC / "models/bindingdb_l2/l2_model.joblib")
RET = "<validated-workspace>/评分_work_package/评分/models/bindingdb_l2_retrieval/l2_model.joblib"
p_old = evaluate(OLD, "OLD")
p_ret = evaluate(RET, "RET")
# 排序一致性(两模型对同一库排序的 Spearman-ish): 用 AUC of one vs other 的排序
from scipy.stats import spearmanr
rho, _ = spearmanr(p_old, p_ret)
print(f"[compare] NA 上 两模型打分 Spearman rho = {rho:.4f}")

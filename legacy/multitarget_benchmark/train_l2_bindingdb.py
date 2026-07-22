# -*- coding: utf-8 -*-
"""
训练 BindingDB L2 靶点感知模型, 保存权重供 MoleculeScorer 调用。

- 特征: 自洽 BindingDBFeature (Morgan512 + 8描述符 + 256靶点哈希 = 776)
- 数据: data_lake/bindingdb/aligned_model_input/bindingdb_202606_target_match_examples.csv
- 产出: 评分_work_package/评分/models/bindingdb_l2/l2_model.joblib
        (bundle: {"mlp":..., "logreg":..., "meta":...})
"""
import os, sys, json, time
from pathlib import Path

import numpy as np

REPO = Path("<validated-workspace>")
EXAMPLES = REPO / "data_lake/bindingdb/aligned_model_input/bindingdb_202606_target_match_examples.csv"
sys.path.insert(0, str(REPO / "评分_work_package/评分"))
from l2_bindingdb import BindingDBFeature, MODELS_DIR

def main():
    t0 = time.time()
    feat = BindingDBFeature()
    MODELS_DIR.mkdir(parents=True, exist_ok=True)

    print(f"[{time.time()-t0:.0f}s] 读取 examples.csv ...")
    import csv
    rows = []
    with open(EXAMPLES) as f:
        r = csv.DictReader(f)
        for row in r:
            rows.append((row["canonical_smiles"], row["target_text"], int(float(row["label"]))))

    # 去重: 分子特征按 smiles, 靶点特征按 target_text
    smiles_set, target_set = set(), set()
    for smi, txt, _ in rows:
        smiles_set.add(smi); target_set.add(txt)
    print(f"[{time.time()-t0:.0f}s] 唯一分子={len(smiles_set)}  唯一靶点文本={len(target_set)}")

    print(f"[{time.time()-t0:.0f}s] 计算分子特征 (512+8) ...")
    mol_feat = {}
    fail = 0
    for i, smi in enumerate(smiles_set):
        mf = feat.mol_features(smi)
        if mf is None:
            fail += 1
        else:
            mol_feat[smi] = mf
    print(f"[{time.time()-t0:.0f}s] 分子特征完成, 失败={fail}")

    print(f"[{time.time()-t0:.0f}s] 计算靶点哈希 (256) ...")
    tgt_feat = {txt: feat.target_features(txt) for txt in target_set}
    print(f"[{time.time()-t0:.0f}s] 靶点特征完成")

    # 组装 X, y (跳过分子解析失败的样本)
    Xlist, ylist = [], []
    dropped = 0
    for smi, txt, lab in rows:
        mf = mol_feat.get(smi)
        if mf is None:
            dropped += 1; continue
        tf = tgt_feat[txt]
        Xlist.append(np.concatenate([mf, tf]))
        ylist.append(lab)
    X = np.array(Xlist, dtype=np.float32)
    y = np.array(ylist, dtype=np.int64)
    print(f"[{time.time()-t0:.0f}s] 组装 X={X.shape}, y pos={int(y.sum())} neg={int((1-y).sum())}, dropped={dropped}")

    # 训练 / 评估
    from sklearn.model_selection import train_test_split
    from sklearn.linear_model import LogisticRegression
    from sklearn.neural_network import MLPClassifier
    from sklearn.metrics import roc_auc_score

    Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=0.2, random_state=0, stratify=y)

    print(f"[{time.time()-t0:.0f}s] 训练 LogisticRegression ...")
    lr = LogisticRegression(max_iter=400, class_weight="balanced", C=1.0)
    lr.fit(Xtr, ytr)
    auc_lr = roc_auc_score(yte, lr.predict_proba(Xte)[:, 1])

    print(f"[{time.time()-t0:.0f}s] 训练 MLP(256,128) ...")
    mlp = MLPClassifier(hidden_layer_sizes=(256, 128), max_iter=40, alpha=1e-4,
                        learning_rate_init=1e-3, batch_size=1024, random_state=0,
                        early_stopping=True, validation_fraction=0.1)
    mlp.fit(Xtr, ytr)
    auc_mlp = roc_auc_score(yte, mlp.predict_proba(Xte)[:, 1])

    print(f"[{time.time()-t0:.0f}s] 留出 AUC: LR={auc_lr:.4f}  MLP={auc_mlp:.4f}")

    # 保存
    try:
        import joblib
        bundle = {"mlp": mlp, "logreg": lr,
                  "meta": {"input_dim": int(X.shape[1]), "auc_logreg": round(float(auc_lr),4),
                           "auc_mlp": round(float(auc_mlp),4),
                           "n_train": int(Xtr.shape[0]), "n_test": int(Xte.shape[0]),
                           "feature": "Morgan r2 512 + 8 RDKit desc + 256 target FeatureHash",
                           "trained_at": time.strftime("%Y-%m-%dT%H:%M:%S")}}
        joblib.dump(bundle, MODELS_DIR / "l2_model.joblib")
        print(f"[{time.time()-t0:.0f}s] 已保存 -> {MODELS_DIR/'l2_model.joblib'}")
    except Exception as e:
        print("保存失败:", e); raise

    json.dump(bundle["meta"], open(MODELS_DIR / "l2_params.json", "w"), indent=1)
    print("done.")

if __name__ == "__main__":
    main()

# -*- coding: utf-8 -*-
"""
A. 检索式负采样重训 L2 (治本: 对齐训练/部署分布)

问题根因 (见 diagnose_na_recall_report.md / bigrun_hash_report.md):
  原模型训练负例来自 examples.csv 的 label=0 (测定非结合物),
  但部署是在 100k 真实药库上做检索排序。训练/部署分布错配 ->
  测得域留出 AUC 0.92, 检索域大通量中位仅 0.689 (NA 甚至 0.485)。

修复:
  正例  = examples.csv 测定活性 (label=1, 34090)
  负例  = 从 100k 真实药库 (lib_feats) 按靶点采样, 排除该靶自身活性
          -> 负例 = 部署分布下的"药物样但不结合该靶"难负例
  这样训练负例 == 部署背景, 让模型直接优化检索排序任务。

产出: 评分_work_package/评分/models/bindingdb_l2_retrieval/l2_model.joblib
      (与 Layer2BindingDB 兼容的 bundle: {mlp, logreg, meta})
"""
import os, sys, json, time, csv
from pathlib import Path
import numpy as np

REPO = Path("<validated-workspace>")
EXAMPLES = REPO / "data_lake/bindingdb/aligned_model_input/bindingdb_202606_target_match_examples.csv"
LIB_FEATS = Path("<external-library-cache>/lib_feats.npy")
LIB_SMILES = Path("<external-library-cache>/lib_smiles.txt")
OUT_DIR = REPO / "评分_work_package/评分/models/bindingdb_l2_retrieval"
K = int(os.environ.get("RETRIEVAL_K", "2"))          # 每正例采样负例数 (1:2 近平衡)
SEED = int(os.environ.get("RETRIEVAL_SEED", "0"))

sys.path.insert(0, str(REPO / "评分_work_package/评分"))
from l2_bindingdb import BindingDBFeature

def main():
    t0 = time.time()
    feat = BindingDBFeature()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"[{time.time()-t0:.0f}s] 读库特征 ...")
    lib_feats = np.load(LIB_FEATS)                    # (100000, 520)
    lib_smiles = [l.strip() for l in open(LIB_SMILES)]
    smi2idx = {s: i for i, s in enumerate(lib_smiles)}
    Nlib = len(lib_smiles)
    print(f"[{time.time()-t0:.0f}s] 库 {Nlib} 分子, feat {lib_feats.shape}")

    print(f"[{time.time()-t0:.0f}s] 读 examples.csv 正例 ...")
    pos = {}                                          # target_text -> [smiles]
    null_cnt = 0
    with open(EXAMPLES) as f:
        for r in csv.DictReader(f):
            smi = (r.get("canonical_smiles") or "").strip()
            if not smi or smi.lower() == "null":
                null_cnt += 1
                continue
            if r["label"] == "1.0":
                pos.setdefault(r["target_text"], []).append(smi)
    print(f"[{time.time()-t0:.0f}s] 跳过 null/空 SMILES {null_cnt} 行")
    print(f"[{time.time()-t0:.0f}s] 正例靶点 {len(pos)}, 正例总数 {sum(len(v) for v in pos.values())}")

    tgt_hash = {}
    def gethash(t):
        if t not in tgt_hash:
            tgt_hash[t] = feat.target_features(t)
        return tgt_hash[t]

    # 分子特征: 库内优先用 lib_feats (与部署一致), 库外现算
    def getmol(smi):
        idx = smi2idx.get(smi)
        return lib_feats[idx] if idx is not None else feat.mol_features(smi)

    # 每靶自身活性索引集合 (负例需排除)
    act_idx = {}
    for t, smis in pos.items():
        act_idx[t] = set(smi2idx[s] for s in smis if s in smi2idx)

    rng = np.random.default_rng(SEED)
    Xlist, ylist = [], []
    n_act_out_lib = 0
    for ti, (t, smis) in enumerate(pos.items()):
        th = gethash(t)
        excl = act_idx.get(t, set())
        # 负例候选: 全库排除该靶活性
        mask = np.ones(Nlib, dtype=bool)
        if excl:
            mask[np.array(list(excl), dtype=int)] = False
        cand = np.where(mask)[0]
        n_neg_want = min(K * len(smis), len(cand))
        neg_idx = rng.choice(cand, size=n_neg_want, replace=False)
        # 正例
        for s in smis:
            mf = getmol(s)
            if mf is None:
                n_act_out_lib += 1
                continue
            Xlist.append(np.concatenate([mf, th])); ylist.append(1)
        # 负例 (药库分子)
        for i in neg_idx:
            Xlist.append(np.concatenate([lib_feats[i], th])); ylist.append(0)
        if (ti + 1) % 200 == 0:
            print(f"[{time.time()-t0:.0f}s] 处理靶点 {ti+1}/{len(pos)}  累计样本 {len(Xlist)}")

    X = np.array(Xlist, dtype=np.float32)
    y = np.array(ylist, dtype=np.int64)
    print(f"[{time.time()-t0:.0f}s] 组装完成 X={X.shape}  pos={int(y.sum())} neg={int((1-y).sum())} 库外活性(跳过)={n_act_out_lib}")

    from sklearn.model_selection import train_test_split
    from sklearn.linear_model import LogisticRegression
    from sklearn.neural_network import MLPClassifier
    from sklearn.metrics import roc_auc_score

    Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=0.15, random_state=SEED, stratify=y)
    print(f"[{time.time()-t0:.0f}s] 训练集 {Xtr.shape} / 测试集 {Xte.shape}")

    print(f"[{time.time()-t0:.0f}s] 训练 LogisticRegression (balanced) ...")
    lr = LogisticRegression(max_iter=600, class_weight="balanced", C=1.0)
    lr.fit(Xtr, ytr)
    auc_lr = roc_auc_score(yte, lr.predict_proba(Xte)[:, 1])

    print(f"[{time.time()-t0:.0f}s] 训练 MLP(256,128) ...")
    mlp = MLPClassifier(hidden_layer_sizes=(256, 128), max_iter=40, alpha=1e-4,
                        learning_rate_init=1e-3, batch_size=1024, random_state=SEED,
                        early_stopping=True, validation_fraction=0.1)
    mlp.fit(Xtr, ytr)
    auc_mlp = roc_auc_score(yte, mlp.predict_proba(Xte)[:, 1])
    print(f"[{time.time()-t0:.0f}s] 留出 AUC: LR={auc_lr:.4f}  MLP={auc_mlp:.4f}")

    import joblib
    bundle = {"mlp": mlp, "logreg": lr,
              "meta": {"input_dim": int(X.shape[1]),
                       "auc_logreg": round(float(auc_lr), 4),
                       "auc_mlp": round(float(auc_mlp), 4),
                       "n_train": int(Xtr.shape[0]), "n_test": int(Xte.shape[0]),
                       "neg_scheme": f"retrieval_from_100k_library K={K}",
                       "feature": "Morgan r2 512 + 8 RDKit desc + 256 target FeatureHash",
                       "trained_at": time.strftime("%Y-%m-%dT%H:%M:%S")}}
    joblib.dump(bundle, OUT_DIR / "l2_model.joblib")
    json.dump(bundle["meta"], open(OUT_DIR / "l2_params.json", "w"), indent=1)
    print(f"[{time.time()-t0:.0f}s] 已保存 -> {OUT_DIR/'l2_model.joblib'}  done.")

if __name__ == "__main__":
    main()

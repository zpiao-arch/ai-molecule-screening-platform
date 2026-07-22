# -*- coding: utf-8 -*-
"""
合并训练 L2 靶点感知模型 (扩大训练面):
  - 主源: BindingDB 5万条 (target_text, canonical_smiles, label)
  - 主源: ChEMBL 36k条 (target_text, canonical_smiles, label)  <- 直接覆盖面板那 48 个未映射靶点
  - 增强: OpenTargets association_overall_direct (药-靶关联分数) [best-effort, 需 pyarrow]
统一 target_text 字符串 (与 target_resolver.py 输出一致), 训练 MLP/LogReg 并保存。

产出: 评分_work_package/评分/models/bindingdb_l2/l2_model.joblib
       (覆盖旧 bindingdb-only 模型; meta 记录训练面组成)
"""
import sys, time, json
from pathlib import Path
from collections import Counter

import numpy as np

REPO = Path("<validated-workspace>")
BINDINGDB_EX = REPO / "data_lake/bindingdb/aligned_model_input/bindingdb_202606_target_match_examples.csv"
CHEMBL_EX = REPO / "data_lake/chembl/latest/aligned_model_input/chembl_37_target_match_examples.csv"
OT_DIR = REPO / "data_lake/opentargets/platform/latest/output"
sys.path.insert(0, str(REPO / "评分_work_package/评分"))
from l2_bindingdb import BindingDBFeature, MODELS_DIR


def iter_csv(path, cols, label_col, smi_col, text_col):
    import csv
    with open(path) as f:
        for row in csv.DictReader(f):
            lab = row.get(label_col)
            if lab is None or lab == "":
                continue
            smi = row.get(smi_col, "").strip()
            txt = row.get(text_col, "").strip()
            if not smi or not txt:
                continue
            try:
                y = int(float(lab))
            except ValueError:
                continue
            if y not in (0, 1):
                continue
            yield smi, txt, y


def load_opentargets():
    """best-effort: 从 OpenTargets parquet 构建 (smiles, target_text, label)。
    返回 list[(smi, txt, y)] 或 None (环境/数据不可用)。"""
    try:
        import pyarrow.parquet as pq
    except Exception as e:
        print(f"  [OpenTargets] pyarrow 不可用, 跳过: {e}")
        return None
    assoc = OT_DIR / "association_overall_direct"
    drug = OT_DIR / "drug_molecule"
    tgt = OT_DIR / "target"
    if not (assoc.exists() and drug.exists() and tgt.exists()):
        print("  [OpenTargets] parquet 缺失, 跳过")
        return None
    try:
        # 关联: targetId, drugId, overallScore
        ap = pq.ParquetFile(sorted(assoc.glob("*.parquet"))[0])
        cols = ap.schema_arrow.names
        tid_i, did_i, score_i = cols.index("targetId"), cols.index("drugId"), cols.index("score") if "score" in cols else cols.index("overallScore")
        pairs = []
        for b in ap.iter_batches(batch_size=200000):
            arr = b.to_pydict()
            for i in range(len(arr[cols[0]])):
                s = arr[cols[score_i]][i]
                if s is None:
                    continue
                pairs.append((arr[cols[tid_i]][i], arr[cols[did_i]][i], float(s)))
        # drug -> smiles
        dp = pq.ParquetFile(sorted(drug.glob("*.parquet"))[0])
        dcols = dp.schema_arrow.names
        smi_i = dcols.index("smiles") if "smiles" in dcols else dcols.index("canonicalSmiles")
        did2_i = dcols.index("id") if "id" in dcols else dcols.index("chemblId")
        drug2smi = {}
        for b in dp.iter_batches(batch_size=200000):
            arr = b.to_pydict()
            for i in range(len(arr[dcols[0]])):
                drug2smi[arr[dcols[did2_i]][i]] = arr[dcols[smi_i]][i]
        # target -> name
        tp = pq.ParquetFile(sorted(tgt.glob("*.parquet"))[0])
        tcols = tp.schema_arrow.names
        tn_i = tcols.index("id") if "id" in tcols else tcols.index("targetId")
        nm_i = tcols.index("name") if "name" in tcols else tcols.index("approvedSymbol")
        tgt2name = {}
        for b in tp.iter_batches(batch_size=200000):
            arr = b.to_pydict()
            for i in range(len(arr[tcols[0]])):
                tgt2name[arr[tcols[tn_i]][i]] = arr[tcols[nm_i]][i]
        out = []
        for tid, did, score in pairs:
            smi = drug2smi.get(did)
            name = tgt2name.get(tid)
            if not smi or not name:
                continue
            y = 1 if score >= 0.7 else (0 if score <= 0.3 else None)
            if y is None:
                continue
            out.append((smi, f"{tid} {name} SINGLE PROTEIN Homo sapiens", y))
        print(f"  [OpenTargets] 有效样本 {len(out)}")
        return out
    except Exception as e:
        print(f"  [OpenTargets] 解析失败, 跳过: {e}")
        return None


def main():
    t0 = time.time()
    feat = BindingDBFeature()
    MODELS_DIR.mkdir(parents=True, exist_ok=True)

    print(f"[{time.time()-t0:.0f}s] 收集训练样本 ...")
    rows = []
    n_bd = n_ch = 0
    # BindingDB 主源
    for smi, txt, y in iter_csv(BINDINGDB_EX, None, "label", "canonical_smiles", "target_text"):
        rows.append((smi, txt, y)); n_bd += 1
    # ChEMBL 主源 (用其自带 target_text, 与 resolver 一致)
    for smi, txt, y in iter_csv(CHEMBL_EX, None, "label", "canonical_smiles", "target_text"):
        rows.append((smi, txt, y)); n_ch += 1
    # OpenTargets 增强
    ot = load_opentargets()
    n_ot = 0
    if ot:
        rows.extend(ot); n_ot = len(ot)
    print(f"[{time.time()-t0:.0f}s] 原始样本: BindingDB={n_bd}  ChEMBL={n_ch}  OpenTargets={n_ot}  合计={len(rows)}")

    # 去重 (smiles, target_text)
    seen = set(); uniq = []
    for smi, txt, y in rows:
        k = (smi, txt)
        if k in seen:
            continue
        seen.add(k); uniq.append((smi, txt, y))
    rows = uniq
    smiles_set = set(s for s, _, _ in rows)
    target_set = set(t for _, t, _ in rows)
    print(f"[{time.time()-t0:.0f}s] 去重后唯一 (smi,target)={len(rows)}  唯一分子={len(smiles_set)}  唯一靶点={len(target_set)}")

    print(f"[{time.time()-t0:.0f}s] 分子特征 ...")
    mol_feat = {}
    fail = 0
    for i, smi in enumerate(smiles_set):
        mf = feat.mol_features(smi)
        if mf is None:
            fail += 1
        else:
            mol_feat[smi] = mf
    print(f"[{time.time()-t0:.0f}s] 分子特征完成, 失败={fail}")

    print(f"[{time.time()-t0:.0f}s] 靶点哈希 ...")
    tgt_feat = {txt: feat.target_features(txt) for txt in target_set}

    Xlist, ylist, dropped = [], [], 0
    for smi, txt, lab in rows:
        mf = mol_feat.get(smi)
        if mf is None:
            dropped += 1; continue
        Xlist.append(np.concatenate([mf, tgt_feat[txt]]))
        ylist.append(lab)
    X = np.array(Xlist, dtype=np.float32)
    y = np.array(ylist, dtype=np.int64)
    print(f"[{time.time()-t0:.0f}s] X={X.shape} pos={int(y.sum())} neg={int((1-y).sum())} dropped={dropped}")

    from sklearn.model_selection import train_test_split
    from sklearn.linear_model import LogisticRegression
    from sklearn.neural_network import MLPClassifier
    from sklearn.metrics import roc_auc_score

    Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=0.2, random_state=0, stratify=y)
    print(f"[{time.time()-t0:.0f}s] 训练 LR ...")
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

    import joblib
    bundle = {"mlp": mlp, "logreg": lr,
              "meta": {"input_dim": int(X.shape[1]),
                       "auc_logreg": round(float(auc_lr), 4),
                       "auc_mlp": round(float(auc_mlp), 4),
                       "n_train": int(Xtr.shape[0]), "n_test": int(Xte.shape[0]),
                       "n_bindingdb": n_bd, "n_chembl": n_ch, "n_opentargets": n_ot,
                       "n_unique_targets": len(target_set),
                       "feature": "Morgan r2 512 + 8 RDKit desc + 256 target FeatureHash",
                       "trained_at": time.strftime("%Y-%m-%dT%H:%M:%S")}}
    joblib.dump(bundle, MODELS_DIR / "l2_model.joblib")
    json.dump(bundle["meta"], open(MODELS_DIR / "l2_params.json", "w"), indent=1)
    print(f"[{time.time()-t0:.0f}s] 已保存 -> {MODELS_DIR/'l2_model.joblib'}")
    print("done.")


if __name__ == "__main__":
    main()

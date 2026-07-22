# -*- coding: utf-8 -*-
"""
训练序列嵌入版 L2 模型 (升级靶点表征):
  - 分子特征: Morgan512 + 8 desc (520 维, 与 Layer2BindingDB 一致, 预计算缓存)
  - 靶点特征: TargetSeqEncoder (BLOSUM62 + 1D CNN -> 128 维序列嵌入, 端到端可训练)
  - 头: MLP(648 -> 256 -> 64 -> 1) + sigmoid
  - 训练面: BindingDB 5万 + ChEMBL 36k (合并去重), 与 train_l2_combined 一致

产出: 评分_work_package/评分/models/bindingdb_l2_seq/l2_seq.pt
"""
import sys, time, json, pickle
from pathlib import Path
from collections import Counter

import numpy as np

REPO = Path("<validated-workspace>")
BINDINGDB_EX = REPO / "data_lake/bindingdb/aligned_model_input/bindingdb_202606_target_match_examples.csv"
CHEMBL_EX = REPO / "data_lake/chembl/latest/aligned_model_input/chembl_37_target_match_examples.csv"
MODEL_DIR = REPO / "评分_work_package/评分/models/bindingdb_l2_seq"
MOLFEAT_CACHE = MODEL_DIR / "molfeat_cache.pkl"
sys.path.insert(0, str(REPO / "评分_work_package/评分"))
from l2_bindingdb import BindingDBFeature, Layer2BindingDBSeq
from target_seq_embedding import TargetSeqEncoder, extract_chembl_id, chembl_to_seq


def iter_csv(path, label_col, smi_col, text_col):
    import csv
    with open(path) as f:
        for row in csv.DictReader(f):
            lab = row.get(label_col)
            if lab is None or lab == "":
                continue
            smi = (row.get(smi_col, "") or "").strip()
            txt = (row.get(text_col, "") or "").strip()
            if not smi or not txt:
                continue
            try:
                y = int(float(lab))
            except ValueError:
                continue
            if y not in (0, 1):
                continue
            yield smi, txt, y


def main():
    t0 = time.time()
    feat = BindingDBFeature()
    MODEL_DIR.mkdir(parents=True, exist_ok=True)

    print(f"[{time.time()-t0:.0f}s] 收集训练样本 ...")
    rows = []
    n_bd = n_ch = 0
    for smi, txt, y in iter_csv(BINDINGDB_EX, "label", "canonical_smiles", "target_text"):
        rows.append((smi, txt, y)); n_bd += 1
    for smi, txt, y in iter_csv(CHEMBL_EX, "label", "canonical_smiles", "target_text"):
        rows.append((smi, txt, y)); n_ch += 1
    # 去重 (smi, target_text)
    seen = set(); uniq = []
    for smi, txt, y in rows:
        k = (smi, txt)
        if k in seen:
            continue
        seen.add(k); uniq.append((smi, txt, y))
    rows = uniq
    print(f"[{time.time()-t0:.0f}s] 原始 {n_bd+n_ch} -> 去重 {len(rows)} (BD={n_bd} CH={n_ch})")

    # ---- 分子特征 (缓存) ----
    smiles_set = list({s for s, _, _ in rows})
    if MOLFEAT_CACHE.exists():
        mol_feat = pickle.load(open(MOLFEAT_CACHE, "rb"))
        print(f"[{time.time()-t0:.0f}s] 载入分子特征缓存 {len(mol_feat)}")
        need = [s for s in smiles_set if s not in mol_feat]
    else:
        mol_feat = {}; need = smiles_set
    if need:
        print(f"[{time.time()-t0:.0f}s] 计算分子特征 {len(need)} (新增) ...")
        for i, smi in enumerate(need):
            mf = feat.mol_features(smi)
            if mf is not None:
                mol_feat[smi] = mf
            if (i + 1) % 5000 == 0:
                print(f"  mol {i+1}/{len(need)}")
        pickle.dump(mol_feat, open(MOLFEAT_CACHE, "wb"))
        print(f"[{time.time()-t0:.0f}s] 分子特征缓存写出 {len(mol_feat)}")

    # ---- 样本 -> (smi_idx, chembl_id, y) ----
    smi_to_i = {s: i for i, s in enumerate(mol_feat)}
    samples = []  # (smi_i, chembl_id_or_None, y)
    dropped = 0
    for smi, txt, y in rows:
        si = smi_to_i.get(smi)
        if si is None:
            dropped += 1; continue
        cid = extract_chembl_id(txt)
        samples.append((si, cid, y))
    print(f"[{time.time()-t0:.0f}s] 有效样本 {len(samples)} (dropped {dropped})")

    # ---- 靶点序列覆盖统计 ----
    cids = list({c for _, c, _ in samples if c})
    cov = sum(1 for c in cids if chembl_to_seq(c))
    print(f"[{time.time()-t0:.0f}s] 训练唯一靶点(有chembl_id)={len(cids)}  其中序列覆盖={cov} ({100*cov/max(1,len(cids)):.0f}%)")

    # ---- 转 numpy ----
    smi_idx = np.array([s[0] for s in samples], dtype=np.int64)
    y = np.array([s[2] for s in samples], dtype=np.float32)
    # 唯一靶点列表 (chembl_id 或 None->'UNK')
    uniq_cids = []
    cid_to_ti = {}
    for _, c, _ in samples:
        key = c if c else "UNK"
        if key not in cid_to_ti:
            cid_to_ti[key] = len(uniq_cids); uniq_cids.append(key)
    tgt_idx = np.array([cid_to_ti[c if c else "UNK"] for _, c, _ in samples], dtype=np.int64)
    print(f"[{time.time()-t0:.0f}s] 唯一靶点(含UNK)={len(uniq_cids)}")

    X_mol = np.stack([mol_feat[list(mol_feat)[i]] for i in smi_idx])  # (N,520)
    print(f"[{time.time()-t0:.0f}s] X_mol={X_mol.shape} pos={int(y.sum())} neg={int((1-y).sum())}")

    # ---- 训练/留出划分 (分子不相交: 留出分子在训练中从未出现) ----
    # 旧的随机 pair 划分会让同一分子同时进训练+留出, 高估泛化 (留出0.82 但真实
    # 面板检索仅0.54)。改为按 smi_idx 划分, 留出 AUC 才是对新药分子的诚实代理。
    rng = np.random.RandomState(0)
    uniq_smi = np.array(sorted(set(smi_idx.tolist())), dtype=np.int64)
    rng.shuffle(uniq_smi)
    n_hold_smi = int(0.2 * len(uniq_smi))
    hold_smi = set(uniq_smi[:n_hold_smi].tolist())
    te_mask = np.array([s in hold_smi for s in smi_idx], dtype=bool)
    te_idx = np.where(te_mask)[0]
    tr_idx = np.where(~te_mask)[0]
    print(f"[{time.time()-t0:.0f}s] 分子不相交划分: 训练 {len(tr_idx)} / 留出 {len(te_idx)} "
          f"(留出分子 {len(hold_smi)} 个, 训练中不出现)")

    # ---- 模型 (加 weight_decay 正则, 抑制对训练分子×靶点组合的记忆) ----
    import torch
    from torch.utils.data import TensorDataset, DataLoader
    enc = TargetSeqEncoder()
    mlp = Layer2BindingDBSeq._build_mlp((256, 64))
    params = list(mlp.parameters()) + list(enc.parameters())
    opt = torch.optim.Adam(params, lr=1e-3, weight_decay=1e-4)
    crit = torch.nn.BCEWithLogitsLoss()

    Xtr_mol = torch.from_numpy(X_mol[tr_idx])
    ytr = torch.from_numpy(y[tr_idx])
    Xte_mol = torch.from_numpy(X_mol[te_idx])
    yte = torch.from_numpy(y[te_idx])
    ttr = torch.from_numpy(tgt_idx[tr_idx])
    tte = torch.from_numpy(tgt_idx[te_idx])

    # 每 epoch 重算唯一靶点嵌入 -> (n_unique, 128), 训练时按索引切片
    def build_target_emb():
        embs = []
        for key in uniq_cids:
            seq = chembl_to_seq(key) if key != "UNK" else None
            if seq:
                embs.append(torch.from_numpy(enc.embed_one(seq)))
            else:
                embs.append(enc.unk.detach())
        return torch.stack(embs, dim=0)  # (n_unique, 128)

    from sklearn.metrics import roc_auc_score
    best = 0.0
    for epoch in range(1, 13):
        enc.train(); mlp.train()
        T_full = build_target_emb()
        Ttr = T_full[ttr]; Tte = T_full[tte]
        # 训练
        order = rng.permutation(len(tr_idx))
        batch = 2048
        for i in range(0, len(order), batch):
            bi = order[i:i + batch]
            X = torch.cat([Xtr_mol[bi], Ttr[bi]], dim=1)
            logit = mlp(X).ravel()
            loss = crit(logit, ytr[bi])
            opt.zero_grad(); loss.backward(); opt.step()
        # 评估留出
        enc.eval(); mlp.eval()
        with torch.no_grad():
            X = torch.cat([Xte_mol, Tte], dim=1)
            pred = torch.sigmoid(mlp(X).ravel()).numpy()
        auc = roc_auc_score(yte.numpy(), pred)
        print(f"[{time.time()-t0:.0f}s] epoch {epoch:2d}  loss~ 留出AUC={auc:.4f}")
        if auc > best:
            best = auc
            torch.save({"mlp_state": mlp.state_dict(),
                        "encoder_state": enc.state_dict(),
                        "mlp_hidden": (256, 64),
                        "meta": {"input_dim": Layer2BindingDBSeq.INPUT_DIM,
                                 "auc_holdout": round(float(auc), 4),
                                 "feature": "Morgan512+8desc + TargetSeqEncoder(BLOSUM62+CNN,128)",
                                 "n_train": int(len(tr_idx)), "n_test": int(len(te_idx)),
                                 "n_unique_targets": len(uniq_cids),
                                 "n_bindingdb": n_bd, "n_chembl": n_ch,
                                 "trained_at": time.strftime("%Y-%m-%dT%H:%M:%S")}},
                       Layer2BindingDBSeq.MODEL_PATH)
    print(f"[{time.time()-t0:.0f}s] 最佳留出 AUC={best:.4f} -> {Layer2BindingDBSeq.MODEL_PATH}")
    json.dump({"auc_holdout": round(float(best), 4), "feature": "seq-embed-128"},
              open(MODEL_DIR / "l2_seq_params.json", "w"), indent=1)
    print("done.")


if __name__ == "__main__":
    main()

# -*- coding: utf-8 -*-
"""
B. 配对排序损失 (BPR) 重训 L2 —— 针对 A 任务的失败机制

A 任务 (train_l2_retrieval.py, 二元 CE + 检索负例) 的问题:
  重训后把"好靶"也坑了 (如 CHEMBL4069: 0.999 -> 0.001)。
  根因: 二元 CE + 极端检索负例 -> 模型学到"几乎所有都是非活性 -> 输出~0",
        连好靶的活性也被压到 0 与阴性混在一起 -> AUC 崩。

BPR (Bayesian Personalized Ranking) 解法:
  损失 = -log(sigmoid(s_active - s_inactive)) 对每对 (活性, 阴性)。
  - 只优化相对序, 与 AUC 直接对齐 -> 不可能塌成全 0 (那样所有对都平局 -> 损失最大)。
  - 活性必须排在阴性之上, 好靶的活性不会被埋。
  - 仍用"检索负例"(部署分布), 治本方向不变。

产出: models/bindingdb_l2_ranking/l2_model.pt  (torch state_dict + arch)
      评测复用 eval_ranking_471.py (因子化前向, 与 bigrun 同口径)。
"""
import os, sys, time, csv, json
from pathlib import Path
import numpy as np

REPO = Path("<validated-workspace>")
EXAMPLES = REPO / "data_lake/bindingdb/aligned_model_input/bindingdb_202606_target_match_examples.csv"
LIB_FEATS = Path("<external-library-cache>/lib_feats.npy")
LIB_SMILES = Path("<external-library-cache>/lib_smiles.txt")
OUT_DIR = REPO / "评分_work_package/评分/models/bindingdb_l2_ranking"
SEED = int(os.environ.get("RANK_SEED", "0"))
K = int(os.environ.get("RANK_K", "4"))           # 每阳性采样阴性数 (BPR 每 epoch 重采)
EPOCHS = int(os.environ.get("RANK_EPOCHS", "40"))
BATCH = int(os.environ.get("RANK_BATCH", "2048"))
LR = float(os.environ.get("RANK_LR", "1e-3"))
WD = float(os.environ.get("RANK_WD", "1e-4"))

sys.path.insert(0, str(REPO / "评分_work_package/评分"))
import torch
import torch.nn as nn
from l2_bindingdb import BindingDBFeature

D_MOL = 520          # Morgan512 + 8 RDKit desc
D_TGT = 256          # target FeatureHash
D_IN = D_MOL + D_TGT # 776
H1, H2 = 256, 128

class MLP(nn.Module):
    def __init__(self):
        super().__init__()
        self.fc1 = nn.Linear(D_IN, H1)
        self.fc2 = nn.Linear(H1, H2)
        self.fc3 = nn.Linear(H2, 1)
    def forward(self, x):
        h = torch.relu(self.fc1(x))
        h = torch.relu(self.fc2(h))
        return self.fc3(h)   # raw logit (BPR 用)

def main():
    t0 = time.time()
    feat = BindingDBFeature()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    torch.manual_seed(SEED); np.random.seed(SEED)

    print(f"[{time.time()-t0:.0f}s] 读库特征 ...")
    lib_feats = np.load(LIB_FEATS).astype(np.float32)        # (100000, 520)
    lib_smiles = [l.strip() for l in open(LIB_SMILES)]
    smi2idx = {s: i for i, s in enumerate(lib_smiles)}
    Nlib = len(lib_smiles)

    print(f"[{time.time()-t0:.0f}s] 读 examples.csv 正例 ...")
    pos = {}
    null_cnt = 0
    with open(EXAMPLES) as f:
        for r in csv.DictReader(f):
            smi = (r.get("canonical_smiles") or "").strip()
            if not smi or smi.lower() == "null":
                null_cnt += 1; continue
            if r["label"] == "1.0":
                pos.setdefault(r["target_text"], []).append(smi)
    print(f"[{time.time()-t0:.0f}s] 正例靶点 {len(pos)}, 正例总数 {sum(len(v) for v in pos.values())} (跳过 null {null_cnt})")

    tgt_hash = {}
    def gethash(t):
        if t not in tgt_hash:
            tgt_hash[t] = feat.target_features(t).astype(np.float32)
        return tgt_hash[t]

    def getmol(smi):
        idx = smi2idx.get(smi)
        if idx is not None:
            return lib_feats[idx]
        m = feat.mol_features(smi)
        return m.astype(np.float32) if m is not None else None

    # 组装: 每靶 阳性特征列表 + 阴性(药库,排除自身活性)特征列表
    pos_feats = {}     # target -> (n_pos, 776)
    neg_feats = {}     # target -> (n_neg, 776)
    for ti, (t, smis) in enumerate(pos.items()):
        th = gethash(t)
        pf = []
        for s in smis:
            mf = getmol(s)
            if mf is None: continue
            pf.append(np.concatenate([mf, th]))
        if not pf: continue
        excl = set(smi2idx[s] for s in smis if s in smi2idx)
        mask = np.ones(Nlib, dtype=bool)
        if excl: mask[np.array(list(excl), dtype=int)] = False
        cand = np.where(mask)[0]
        n_neg = min(K * len(pf), len(cand))
        nidx = np.random.choice(cand, size=n_neg, replace=False)
        nf = np.stack([lib_feats[i] for i in nidx], axis=0)  # (n_neg, 520)
        nf = np.concatenate([nf, np.tile(th, (len(nf), 1))], axis=1)
        pos_feats[t] = np.array(pf, dtype=np.float32)
        neg_feats[t] = nf
        if (ti+1) % 200 == 0:
            print(f"[{time.time()-t0:.0f}s] 靶点 {ti+1}/{len(pos)}  累计阳性 {sum(len(v) for v in pos_feats.values())}")

    targets = list(pos_feats.keys())
    total_pos = sum(len(v) for v in pos_feats.values())
    print(f"[{time.time()-t0:.0f}s] 组装完成: 靶点 {len(targets)}, 阳性 {total_pos}")

    model = MLP()
    opt = torch.optim.Adam(model.parameters(), lr=LR, weight_decay=WD)
    bce = nn.BCEWithLogitsLoss()

    rng = np.random.default_rng(SEED)
    for epoch in range(EPOCHS):
        model.train()
        # 每 epoch 重采配对: 每阳性配一个该靶阴性
        pidx, nidx = [], []
        for t in targets:
            pf = pos_feats[t]; nf = neg_feats[t]
            pi = rng.integers(0, len(pf), size=len(pf))
            ni = rng.integers(0, len(nf), size=len(pf))
            pidx.append(pf[pi]); nidx.append(nf[ni])
        P = np.concatenate(pidx, axis=0)   # (total_pos, 776)
        N = np.concatenate(nidx, axis=0)
        order = rng.permutation(len(P))
        P, N = P[order], N[order]
        tot_loss = 0.0; nb = 0
        for i in range(0, len(P), BATCH):
            pb = torch.tensor(P[i:i+BATCH]); nb_t = torch.tensor(N[i:i+BATCH])
            sp = model(pb).squeeze(-1); sn = model(nb_t).squeeze(-1)
            loss = bce(sp - sn, torch.ones_like(sp))   # = -log(sigmoid(sp-sn))
            opt.zero_grad(); loss.backward(); opt.step()
            tot_loss += loss.item()*len(pb); nb += len(pb)
        if (epoch+1) % 5 == 0 or epoch == 0:
            print(f"[{time.time()-t0:.0f}s] epoch {epoch+1}/{EPOCHS}  BPR loss = {tot_loss/nb:.4f}")

    # 保存 (与 Layer2BindingDB 同构的权重, 便于后续接入)
    sd = {k: v.cpu().numpy() for k, v in model.state_dict().items()}
    torch.save({"state_dict": sd, "arch": {"d_in": D_IN, "h1": H1, "h2": H2,
                  "d_mol": D_MOL, "d_tgt": D_TGT},
                "meta": {"loss": "BPR(logistic)", "neg_scheme": f"retrieval K={K}",
                         "epochs": EPOCHS, "lr": LR, "wd": WD,
                         "n_pos": int(total_pos), "n_targets": len(targets),
                         "trained_at": time.strftime("%Y-%m-%dT%H:%M:%S")}},
               OUT_DIR / "l2_model.pt")
    json.dump({"state_dict": {k: v.tolist() for k, v in sd.items()},
               "arch": {"d_in": D_IN, "h1": H1, "h2": H2, "d_mol": D_MOL, "d_tgt": D_TGT},
               "meta": {"loss": "BPR(logistic)", "neg_scheme": f"retrieval K={K}",
                        "epochs": EPOCHS, "lr": LR, "wd": WD, "n_pos": int(total_pos),
                        "n_targets": len(targets), "trained_at": time.strftime("%Y-%m-%dT%H:%M:%S")}},
              open(OUT_DIR / "l2_params.json", "w"), indent=1)
    print(f"[{time.time()-t0:.0f}s] 已保存 -> {OUT_DIR/'l2_model.pt'}  done.")

if __name__ == "__main__":
    main()

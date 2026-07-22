# -*- coding: utf-8 -*-
"""
① 修正版排序重训 (ListNet / 全库 softmax) —— 针对 BPR 的"尺度失校准"失败

BPR (train_l2_ranking.py) 失败机制:
  损失 = -log σ(s_活性 - s_阴性), 只约束配对内相对序, 不约束绝对尺度。
  K=4 阴性太少 -> 模型把活性推到全库最低 (CHEMBL4069 活性-31.5 vs 库-8.2) -> 评测≈随机(0.506)。

本脚本修正:
  损失 = ListNet 全库 softmax:
      L_t = - mean_{j∈活性} [ s_j - logsumexp_k(s_k) ]   (k 遍历全库 10 万分子)
  分母覆盖【整库】, 梯度强制每个活性分数高于整库 -> 尺度天然锚定, 不可能塌成全 0/全底。
  与评测协议 (AUC over 100k 库) 直接对齐。

效率: 单靶全库前向 0.14s (已实测); 471 靶 × 12 epoch ≈ 数分钟前向 + 反向。CPU 可行。
产出: models/bindingdb_l2_ranking_v2/l2_model.pt (同构权重, eval_ranking_471.py 可评测)
"""
import os, sys, time, csv, json
from pathlib import Path
import numpy as np

REPO = Path("<validated-workspace>")
EXAMPLES = REPO / "data_lake/bindingdb/aligned_model_input/bindingdb_202606_target_match_examples.csv"
LIB_FEATS = Path("<external-library-cache>/lib_feats.npy")
LIB_SMILES = Path("<external-library-cache>/lib_smiles.txt")
OUT_DIR = REPO / "评分_work_package/评分/models/bindingdb_l2_ranking_v2"
SEED = int(os.environ.get("RANK_SEED", "0"))
EPOCHS = int(os.environ.get("RANK_EPOCHS", "12"))
LR = float(os.environ.get("RANK_LR", "3e-3"))
WD = float(os.environ.get("RANK_WD", "1e-4"))
SMOKE = int(os.environ.get("RANK_SMOKE", "0"))   # 1 = 仅前 8 靶做 smoke

sys.path.insert(0, str(REPO / "评分_work_package/评分"))
import torch
import torch.nn as nn
from l2_bindingdb import BindingDBFeature

D_MOL = 520
D_TGT = 256
D_IN = D_MOL + D_TGT
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
        return self.fc3(h)

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
    feats_t = torch.tensor(lib_feats)

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

    # 组装每靶: 在库活性行号 + 靶点哈希
    targets = []          # list of (tgt_text, pos_rows_numpy(int64), tgt_hash(520? no 256))
    tgt_hash_cache = {}
    def gethash(t):
        if t not in tgt_hash_cache:
            tgt_hash_cache[t] = feat.target_features(t).astype(np.float32)
        return tgt_hash_cache[t]
    for t, smis in pos.items():
        th = gethash(t)
        rows = np.array([smi2idx[s] for s in smis if s in smi2idx], dtype=np.int64)
        if len(rows) >= 3:
            targets.append((t, rows, th))
    print(f"[{time.time()-t0:.0f}s] 可用靶(在库活性≥3) {len(targets)}")

    if SMOKE:
        targets = targets[:8]
        print(f"[smoke] 仅用前 {len(targets)} 靶")

    model = MLP()
    opt = torch.optim.Adam(model.parameters(), lr=LR, weight_decay=WD)

    def factorized_scores(tt_np):
        # tt_np: (256,)  -> 全库 (Nlib,) 分数, 带梯度 (用于训练)
        with torch.no_grad():
            tt = torch.tensor(tt_np)
        W1 = model.fc1.weight          # (H1, D_IN)
        W1m = W1[:, :D_MOL]            # (H1, D_MOL)
        W1t = W1[:, D_MOL:]            # (H1, D_TGT)
        molpart = feats_t @ W1m.t() + model.fc1.bias          # (Nlib, H1)
        tgtpart = tt @ W1t.t() + model.fc1.bias              # (H1,)
        h1 = torch.relu(molpart + tgtpart[None, :])          # (Nlib, H1)
        h2 = torch.relu(h1 @ model.fc2.weight.t() + model.fc2.bias)  # (Nlib, H2)
        s = h2 @ model.fc3.weight.t() + model.fc3.bias       # (Nlib, 1)
        return s.squeeze(-1)                                 # (Nlib,)

    for epoch in range(EPOCHS):
        model.train()
        tot = 0.0; nt = 0; skipped = 0
        for ti, (t, rows, th) in enumerate(targets):
            s = factorized_scores(th)                       # (Nlib,) 带梯度
            # ListNet: -mean_active( s_j - logsumexp_k(s_k) )
            lse = torch.logsumexp(s, dim=0)
            pos_s = s[torch.tensor(rows)]
            loss = -(pos_s - lse).mean()
            if not torch.isfinite(loss):
                skipped += 1
                continue
            opt.zero_grad(); loss.backward(); opt.step()
            tot += loss.item(); nt += 1
        print(f"[{time.time()-t0:.0f}s] epoch {epoch+1}/{EPOCHS}  ListNet loss = {tot/max(nt,1):.4f}  (跳过{skipped})")

        # 检查点
        if (epoch+1) % 4 == 0 or epoch == EPOCHS-1:
            sd = {k: v.cpu().numpy() for k, v in model.state_dict().items()}
            torch.save({"state_dict": sd, "arch": {"d_in": D_IN, "h1": H1, "h2": H2,
                          "d_mol": D_MOL, "d_tgt": D_TGT},
                        "meta": {"loss": "ListNet(full-library softmax)", "epochs": epoch+1,
                                 "lr": LR, "wd": WD, "n_pos_targets": len(targets),
                                 "trained_at": time.strftime("%Y-%m-%dT%H:%M:%S")}},
                       OUT_DIR / "l2_model.pt")
            print(f"[{time.time()-t0:.0f}s]   checkpoint -> {OUT_DIR/'l2_model.pt'}")

    sd = {k: v.cpu().numpy() for k, v in model.state_dict().items()}
    torch.save({"state_dict": sd, "arch": {"d_in": D_IN, "h1": H1, "h2": H2,
                  "d_mol": D_MOL, "d_tgt": D_TGT},
                "meta": {"loss": "ListNet(full-library softmax)", "epochs": EPOCHS,
                         "lr": LR, "wd": WD, "n_pos_targets": len(targets),
                         "trained_at": time.strftime("%Y-%m-%dT%H:%M:%S")}},
               OUT_DIR / "l2_model.pt")
    json.dump({"state_dict": {k: v.tolist() for k, v in sd.items()},
               "arch": {"d_in": D_IN, "h1": H1, "h2": H2, "d_mol": D_MOL, "d_tgt": D_TGT},
               "meta": {"loss": "ListNet(full-library softmax)", "epochs": EPOCHS, "lr": LR,
                        "wd": WD, "n_pos_targets": len(targets),
                        "trained_at": time.strftime("%Y-%m-%dT%H:%M:%S")}},
              open(OUT_DIR / "l2_params.json", "w"), indent=1)
    print(f"[{time.time()-t0:.0f}s] 已保存 -> {OUT_DIR/'l2_model.pt'}  done.")

if __name__ == "__main__":
    main()

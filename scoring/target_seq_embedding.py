# -*- coding: utf-8 -*-
"""
靶点蛋白序列嵌入 (结构性靶点表征, 替代原 256 维 target_text 哈希)。

为什么做这个升级:
  原 L2 用 target_text 的 256 维 FeatureHash 做靶点表征 —— 不同靶点的向量
  高度相似、且不含任何生物学结构信息, 是 AUC 偏低/反向的根因(见
  auc_root_cause_report.md)。

本模块提供离线可用的"结构性靶点嵌入":
  - 残基编码: BLOSUM62 替代矩阵(20 维/残基, 固定, 含进化替代信息)
  - 序列编码器: 1D CNN (BLOSUM序列 -> 局部基序 -> 128 维靶点嵌入), 可随 L2 端到端训练
  - 未见于训练/序列缺失的靶点 -> 可学习的 [UNK] 嵌入

ESM-2 权重本环境无法离线获取(缓存为空、无外网), 故用 BLOSUM62+CNN 作为
最佳离线替代; 该接口 pluggable, 未来若有 ESM-2/ProtBERT 权重可直接替换
`TargetSeqEncoder` 的残基编码层。

依赖: torch (本地可用 2.8.0)
"""
from __future__ import annotations
import re
import os
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn

try:
    from .asset_integrity import verify_asset
except ImportError:  # PYTHONPATH=scoring compatibility.
    from asset_integrity import verify_asset

REPO = Path(os.environ.get("FOUR_LEVEL_ROOT", Path(__file__).resolve().parents[1])).resolve()

# ---------------------------------------------------------------------------
# BLOSUM62 (标准 20x20, 顺序 ARNDCQEGHILKMFPSTWYV)
# ---------------------------------------------------------------------------
BLOSUM62_ORDER = "ARNDCQEGHILKMFPSTWYV"
_B62 = [
    [ 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0],
    [-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3],
    [-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3],
    [-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3],
    [ 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1],
    [-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2],
    [-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2],
    [ 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3],
    [-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3],
    [-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3],
    [-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1],
    [-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2],
    [-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1],
    [-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1],
    [-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2],
    [ 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2],
    [ 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0],
    [-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3],
    [-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1],
    [ 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4],
]
# 索引: 20 标准 AA + X(未知, 全0) = 21 类
_RES2IDX = {a: i for i, a in enumerate(BLOSUM62_ORDER)}
_RES2IDX['X'] = 20
_RES2IDX['B'] = _RES2IDX['N']
_RES2IDX['Z'] = _RES2IDX['Q']
_RES2IDX['U'] = 20
_RES2IDX['O'] = 20
_RES2IDX['*'] = 20
_RES2IDX['-'] = 20

# BLOSUM 嵌入矩阵 (21 x 20), 最后一行 X 全 0
_B62_MAT = np.array(_B62 + [[0]*20], dtype=np.float32)

MAX_LEN = 1500          # 覆盖面板最长序列(1504)
EMBED_DIM = 128         # 靶点嵌入维度
RES_DIM = 20            # 残基编码维度 (= BLOSUM62 行)


def aa_to_indices(seq: str) -> np.ndarray:
    """AA 序列 -> 整数索引数组 (截断/填充到 MAX_LEN)。"""
    s = re.sub(r"[^A-Za-z]", "", seq.upper())
    idx = np.array([_RES2IDX.get(c, 20) for c in s[:MAX_LEN]], dtype=np.int64)
    if len(idx) == 0:
        idx = np.array([20], dtype=np.int64)
    return idx


# ---------------------------------------------------------------------------
# 序列字典: chembl_id -> uniprot -> 序列
# ---------------------------------------------------------------------------
_CHEMBL_TO_UNIPROT = None
_UNIPROT_TO_SEQ = None


def build_seq_dicts(fasta_gz=None, mapping_txt=None):
    """从 ChEMBL FASTA + chembl->uniprot 映射构建序列字典。
    返回 (chembl_to_uniprot:dict, uniprot_to_seq:dict)。"""
    fasta_gz = fasta_gz or os.environ.get("FOUR_LEVEL_CHEMBL_FASTA", str(REPO / "data" / "chembl_37.fa.gz"))
    mapping_txt = mapping_txt or os.environ.get("FOUR_LEVEL_CHEMBL_UNIPROT_MAPPING", str(REPO / "data" / "chembl_uniprot_mapping.txt"))
    cmap, useq = {}, {}
    # uniprot -> seq
    import gzip
    cur = None
    with gzip.open(fasta_gz, "rt") as f:
        for line in f:
            if line.startswith(">"):
                m = re.search(r"\[([A-Z0-9]+)\]", line)
                cur = m.group(1) if m else None
                if cur not in useq:
                    useq[cur] = ""
            elif cur is not None:
                useq[cur] += line.strip()
    # chembl -> uniprot
    with open(mapping_txt) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                cmap[parts[1]] = parts[0]
    return cmap, useq


def get_seq_dicts():
    """懒加载并缓存序列字典。"""
    global _CHEMBL_TO_UNIPROT, _UNIPROT_TO_SEQ
    if _CHEMBL_TO_UNIPROT is None:
        _CHEMBL_TO_UNIPROT, _UNIPROT_TO_SEQ = build_seq_dicts()
    return _CHEMBL_TO_UNIPROT, _UNIPROT_TO_SEQ


def chembl_to_seq(chembl_id: str):
    """chembl_id -> 蛋白序列 (str) 或 None。"""
    cmap, useq = get_seq_dicts()
    u = cmap.get(chembl_id)
    if u and useq.get(u) and len(useq[u]) > 10:
        return useq[u]
    return None


def extract_chembl_id(target_text: str):
    """从 target_text 提取第一个 CHEMBLxxxx (训练行 target_text 含 chembl id)。"""
    if not target_text:
        return None
    m = re.search(r"CHEMBL\d+", target_text)
    return m.group(0) if m else None


# ---------------------------------------------------------------------------
# 靶点序列编码器 (PyTorch, 端到端可训练)
# ---------------------------------------------------------------------------
class TargetSeqEncoder(nn.Module):
    """BLOSUM62 残基编码 + 1D CNN -> 128 维靶点嵌入。

    - 残基编码层: 固定 BLOSUM62 (不训练, 提供进化信息先验)
    - CNN: 学习局部残基基序 (潜在结合口袋特征)
    - [UNK]: 序列缺失时的可学习嵌入
    """

    def __init__(self, embed_dim: int = EMBED_DIM, max_len: int = MAX_LEN):
        super().__init__()
        self.embed_dim = embed_dim
        self.max_len = max_len
        # 固定 BLOSUM 嵌入 (21 类 -> 20 维)
        self.res_embed = nn.Embedding(21, RES_DIM)
        self.res_embed.weight.data.copy_(torch.from_numpy(_B62_MAT))
        self.res_embed.weight.requires_grad = False
        # 1D CNN: (B, RES_DIM, L) -> (B, 64, L) -> (B, 128, L)
        self.conv1 = nn.Conv1d(RES_DIM, 64, kernel_size=5, padding=2)
        self.conv2 = nn.Conv1d(64, 128, kernel_size=5, padding=2)
        self.relu = nn.ReLU()
        self.unk = nn.Parameter(torch.zeros(embed_dim))  # 序列缺失时

    def forward(self, seqs):
        """seqs: list[str] (AA 序列) -> (B, embed_dim) 浮点嵌入。"""
        device = next(self.parameters()).device
        idxs = [aa_to_indices(s) for s in seqs]
        # pad 到同长
        L = max(len(i) for i in idxs)
        L = min(L, self.max_len)
        batch = np.zeros((len(idxs), L), dtype=np.int64)
        for i, ix in enumerate(idxs):
            batch[i, :min(len(ix), L)] = ix[:L]
        x = torch.from_numpy(batch).to(device)            # (B, L)
        e = self.res_embed(x)                             # (B, L, 20)
        e = e.transpose(1, 2)                             # (B, 20, L)
        h = self.relu(self.conv1(e))                      # (B, 64, L)
        h = self.relu(self.conv2(h))                      # (B, 128, L)
        out = torch.max(h, dim=2).values                  # 全局最大池化 (B, 128)
        return out

    def embed_one(self, seq: str) -> np.ndarray:
        """单条序列 -> (embed_dim,) numpy。"""
        if not seq:
            return self.unk.detach().cpu().numpy()
        with torch.no_grad():
            v = self.forward([seq])
        return v.squeeze(0).cpu().numpy()

    def embed_chembl(self, chembl_id: str) -> np.ndarray:
        """chembl_id -> 嵌入 (序列缺失用 UNK)。"""
        seq = chembl_to_seq(chembl_id)
        return self.embed_one(seq) if seq else self.unk.detach().cpu().numpy()

    def save(self, path):
        torch.save(self.state_dict(), path)

    @classmethod
    def load(cls, path, **kw):
        m = cls(**kw)
        verify_asset(path)
        m.load_state_dict(torch.load(path, map_location="cpu"))
        m.eval()
        return m


if __name__ == "__main__":
    enc = TargetSeqEncoder()
    s = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAVGVSIEGGIQEVKPEQ"
    v = enc.embed_one(s)
    print("encoder output dim:", v.shape, "norm=%.3f" % float(np.linalg.norm(v)))
    print("CHEMBL205 seq embed dim:", enc.embed_chembl("CHEMBL205").shape)
    print("UNK (no seq) dim:", enc.embed_chembl("CHEMBL_NONE").shape)

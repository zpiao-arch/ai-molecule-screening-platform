# -*- coding: utf-8 -*-
"""
Layer 2 (升级版): BindingDB 靶点感知结合亲和力预测

取代原 Layer2DeepPurpose (依赖未安装的 DeepPurpose, 运行时返回 docking_normalized=0.0,
导致 L2 槽位(权重0.50)失效 → 闭环"瞎选")。

本模块用本地 BindingDB 已对齐的 5 万条 (药物SMILES, 靶点文本, 结合标签) 训练一个
靶点感知分类器 (MLP / 逻辑回归), 在 MoleculeScorer 的 L2 槽位提供真实可运行的结合打分。

特征定义 (复刻对齐报告 recipe, 776 维, 但 target 哈希改为本模块自洽可复现版本):
  - 分子: RDKit Morgan fingerprint radius=2 fpSize=512  (512 维, 0/1)
  - 分子: 8 个 RDKit 描述符 (原始值): MolWt, MolLogP, TPSA,
          NumHDonors, NumHAcceptors, NumRotatableBonds, NumAromaticRings, QED
  - 靶点: 靶点文本 -> sklearn FeatureHasher(n_features=256, input_type='string') (256 维, 整数计数)
  => 512 + 8 + 256 = 776 维

关键: 训练与推理共用同一个 BindingDBFeature, 保证特征空间一致。
"""
from __future__ import annotations
import os
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

try:
    from .asset_integrity import verify_asset
except ImportError:  # Direct script/PYTHONPATH=scoring compatibility.
    from asset_integrity import verify_asset

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
    from sklearn.feature_extraction import FeatureHasher
    _HAVE_RDKIT = True
except Exception:  # pragma: no cover
    _HAVE_RDKIT = False

MODELS_DIR = Path(__file__).parent / "models" / "bindingdb_l2"
_COMPAT_MODEL_PATH = MODELS_DIR / "l2_model_sklearn_1_7_2.joblib"
MODEL_PATH = _COMPAT_MODEL_PATH if _COMPAT_MODEL_PATH.is_file() else MODELS_DIR / "l2_model.joblib"
PARAMS_PATH = MODELS_DIR / "l2_params.json"


class BindingDBFeature:
    """776 维 (药物SMILES, 靶点文本) 特征, 训练/推理共用。"""

    FP_SIZE = 512
    N_DESC = 8
    N_TARGET = 256
    DIM = FP_SIZE + N_DESC + N_TARGET  # 776

    # 8 个描述符 (顺序固定)
    _DESC_NAMES = [
        "MolWt", "MolLogP", "TPSA",
        "NumHDonors", "NumHAcceptors", "NumRotatableBonds",
        "NumAromaticRings", "QED",
    ]

    def __init__(self):
        if not _HAVE_RDKIT:
            raise RuntimeError("RDKit 不可用, 无法计算 BindingDB L2 特征")
        self._hasher = FeatureHasher(n_features=self.N_TARGET, input_type="string")

    # ---------- 分子特征 ----------
    def mol_features(self, smiles: str) -> Optional[np.ndarray]:
        """返回 520 维分子特征; 解析失败返回 None。"""
        try:
            m = Chem.MolFromSmiles(smiles)
            if m is None:
                return None
            fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(
                m, radius=2, nBits=self.FP_SIZE)
            # ToBitString 返回 "0101..." 字符串 -> 逐字符转 0/1 (512,)
            bits = (np.frombuffer(fp.ToBitString().encode("ascii"),
                                 dtype=np.uint8) - 0x30).astype(np.float32)  # (512,)
            desc = np.array([
                Descriptors.MolWt(m),
                Descriptors.MolLogP(m),
                rdMolDescriptors.CalcTPSA(m),
                rdMolDescriptors.CalcNumHBD(m),
                rdMolDescriptors.CalcNumHBA(m),
                rdMolDescriptors.CalcNumRotatableBonds(m),
                rdMolDescriptors.CalcNumAromaticRings(m),
                Descriptors.qed(m),
            ], dtype=np.float32)
            return np.concatenate([bits, desc])  # (520,)
        except Exception:
            return None

    # ---------- 靶点特征 ----------
    def target_features(self, target_text: str) -> np.ndarray:
        """返回 256 维靶点哈希特征 (确定性)。"""
        toks = [t for t in str(target_text).split() if t]
        if not toks:
            toks = ["<unknown>"]
        X = self._hasher.transform([toks])
        return np.asarray(X.toarray(), dtype=np.float32).ravel()  # (256,)

    # ---------- 组合 ----------
    def features(self, smiles: str, target_text: str) -> Optional[np.ndarray]:
        mf = self.mol_features(smiles)
        if mf is None:
            return None
        tf = self.target_features(target_text)
        return np.concatenate([mf, tf])  # (776,)


def _sigmoid(x):
    return 1.0 / (1.0 + np.exp(-x))


class Layer2BindingDB:
    """MoleculeScorer 的 L2 槽位实现: 给定 (SMILES, 靶点文本) -> P(结合) ∈ [0,1]。

    用法:
        l2 = Layer2BindingDB()
        r = l2.score("CC(=O)Oc1ccccc1C(=O)O",
                     target_text="Carbonic anhydrase 2 Homo sapiens P00918 ...")
        r["docking_normalized"]  # = P(结合)
    """

    def __init__(self, model_path: str = None, params_path: str = None,
                 prefer: str = "mlp"):
        self.model_path = Path(model_path or MODEL_PATH)
        self.params_path = Path(params_path or PARAMS_PATH)
        self.prefer = prefer
        self._feat = BindingDBFeature()
        self._model = None
        self._model_kind = None
        self._meta = None
        self._cache: Dict[tuple, float] = {}

    # ---- 惰性加载模型 ----
    def _ensure_model(self):
        if self._model is not None:
            return
        if not self.model_path.exists():
            raise FileNotFoundError(
                f"未找到 L2 模型权重: {self.model_path}\n"
                f"请先运行训练脚本生成 (见 scientific_validation/.../train_l2_bindingdb.py)")
        try:
            import joblib
        except Exception as e:
            raise RuntimeError(f"joblib 不可用, 无法加载 L2 模型: {e}")
        verify_asset(self.model_path)
        bundle = joblib.load(self.model_path)
        # bundle = {"mlp": clf_or_None, "logreg": clf_or_None, "meta": {...}}
        self._meta = bundle.get("meta", {})
        mlp = bundle.get("mlp")
        logreg = bundle.get("logreg")
        if self.prefer == "mlp" and mlp is not None:
            self._model, self._model_kind = mlp, "mlp"
        elif logreg is not None:
            self._model, self._model_kind = logreg, "logreg"
        elif mlp is not None:
            self._model, self._model_kind = mlp, "mlp"
        else:
            raise RuntimeError("模型权重中 mlp/logreg 均为空")

    @property
    def model_kind(self) -> str:
        self._ensure_model()
        return self._model_kind

    @staticmethod
    def normalize(raw_prob: float) -> float:
        """P(结合) 已是 0-1, 直接裁剪。"""
        return round(float(min(1.0, max(0.0, raw_prob))), 4)

    def score(self, smiles: str, target_text: str = "",
              mol_id: str = "mol") -> Dict:
        key = (smiles, target_text)
        if key in self._cache:
            p = self._cache[key]
        else:
            try:
                self._ensure_model()
                x = self._feat.features(smiles, target_text)
                if x is None:
                    return {"docking_score_kcal_mol": "",
                            "docking_normalized": 0.0,
                            "docking_status": "failed:invalid_smiles",
                            "docking_method": "BindingDB-L2"}
                x = x.reshape(1, -1)
                proba = self._model.predict_proba(x)[0]
                p = float(proba[1]) if proba.shape[0] > 1 else float(proba[0])
                p = self.normalize(p)
                self._cache[key] = p
            except Exception as e:
                return {"docking_score_kcal_mol": "",
                        "docking_normalized": 0.0,
                        "docking_status": f"failed:{e}",
                        "docking_method": "BindingDB-L2"}
        return {"docking_score_kcal_mol": "",
                "docking_normalized": p,
                "docking_status": "success",
                "docking_method": f"BindingDB-L2-{self.model_kind}"}

    def predict_proba(self, X: np.ndarray) -> np.ndarray:
        """批量推理: X 为 (n, 776) 特征矩阵, 返回 (n,) 的 P(结合)。"""
        self._ensure_model()
        X = np.asarray(X, dtype=np.float32)
        proba = self._model.predict_proba(X)
        return proba[:, 1] if proba.shape[1] > 1 else proba[:, 0]


# ===========================================================================
# 升级版 L2: 结构性靶点序列嵌入 (BLOSUM62+CNN) 替代 256 维 target_text 哈希
# 见 target_seq_embedding.py ; 模型用 PyTorch 端到端训练。
# ===========================================================================
class Layer2BindingDBSeq:
    """MoleculeScorer 的 L2 槽位实现 (序列嵌入版)。

    特征: 分子 520 维 (Morgan512+8desc) + 靶点 128 维序列嵌入 (TargetSeqEncoder)
          = 648 维 -> MLP(648->256->64->1) -> sigmoid -> P(结合)

    相比 Layer2BindingDB (256 维 target_text FeatureHash):
      - 靶点是真实蛋白序列的结构/进化表征, 不同靶点嵌入显著可分
      - 直接打在 "弱靶点条件化" 根因上 (见 auc_root_cause_report.md)
    """

    SEQ_EMBED_DIM = 128
    MOL_DIM = 520
    INPUT_DIM = MOL_DIM + SEQ_EMBED_DIM  # 648
    MODEL_DIR = Path(__file__).parent / "models" / "bindingdb_l2_seq"
    MODEL_PATH = MODEL_DIR / "l2_seq.pt"

    def __init__(self, model_path=None, encoder=None):
        self.model_path = Path(model_path or self.MODEL_PATH)
        self.MODEL_DIR.mkdir(parents=True, exist_ok=True)
        self._feat = BindingDBFeature()
        self._encoder = encoder
        self._mlp = None
        self._cache = {}
        self._l2_norm = None  # 训练时记录的输入归一化 (可选)

    # ---- 惰性加载 ----
    def _ensure(self):
        if self._mlp is not None:
            return
        import torch
        if not self.model_path.exists():
            raise FileNotFoundError(
                f"未找到序列版 L2 模型: {self.model_path}\n"
                f"请先运行 train_l2_seq.py")
        if self._encoder is None:
            from target_seq_embedding import TargetSeqEncoder
            self._encoder = TargetSeqEncoder()
        verify_asset(self.model_path)
        bundle = torch.load(self.model_path, map_location="cpu")
        self._mlp = self._build_mlp(tuple(bundle.get("mlp_hidden", (256, 64))))
        self._mlp.load_state_dict(bundle["mlp_state"])
        self._encoder.load_state_dict(bundle["encoder_state"])
        self._mlp.eval()
        self._encoder.eval()

    @staticmethod
    def _build_mlp(hidden=(256, 64)):
        import torch.nn as nn
        layers = [nn.Linear(Layer2BindingDBSeq.INPUT_DIM, hidden[0]), nn.ReLU(),
                  nn.Linear(hidden[0], hidden[1]), nn.ReLU(),
                  nn.Linear(hidden[1], 1)]
        return nn.Sequential(*layers)

    def mol_features(self, smiles):
        return self._feat.mol_features(smiles)

    @staticmethod
    def _chembl_from(text):
        if not text:
            return None
        from target_seq_embedding import extract_chembl_id
        return extract_chembl_id(text)

    def score(self, smiles, chembl_id=None, target_text=None, mol_id="mol"):
        import torch, numpy as np
        key = (smiles, chembl_id, target_text)
        if key in self._cache:
            p = self._cache[key]
        else:
            try:
                self._ensure()
                mf = self.mol_features(smiles)
                if mf is None:
                    return {"docking_score_kcal_mol": "", "docking_normalized": 0.0,
                            "docking_status": "failed:invalid_smiles",
                            "docking_method": "BindingDB-L2-SeqCNN"}
                if chembl_id is None:
                    chembl_id = self._chembl_from(target_text)
                if self._encoder is None:
                    from target_seq_embedding import TargetSeqEncoder
                    self._encoder = TargetSeqEncoder()
                te = self._encoder.embed_chembl(chembl_id or "")
                x = np.concatenate([mf, te]).astype(np.float32)
                with torch.no_grad():
                    t = torch.from_numpy(x).unsqueeze(0)
                    logit = self._mlp(t).item()
                    p = float(1.0 / (1.0 + np.exp(-logit)))
                p = round(float(min(1.0, max(0.0, p))), 4)
                self._cache[key] = p
            except Exception as e:
                return {"docking_score_kcal_mol": "", "docking_normalized": 0.0,
                        "docking_status": f"failed:{e}",
                        "docking_method": "BindingDB-L2-SeqCNN"}
        return {"docking_score_kcal_mol": "", "docking_normalized": p,
                "docking_status": "success", "docking_method": "BindingDB-L2-SeqCNN"}

    def predict_matrix(self, molfeat_matrix, chembl_id, batch_size=4096):
        """批量推理 (供 HTS 全库闭环用)。

        molfeat_matrix: (N, 520) 分子特征 (已预计算缓存)
        chembl_id: 单靶点 chembl_id (该靶点所有分子共享同一靶点嵌入)
        返回 (N,) 的 P(结合)。
        """
        import torch, numpy as np
        self._ensure()
        if self._encoder is None:
            from target_seq_embedding import TargetSeqEncoder
            self._encoder = TargetSeqEncoder()
        te = self._encoder.embed_chembl(chembl_id or "")
        te = np.asarray(te, dtype=np.float32)
        M = np.asarray(molfeat_matrix, dtype=np.float32)
        X = np.concatenate([M, np.tile(te, (M.shape[0], 1))], axis=1)  # (N, 648)
        out = []
        with torch.no_grad():
            for i in range(0, X.shape[0], batch_size):
                xb = torch.from_numpy(X[i:i + batch_size])
                logits = self._mlp(xb).numpy().ravel()
                out.append(1.0 / (1.0 + np.exp(-logits)))
        return np.concatenate(out).astype(np.float32)

# -*- coding: utf-8 -*-
"""
hts_closed_loop 的多进程 worker (独立模块, 避免子进程重跑主脚本顶层代码)。
每个 worker 进程计算单分子 (L1, L3, molfeat向量); 任何异常/卡死返回 None,
由主进程按 35s 超时放弃。
"""
import sys
from pathlib import Path
import numpy as np

SC = Path(__file__).resolve().parent.parent.parent / "评分_work_package" / "评分"


def _worker_l1l3(smi):
    try:
        sys.path.insert(0, str(SC))
        from scoring import Layer1Scorer, Layer3Scorer
        from l2_bindingdb import BindingDBFeature
        feat = BindingDBFeature()
        mf = feat.mol_features(smi)
        if mf is None:
            mf = np.zeros(feat.DIM - feat.N_TARGET, dtype=np.float32)
        l1 = Layer1Scorer(); l3 = Layer3Scorer()
        try:
            v1 = float(l1.score(smi).get("layer1_score") or 0.0)
        except Exception:
            v1 = 0.0
        try:
            a = l3.score(smi).get("admet_score")
            v3 = float(a) if a is not None else 0.5
        except Exception:
            v3 = 0.5
        return (v1, v3, mf)
    except Exception:
        return None

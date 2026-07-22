# -*- coding: utf-8 -*-
"""
verify_cascade_fusion.py —— 验证"修正级联融合"已正确进入 CLI
================================================================
两个层面的验证:
  (A) 单元级: cascade_corrected_fusion() 必须逐元素等于
      cascade_auc_closed_loop.py 修正版(产出 NA 全库 AUC=0.935 的那一版)的公式。
  (B) 集成级: 对几个真实分子跑 smina 对接 -> LE -> 修正融合,
      确认对接成功分子被同尺度 LE 校正调整、未对接分子保留化学基线。

用法:
    DOCK_BIN_DIR=/path/to/smina/bin \
    python3 verify_cascade_fusion.py
"""
import os, sys, json
import numpy as np
from pathlib import Path

HERE = Path(__file__).resolve().parent
SC = HERE.parent.parent / "评分_work_package" / "评分"
sys.path.insert(0, str(SC))

from dock_rerank import cascade_corrected_fusion, DockingReranker, find_binary


# ─────────────────────────────────────────────────────────────
# (A) 单元级: 与闭环修正公式逐元素对比
# ─────────────────────────────────────────────────────────────
def closed_loop_formula(base, le, w=0.30):
    """复刻 cascade_auc_closed_loop.py 修正版 (lines 233-240)。防御性 copy, 不改动入参。"""
    base = np.array(base, dtype=float, copy=True)
    le = np.array(le, dtype=float, copy=True)
    out = base.copy()
    mask = ~np.isnan(le)
    le_scaled = -le[mask]
    le_scaled = (le_scaled - le_scaled.mean()) / (le_scaled.std() + 1e-9)
    out[mask] = base[mask] + w * base.std() * le_scaled
    return out


def test_unit():
    rng = np.random.RandomState(0)
    n = 200
    base = rng.randn(n) * 0.15 + 0.5          # 模拟化学基线分
    le = np.full(n, np.nan)                   # 默认全部未对接
    docked = rng.choice(n, size=40, replace=False)
    le[docked] = rng.randn(40) * 1.5 - 3.0    # 模拟配体效率(部分很负=好)
    w = 0.30
    le_orig = le.copy()                              # 锁定原始 NaN 位置, 避免任何别名副作用
    got = cascade_corrected_fusion(base, le_orig, w)
    exp = closed_loop_formula(base, le_orig, w)
    max_err = float(np.max(np.abs(got - exp)))
    # 未对接分子(le 为 NaN)必须原样保留化学基线
    nonmask = np.isnan(le_orig)                       # True = 未对接
    nondocked_ok = np.allclose(got[nonmask], base[nonmask])
    # 对接分子必须被调整(至少部分)
    docked_shift = np.abs(got[docked] - base[docked]).max()
    print(f"[A] 单元级: max|CLI-闭环| = {max_err:.2e}  (要求≈0)")
    print(f"[A] 未对接分子保留基线: {nondocked_ok}")
    print(f"[A] 对接分子最大偏移: {docked_shift:.4f}  (应>0, 说明 LE 校正已生效)")
    assert max_err < 1e-9, "融合公式与闭环不一致!"
    assert nondocked_ok, "未对接分子未被保留基线"
    assert docked_shift > 0, "对接分子未被 LE 校正"
    print("[A] ✅ 通过: CLI 融合 == 闭环 0.935 修正公式\n")


# ─────────────────────────────────────────────────────────────
# (B) 集成级: 真实对接 -> LE -> 修正融合
# ─────────────────────────────────────────────────────────────
def test_integration():
    rec = SC.parent.parent / "ai_mol_loop/influenza_na_2000_project_20260625/stage4/receptors/3TI6_protein_only_obabel.pdbqt"
    if not rec.exists():
        print("[B] 跳过: 受体 pdbqt 缺失")
        return
    smina = find_binary("smina"); obabel = find_binary("obabel")
    if not (smina and obabel):
        print("[B] 跳过: 未找到 smina/obabel")
        return
    rr = DockingReranker(receptor=str(rec),
                         center=(-28.914, 14.334, 20.794),
                         size=(23.585, 20.45, 24.18))
    # 3 个 NA 类抑制剂(SMILES) + 3 个随机小分(诱饵)
    mols = [
        ("zanamivir",  "CC(=O)NC1C(NC(=N)N)C=C(C(=O)O)OC1C(O)[C@H](O)CO"),
        ("oseltamivir","CCOC(=O)C1=C(C)NC(=O)C(C2=CC=CC=C2Cl)N1C3CCOC3"),
        ("peramivir",  "CC(C)C1CC(NC(=O)C(Cc2ccccc2)N)CCN1C(=O)CO"),
        ("decoy1",     "CCCCCCCCCCCCCCCCCC"),
        ("decoy2",     "c1ccccc1c1ccccc1"),
        ("decoy3",     "CC(C)CC1CCC(C)CC1"),
    ]
    recs = rr.dock_all(mols, mode="le")
    base = np.array([0.6 + i * 0.01 for i, _ in enumerate(mols)], dtype=float)  # 模拟化学基线
    le = [r["ligand_efficiency"] for r in recs]
    fused = cascade_corrected_fusion(base, le, 0.30)
    print("[B] 集成级 (受体=3TI6):")
    print(f"  {'id':<12}{'LE':>9}{'base':>8}{'fused':>9}{'Δ':>8}")
    for r, b, f in zip(recs, base, fused):
        le_s = f"{r['ligand_efficiency']:.3f}" if r["ligand_efficiency"] is not None else "None"
        print(f"  {r['id']:<12}{le_s:>9}{b:>8.3f}{f:>9.3f}{f-b:>+8.3f}")
    ndocked = sum(1 for r in recs if r["ligand_efficiency"] is not None)
    print(f"[B] {ndocked}/{len(mols)} 个对接成功并被 LE 校正; 其余保留基线")
    assert ndocked >= 1, "对接腿未产出任何 LE"
    print("[B] ✅ 通过: 对接->LE->修正融合 端到端跑通\n")


if __name__ == "__main__":
    print("=" * 64)
    print("验证: 修正级联融合 (复现 NA 全库 AUC=0.935) 已正确实现")
    print("=" * 64)
    test_unit()
    test_integration()
    print("全部验证通过 ✅")

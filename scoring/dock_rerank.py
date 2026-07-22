# -*- coding: utf-8 -*-
"""
对接精排模块 (Docking Rerank) — smina(AutoDock Vina 打分) + 配体效率(LE)校正
================================================================================
来自 ② 验证实验的核心发现:
  - 原始 Vina affinity 因「分子越大/越亲脂 → affinity 越负」的系统性偏差而
    **反富集** (affinity vs 重原子数 r = -0.70; 诱饵 -7.90 反而优于活性 -7.55)。
  - 配体效率 LE = affinity / heavy_atom_count 校正后,
    NA 面板富集 AUC 由 0.39 翻转到 **0.756**, recall@30 由 0.06 → 0.688。

因此本模块默认 dock_mode='le' (ligand efficiency); raw affinity 仅作为诊断列输出。

依赖 (离线): 在 `DOCK_BIN_DIR` 指定目录或系统 `PATH` 中提供：
  - smina   (刚性/柔性对接, Vina 打分)
  - obabel  (SMILES -> 3D pdbqt, --gen3d -p 7.4 生理 pH 质子化)

典型用法:
    from dock_rerank import DockingReranker
    rr = DockingReranker(receptor="3TI6.pdbqt",
                         center=(-28.914,14.334,20.794),
                         size=(23.585,20.45,24.18))
    out = rr.dock_all([("mol1","CCO"), ("mol2","c1ccccc1")], mode="le")
    # out: [{id, smiles, affinity, heavy_atoms, ligand_efficiency, dock_rerank_rank, status}, ...]
"""
import hashlib
import os
import sys
import shutil
import subprocess
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import List, Tuple, Optional, Dict

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

RDLogger.DisableLog("rdApp.*")

REPO = Path(os.environ.get("FOUR_LEVEL_ROOT", Path(__file__).resolve().parents[1])).resolve()
_DEFAULT_SMINA_ENV = Path(os.environ.get(
    "DOCK_BIN_DIR",
    REPO / "bin",
)).resolve()


# ─────────────────────────────────────────────────────────────
# 二进制探测
# ─────────────────────────────────────────────────────────────
def find_binary(name: str) -> Optional[str]:
    """优先: <NAME>_BIN / DOCK_BIN_DIR 环境变量; 再已知 micromamba env; 再 PATH。"""
    direct = os.environ.get(f"{name.upper()}_BIN")
    if direct:
        configured = Path(direct).expanduser()
        candidate = configured / name if configured.is_dir() else configured
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate)
    bin_dir = os.environ.get("DOCK_BIN_DIR")
    if bin_dir:
        candidate = Path(bin_dir).expanduser() / name
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate)
    p = _DEFAULT_SMINA_ENV / name
    if p.is_file() and os.access(p, os.X_OK):
        return str(p)
    return shutil.which(name)


@contextmanager
def protonated_receptor(
    receptor: str | Path,
    obabel_bin: str,
    *,
    timeout: int = 300,
):
    """Create a short-lived pH 7.4 receptor and reject stale/partial output."""
    source = Path(receptor).expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError(f"receptor not found: {source}")
    with tempfile.TemporaryDirectory(prefix="four-level-receptor-") as directory:
        output = Path(directory) / "receptor.pdb"
        output.unlink(missing_ok=True)
        try:
            result = subprocess.run(
                [obabel_bin, str(source), "-O", str(output), "-p", "7.4"],
                capture_output=True,
                text=True,
                timeout=timeout,
            )
        except subprocess.TimeoutExpired as exc:
            raise RuntimeError("Open Babel receptor protonation timed out") from exc
        if result.returncode != 0 or not output.is_file() or output.stat().st_size == 0:
            detail = (result.stderr or result.stdout or "no diagnostic output").strip()[:500]
            raise RuntimeError(f"Open Babel receptor protonation failed: {detail}")
        yield str(output)


def _heavy_atoms(smiles: str) -> Optional[int]:
    try:
        m = Chem.MolFromSmiles(smiles)
        return Descriptors.HeavyAtomCount(m) if m else None
    except Exception:
        return None


# ─────────────────────────────────────────────────────────────
# 单步: 配体预处理 + 对接
# ─────────────────────────────────────────────────────────────
def prep_ligand(smiles: str, out_pdbqt: str, obabel_bin: str, pH: float = 7.4,
                timeout: int = 120) -> bool:
    """SMILES -> 3D pdbqt (生理 pH 质子化)。成功返回 True。"""
    output = Path(out_pdbqt)
    output.unlink(missing_ok=True)
    try:
        r = subprocess.run(
            [obabel_bin, f"-:{smiles}", "-O", out_pdbqt, "--gen3d", "-p", str(pH)],
            capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return False
    return r.returncode == 0 and output.is_file() and output.stat().st_size > 0


def dock(smiles_pdbqt: str, receptor: str, center: Tuple[float, float, float],
         size: Tuple[float, float, float], smina_bin: str,
         exhaustiveness: int = 8, cpu: int = 4, num_modes: int = 3,
         seed: int = 42, flexres: str = "", timeout: int = 300) -> Optional[float]:
    """smina 刚性(可选柔性)对接, 返回最优 pose 的 affinity(kcal/mol, 越负越好)。"""
    cmd = [smina_bin, "--receptor", receptor, "--ligand", smiles_pdbqt,
           "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
           "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2]),
           "--exhaustiveness", str(exhaustiveness), "--cpu", str(cpu),
           "--num_modes", str(num_modes), "--seed", str(seed)]
    if flexres:
        cmd += ["--flexres", flexres]
    out_pose = smiles_pdbqt.replace(".pdbqt", "_pose.pdbqt")
    cmd += ["--out", out_pose]
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return None
    best = None
    for line in p.stdout.splitlines():
        s = line.split()
        if len(s) >= 2 and s[0] == "1":      # 第一个 pose 即最优
            try:
                best = float(s[1]); break
            except ValueError:
                pass
    return best


# ─────────────────────────────────────────────────────────────
# 批量对接精排器
# ─────────────────────────────────────────────────────────────
class DockingReranker:
    def __init__(self, receptor: str,
                 center: Tuple[float, float, float] = (-28.914, 14.334, 20.794),
                 size: Tuple[float, float, float] = (23.585, 20.45, 24.18),
                 smina_bin: Optional[str] = None, obabel_bin: Optional[str] = None,
                 pH: float = 7.4, exhaustiveness: int = 8, cpu: int = 4,
                 num_modes: int = 3, seed: int = 42, flexres: str = "",
                 workdir: Optional[str] = None):
        self.receptor = receptor
        self.center = center
        self.size = size
        self.smina_bin = smina_bin or find_binary("smina")
        self.obabel_bin = obabel_bin or find_binary("obabel")
        self.pH = pH
        self.exhaustiveness = exhaustiveness
        self.cpu = cpu
        self.num_modes = num_modes
        self.seed = seed
        self.flexres = flexres
        self.workdir = workdir or tempfile.mkdtemp(prefix="dock_rerank_")
        os.makedirs(self.workdir, exist_ok=True)
        if not (self.smina_bin and self.obabel_bin):
            raise RuntimeError(
                "未找到 smina/obabel 二进制; 请设置 SMINA_BIN/OBABEL_BIN 或 DOCK_BIN_DIR 环境变量")

    @property
    def available(self) -> bool:
        return bool(self.smina_bin and self.obabel_bin and os.path.exists(self.receptor))

    def dock_one(self, mol_id: str, smiles: str) -> Dict:
        token = hashlib.sha256(f"{mol_id}\0{smiles}".encode("utf-8")).hexdigest()[:20]
        lig = str(Path(self.workdir) / f"ligand_{token}.pdbqt")
        rec = {"id": mol_id, "smiles": smiles,
               "affinity": None, "heavy_atoms": _heavy_atoms(smiles),
               "ligand_efficiency": None, "status": "ok"}
        try:
            if prep_ligand(smiles, lig, self.obabel_bin, self.pH):
                aff = dock(lig, self.receptor, self.center, self.size, self.smina_bin,
                           self.exhaustiveness, self.cpu, self.num_modes, self.seed, self.flexres)
                if aff is None:
                    rec["status"] = "dock_no_score"
                else:
                    rec["affinity"] = round(aff, 3)
            else:
                rec["status"] = "prep_failed"
        except Exception as e:
            rec["status"] = f"err:{type(e).__name__}"
        if rec["affinity"] is not None and rec["heavy_atoms"]:
            rec["ligand_efficiency"] = round(rec["affinity"] / rec["heavy_atoms"], 4)
        return rec

    def dock_all(self, molecules: List[Tuple[str, str]], mode: str = "le",
                 progress_every: int = 10, log_fn=print) -> List[Dict]:
        """对接全部分子并计算 LE, 返回带 dock_rerank_rank 的记录列表。

        mode='le'  -> 按配体效率(越负越好)重排 (默认, 校正尺寸偏差)
        mode='raw' -> 按原始 affinity(越负越好)重排
        dock_rerank_rank: 1 = 最佳。
        """
        if not self.available:
            raise RuntimeError("DockingReranker 不可用: 缺少受体/二进制")
        out = []
        for k, (mid, smi) in enumerate(molecules, 1):
            rec = self.dock_one(mid, smi)
            out.append(rec)
            if log_fn and (k % progress_every == 0 or k == len(molecules)):
                ok = sum(1 for r in out if r["affinity"] is not None)
                log_fn(f"[dock] {k}/{len(molecules)} 成功对接 {ok}")
        # 重排
        key = (lambda r: r["ligand_efficiency"]) if mode == "le" else (lambda r: r["affinity"])
        scored = [r for r in out if key(r) is not None]
        scored.sort(key=key)              # 升序: 最负(最好)在前
        for rank, r in enumerate(scored, 1):
            r["dock_rerank_rank"] = rank
        for r in out:
            r.setdefault("dock_rerank_rank", None)
        return out


# ─────────────────────────────────────────────────────────────
# 便捷函数: 直接对一组分子精排
# ─────────────────────────────────────────────────────────────
def rerank(molecules: List[Tuple[str, str]], receptor: str,
           center: Tuple[float, float, float], size: Tuple[float, float, float],
           mode: str = "le", **kwargs) -> List[Dict]:
    rr = DockingReranker(receptor=receptor, center=center, size=size, **kwargs)
    return rr.dock_all(molecules, mode=mode)


# ─────────────────────────────────────────────────────────────
# 历史级联展示使用的修正融合公式
# ─────────────────────────────────────────────────────────────
def cascade_corrected_fusion(base_scores, le_values, w: float = 0.30):
    """
    修正版级联融合: 以化学基线分为全库基线, 仅对【对接成功】分子叠加同尺度 LE 校正。

    公式来自历史 cascade_auc_closed_loop.py 的人工平衡集展示版本:
        score = base + w · σ(base) · z(−LE)      # 对接成功分子
        score = base                              # 未对接分子, 保留化学基线
    其中 z(−LE) = (−LE − mean) / (std + ε), LE 越负越好故取负; σ(base) 用全库基线标准差,
    使 LE 校正与基线同尺度, 从而【不破坏】化学分的整体排序。

    它修复了“把非对接分子硬设为极小值”的排序破坏，但历史 0.935 只属于
    人工平衡/特定展示集，不能据此声称真实万级池或跨靶点泛化达到 0.935。

    参数:
        base_scores: 与 le_values 对齐的化学基线分 (list/ndarray)
        le_values:   配体效率 LE (=affinity/重原子, 越负越好), 与 base 对齐; 未对接用 None/NaN
        w:           LE 校正权重 (默认 0.30, 对应 W_LE)
    返回 ndarray: 融合后的最终分 (越高越好)
    """
    base = np.array(base_scores, dtype=float, copy=True)            # 防御性 copy, 绝不改动入参
    le = np.array([(np.nan if v is None else float(v)) for v in le_values],
                  dtype=float, copy=True)
    out = base.copy()
    mask = ~np.isnan(le)
    if mask.any():
        lez = -le[mask]                                  # 越负越好 -> 取负
        lez = (lez - lez.mean()) / (lez.std() + 1e-9)
        out[mask] = base[mask] + w * (base.std() + 1e-9) * lez
    return out


if __name__ == "__main__":
    # 自测: 一个已知 NA 抑制剂
    demo = [("zanamivir_like", "CC(=O)NC1C(NC(=N)N)C=C(C(=O)O)OC1C(O)[C@H](O)CO")]
    rec = Path(os.environ.get(
        "FOUR_LEVEL_RECEPTOR",
        Path(__file__).resolve().parent / "receptors" / "3TI6_protein_only_obabel.pdbqt",
    ))
    res = rerank(demo, str(rec), (-28.914, 14.334, 20.794), (23.585, 20.45, 24.18))
    print(res[0])

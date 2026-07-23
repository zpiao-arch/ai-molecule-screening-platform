"""
AI药物分子打分系统
四层评分: 分子质量(L1) → 结合亲和力(L2) → ADMET安全性(L3) → 综合排名(L4)

用法:
    # 作为脚本运行 (默认靶点HIV-1蛋白酶)
    python scoring.py -i candidates.csv -o scores.csv

    # 自定义靶点
    python scoring.py -i candidates.csv -t EGFR

    # 作为模块导入
    from scoring import MoleculeScorer
    scorer = MoleculeScorer()
    result = scorer.score_one("CC(=O)Oc1ccccc1C(=O)O")
"""

import argparse
import csv
import math
import os
import subprocess
import tempfile
from contextlib import nullcontext
from pathlib import Path

from typing import Dict, List, Optional, Tuple

import numpy as np
from rdkit import Chem
from rdkit.Chem import (
    Descriptors,
    Lipinski,
    QED,
    rdFingerprintGenerator,
)
from rdkit.DataStructs import TanimotoSimilarity

try:
    from .asset_integrity import verify_asset
except ImportError:  # Direct script/PYTHONPATH=scoring compatibility.
    from asset_integrity import verify_asset

try:
    from .asset_paths import resolve_asset_paths
except ImportError:  # Direct script/PYTHONPATH=scoring compatibility.
    from asset_paths import resolve_asset_paths


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 工具函数
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def clamp(x: float, lo: float = 0.0, hi: float = 1.0) -> float:
    return max(lo, min(hi, x))


def tanimoto(smi_a: str, smi_b: str) -> float:
    mol_a, mol_b = Chem.MolFromSmiles(smi_a), Chem.MolFromSmiles(smi_b)
    if mol_a is None or mol_b is None:
        return 0.0
    fpg = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    return TanimotoSimilarity(fpg.GetFingerprint(mol_a), fpg.GetFingerprint(mol_b))


def select_cascade_candidates(molecules: List[Tuple[str, str]],
                              l2_scores: Dict[str, float],
                              top_n: int) -> List[Tuple[str, str]]:
    if top_n <= 0:
        return []
    indexed = list(enumerate(molecules))
    indexed.sort(key=lambda item: (-float(l2_scores.get(item[1][0], float("-inf"))), item[0]))
    return [molecule for _, molecule in indexed[:min(top_n, len(indexed))]]


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Layer 1: 分子质量评分
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# 毒性预警 SMARTS
TOXIC_SMARTS = {
    # ── 高优先级: DNA反应性/致癌 ──
    "亚硝胺(N-亚硝基,强致癌)":   "[#7]-[#7]=[#8]",
    "氮芥(烷化剂)":              "N(CCCl)CCCl",
    "硫芥(糜烂性毒剂)":          "SCCCl",
    "环氧化物(烷化)":            "C1OC1",
    "硫酸二甲酯(强烷化剂)":      "COS(=O)(=O)OC",
    "氯甲酸酯(烷化剂)":          "ClC(=O)O[#6]",
    "卤代乙酰(催泪/烷化)":       "[F,Cl,Br,I]CC(=O)",
    "磺酸酯(潜在烷化剂)":        "[SX4](=O)(=O)O[#6]",
    "氮丙啶(高活性烷化剂)":      "C1CN1",
    "β-丙内酯(致癌)":           "O=C1CCO1",
    "卤代醚(致癌,BCME类)":       "[Cl,Br]COC[Cl,Br]",
    "炔丙基卤(高反应性)":        "C#CC[Cl,Br,I]",
    "烯丙基卤(反应性)":          "C=CC[Cl,Br,I]",
    # ── 高优先级: 神经毒/急性毒 ──
    "有机磷酸酯(神经毒)":        "[P](=O)([O,N])([F,Cl,Br,I])[#6]",
    "叠氮基(爆炸性/毒)":         "[#7-]=[#7+]=[#7-]",
    "联吡啶季铵盐(百草枯类)":    "[n+]1ccc(cc1)-c2cc[n+]cc2",
    # ── 中优先级: 器官毒性/致癌 ──
    "Michael受体(αβ不饱和酮)":  "C=C[CX3](=O)",
    "甲醛/活性醛":               "[CX3H2](=O)",
    "醛基_反应性":               "[CX3H1](=O)[#6]",
    "异氰酸酯(呼吸致敏)":        "[NX2]=C=[OX1]",
    "异硫氰酸酯(致敏)":          "[NX2]=C=[SX1]",
    "丙烯酰胺(神经毒)":           "[CH2]=[CH]C(=O)N",
    "醌类(氧化应激)":            "O=c1[cX3]=c([cX3])[cX3]=c([cX3])c1=O",
    "硫脲(肝毒)":                "NC(=S)N",
    "硫代乙酰胺(肝毒)":          "CC(=S)N",
    "肼类(肝毒/致癌)":           "[#7X3H2][#7X3H1,#7X3H2]",
    "芳基肼(致癌)":              "[NX3]([NX3])[c]",
    "苯胺_卤代(致癌)":           "Nc1ccc([Cl,Br])cc1",
    "苯胺_N-烷基(致癌)":         "[#6]N([#6])c1ccccc1",
    "硝基芳香(致突变)":          "[NX3](=O)=Oc1ccccc1",
    "偶氮基(潜在致癌)":          "[c]~[NX2]=[NX2]~[c]",
    # ── 中优先级: 发育/生殖毒 ──
    "邻苯二甲酰亚胺(致畸)":       "O=C1C2=C(C=CC=C2)C(=O)N1",
    # ── 低优先级: 环境持久性/慢毒 ──
    "多环芳烃_4环+(致癌)":       "c1ccc2c(c1)ccc1c2cccc1",
    "多氯联苯(环境毒素)":        "Clc1c(Cl)c(Cl)c(Cl)c(Cl)c1",
    "多氯二噁英/呋喃(TCDD类)":   "Clc1cc2Oc3cc(Cl)c(Cl)cc3Oc2cc1Cl",
    "多卤代烃_3+":               "[CX4]([Cl,Br,I])([Cl,Br,I])[Cl,Br,I]",
    "磺酰卤(反应性)":            "[SX4](=O)(=O)[Cl,Br]",
    "酸酐(反应性)":              "[CX3](=O)O[CX3](=O)",
}
# 严重毒性基团: 命中任意一条 → L3上限0.35
TOXIC_SEVERE = {
    "亚硝胺(N-亚硝基,强致癌)", "氮芥(烷化剂)", "硫芥(糜烂性毒剂)",
    "硫酸二甲酯(强烷化剂)", "氯甲酸酯(烷化剂)",
    "卤代醚(致癌,BCME类)", "有机磷酸酯(神经毒)", "叠氮基(爆炸性/毒)",
    "多氯二噁英/呋喃(TCDD类)", "甲醛/活性醛",
    "联吡啶季铵盐(百草枯类)", "无机氰化物(急性剧毒)",
}
# 重金属原子: 命中任意一个 → L3上限0.25
HEAVY_METALS = {80, 82, 33, 48, 50, 51, 81}  # Hg, Pb, As, Cd, Sn, Sb, Tl
# 碱金属: 与C≡N组合 → 无机氰化物
ALKALI_METALS = {3, 11, 19, 37, 55}  # Li, Na, K, Rb, Cs


class Layer1Scorer:
    """分子质量评分: QED / SA / Lipinski / MW / logP / TPSA / HBD / HBA"""

    WEIGHTS = {
        "qed": 0.20, "sa": 0.15, "lipinski": 0.15,
        "mw": 0.12, "logp": 0.12, "tpsa": 0.10,
        "hbd": 0.08, "hba": 0.08,
    }

    @staticmethod
    def properties(smiles: str) -> Dict:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {"valid": 0, "error": "invalid_smiles"}

        mw      = Descriptors.MolWt(mol)
        logp    = Descriptors.MolLogP(mol)
        tpsa    = Descriptors.TPSA(mol)
        hbd     = Lipinski.NumHDonors(mol)
        hba     = Lipinski.NumHAcceptors(mol)
        rot     = Lipinski.NumRotatableBonds(mol)
        heavy   = mol.GetNumHeavyAtoms()
        rings   = Lipinski.RingCount(mol)
        arom    = Lipinski.NumAromaticRings(mol)
        charge  = Chem.GetFormalCharge(mol)
        csp3    = Descriptors.FractionCSP3(mol)
        qed     = QED.qed(mol)
        bertz   = Descriptors.BertzCT(mol)
        lv      = sum([mw > 500, logp > 5, hbd > 5, hba > 10])

        # SA (片段计数法)
        hetero = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() not in (6, 1))
        halogen = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() in (9, 17, 35, 53))
        sa_pen = (max(0, rings - 3) * 0.15 + max(0, arom - 2) * 0.10
                  + max(0, rot - 6) * 0.08 + max(0, hetero - 6) * 0.06
                  + halogen * 0.05 + max(0, heavy - 35) * 0.02)
        sa = clamp(0.95 - sa_pen)

        return {
            "valid": 1, "mw": round(mw, 2), "logp": round(logp, 3),
            "tpsa": round(tpsa, 2), "hbd": hbd, "hba": hba,
            "rotatable_bonds": rot, "heavy_atoms": heavy,
            "rings": rings, "aromatic_rings": arom,
            "formal_charge": charge, "fraction_csp3": round(csp3, 4),
            "lipinski_violations": lv, "qed": round(qed, 4),
            "sa": round(sa, 4), "bertz": round(bertz, 2),
        }

    @classmethod
    def score(cls, smiles: str, all_smiles: Optional[List[str]] = None) -> Dict:
        props = cls.properties(smiles)
        if not props.get("valid"):
            return {**props, "layer1_score": 0.0, "diversity": 0.0}

        sub = {}
        sub["qed"]     = props["qed"]
        sub["sa"]      = props["sa"]
        sub["lipinski"] = clamp(1.0 - props["lipinski_violations"] / 4.0)

        mw = props["mw"]
        # MW: 150-600满分(药物常见区间), 600-800缓降, <150逐步降
        sub["mw"] = (0.1 if mw < 50 else 0.5 + 0.5 * (mw - 50) / 100 if mw < 150
                     else 1.0 if mw <= 600 else 1.0 - (mw - 600) / 200 * 0.4 if mw <= 800 else 0.2)

        lp = props["logp"]
        # logP: -0.5~4满分, 4~6缓降到0.5, >6继续降
        sub["logp"] = (0.3 if lp < -2 else 0.5 + 0.5 * (lp + 2) / 1.5 if lp < -0.5
                       else 1.0 if lp <= 4 else 1.0 - (lp - 4) / 2 * 0.4 if lp <= 6 else 0.3)

        tp = props["tpsa"]
        # TPSA: 20-140满分, 140-200缓降
        sub["tpsa"] = (0.4 if tp < 20 else 1.0 if tp <= 140
                       else 1.0 - (tp - 140) / 60 * 0.5 if tp <= 200 else 0.3)

        sub["hbd"] = clamp(1.0 - max(0, props["hbd"] - 4) / 5.0)
        # HBA: 1-10满分, 10-16缓降
        sub["hba"] = (0.3 if props["hba"] < 1 else 1.0 if props["hba"] <= 10
                      else clamp(1.0 - (props["hba"] - 10) / 6.0))

        l1 = clamp(sum(sub[k] * cls.WEIGHTS[k] for k in cls.WEIGHTS))

        # 多样性
        div = 0.0
        if all_smiles and len(all_smiles) > 1:
            sims = [tanimoto(smiles, o) for o in all_smiles if o != smiles]
            div = clamp(1.0 - sum(sims) / len(sims)) if sims else 1.0

        return {**props, "layer1_score": round(l1, 4), "diversity": round(div, 4)}


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Layer 2: DeepPurpose 结合亲和力预测 (替代分子对接)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━



# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Layer 3: ADMET预测
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Layer 3: ADMET预测 (DeepChem模型 + RDKit规则混合)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

class Layer3Scorer:
    """ADMET预测: DeepChem训练模型 (Tox21/BBBP/ClinTox/SIDER) + RDKit规则补充"""

    REQUIRED_MODELS = ("tox21", "bbbp", "clintox", "sider")

    def __init__(self, *, strict_backends: bool = False, model_dir: Optional[str | Path] = None):
        self._models = None
        self.strict_backends = bool(strict_backends)
        self.model_dir = Path(model_dir) if model_dir is not None else Path(__file__).parent / "models" / "admet"
        if self.strict_backends:
            self._models = self._load_models()

    @property
    def models(self):
        if self._models is None:
            self._models = self._load_models()
        return self._models

    def _load_models(self) -> Dict:
        import pickle
        models = {}
        failures = []
        for name in self.REQUIRED_MODELS:
            path = self.model_dir / f"{name}.pkl"
            if not path.is_file():
                failures.append(f"missing:{name}")
                continue
            try:
                verify_asset(path)
                with path.open("rb") as f:
                    models[name] = pickle.load(f)
            except Exception as exc:
                failures.append(f"load:{name}:{type(exc).__name__}")
        if self.strict_backends and failures:
            raise RuntimeError(f"L3 ADMET backend unavailable: {', '.join(failures)}")
        return models

    @staticmethod
    def _fp(smiles: str) -> "np.ndarray":
        fpg = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return np.zeros(1024, dtype=np.float32)
        return np.array(fpg.GetFingerprintAsNumPy(mol), dtype=np.float32).reshape(1, -1)

    def predict_tox21(self, smiles: str) -> Dict:
        """Tox21: 12个核受体+应激反应毒性任务"""
        m = self.models.get("tox21")
        if m is None:
            return {}
        fp = self._fp(smiles)
        probs = m["model"].predict_proba(fp) if hasattr(m["model"], "predict_proba") else None
        if probs is None:
            return {}
        # 提取每个任务的毒性概率(类别1)
        results = {}
        for i, task in enumerate(m["meta"]["tasks"]):
            if isinstance(probs, list) and i < len(probs):
                prob = probs[i][0][1] if probs[i].shape[1] > 1 else probs[i][0][0]
                results[task] = round(float(prob), 4)
        return results

    def predict_bbbp(self, smiles: str) -> float:
        """BBB穿透性"""
        m = self.models.get("bbbp")
        if m is None:
            return 0.5
        fp = self._fp(smiles)
        prob = m["model"].predict_proba(fp)[0][1]
        return round(float(prob), 4)

    def predict_clintox(self, smiles: str) -> Dict:
        """临床毒性: FDA批准状态 + 临床试验毒性"""
        m = self.models.get("clintox")
        if m is None:
            return {}
        fp = self._fp(smiles)
        probs = m["model"].predict_proba(fp)
        results = {}
        for i, task in enumerate(m["meta"]["tasks"]):
            results[task] = round(float(probs[i][0][1]), 4) if probs[i].shape[1] > 1 else 0.0
        return results

    def predict_sider(self, smiles: str) -> Dict:
        """SIDER: 27个副作用类别"""
        m = self.models.get("sider")
        if m is None:
            return {}
        fp = self._fp(smiles)
        probs = m["model"].predict_proba(fp)
        results = {}
        for i, task in enumerate(m["meta"]["tasks"]):
            prob = probs[i][0][1] if probs[i].shape[1] > 1 else probs[i][0][0]
            if prob > 0.8:  # 只返回高置信度副作用
                results[task] = round(float(prob), 4)
        return results

    def score(self, smiles: str, _predictions=None) -> Dict:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {"admet_score": 0.0, "toxicity_count": 0, "toxicity_flags": "invalid"}

        # ── DeepChem模型预测 ──
        # Tox21毒性
        if _predictions is None:
            tox21 = self.predict_tox21(smiles)
            bbb_p = self.predict_bbbp(smiles)
            clintox = self.predict_clintox(smiles)
            sider = self.predict_sider(smiles)
        else:
            tox21, bbb_p, clintox, sider = _predictions
        tox21_risk = sum(tox21.values()) / max(len(tox21), 1) if tox21 else 0.5

        # BBB穿透
        ct_tox = clintox.get("CT_TOX", 0.5)

        # SIDER副作用
        sider_risk = sum(sider.values()) / max(len(sider), 1) if sider else 0.0

        # ── RDKit规则补充 ──
        # 溶解度 (ESOL)
        lp, mw, tp = Descriptors.MolLogP(mol), Descriptors.MolWt(mol), Descriptors.TPSA(mol)
        sol = 0.16 - 0.63 * lp - 0.0062 * mw + 0.066 * tp / 100

        # hERG风险 (规则)
        arom, rot = Lipinski.NumAromaticRings(mol), Lipinski.NumRotatableBonds(mol)
        herg = clamp(0.15 * (mw > 350) + 0.1 * (mw > 500) + 0.15 * (lp > 3)
                     + 0.1 * (lp > 5) + 0.15 * (arom >= 2) + 0.1 * (rot > 5))

        # SMARTS毒性基团
        tox_smarts = {}
        for name, smarts in TOXIC_SMARTS.items():
            pat = Chem.MolFromSmarts(smarts)
            if pat and mol.HasSubstructMatch(pat):
                tox_smarts[name] = len(mol.GetSubstructMatches(pat))

        # CYP风险
        hba_val = Lipinski.NumHAcceptors(mol)
        cyp3a4 = clamp(0.2 * (mw > 400) + 0.2 * (lp > 3)
                       + 0.2 * (arom >= 3) + 0.1 * (hba_val > 5))
        cyp2d6 = clamp(0.2 * (lp > 2.5) + 0.2 * (arom >= 2)
                       + 0.1 * (Descriptors.NumAliphaticRings(mol) > 0))

        # ── 重金属/碱金属检测 ──
        has_heavy_metal = any(a.GetAtomicNum() in HEAVY_METALS for a in mol.GetAtoms())
        has_alkali = any(a.GetAtomicNum() in ALKALI_METALS for a in mol.GetAtoms())

        # ── 特殊检测: 无机氰化物 (C≡N + 小分子或碱金属) ──
        cn_mol = Chem.MolFromSmarts("[#6]#[#7]")
        if cn_mol and mol.HasSubstructMatch(cn_mol) and (mw < 100 or has_alkali):
            has_inorganic_cyanide = True
        else:
            has_inorganic_cyanide = False

        # ── SMARTS毒性分级 ──

        # ── SMARTS毒性分级 ──
        n_severe = sum(1 for k in tox_smarts if k in TOXIC_SEVERE)
        n_moderate = len(tox_smarts) - n_severe
        # 无机氰化物等同于严重基团
        if has_inorganic_cyanide:
            n_severe += 1
            tox_smarts["无机氰化物(急性剧毒)"] = 1
        # 严重基团: 每条0.30 | 中等基团: 每条0.22 | 上限0.85
        smarts_penalty = clamp(n_severe * 0.30 + n_moderate * 0.22, 0.0, 0.85)

        # ── 聚合评分 ──
        # 毒性 (SMARTS硬证据权重50% + DeepChem软预测辅助)
        total_tox = (smarts_penalty * 0.50       # 结构警报(硬证据,最高权)
                     + tox21_risk * 0.20          # 细胞毒性(软预测)
                     + ct_tox * 0.15              # 临床毒性(软预测)
                     + sider_risk * 0.15)         # 副作用(软预测)

        # BBB (DeepChem BBBP模型)
        bbb_score = 1.0 - bbb_p if bbb_p < 0.5 else 0.3

        # 安全性 (ClinTox + SIDER + hERG)
        safety = clamp(1.0 - (ct_tox * 0.4
                              + sider_risk * 0.3
                              + herg * 0.3))

        # 吸收性 (溶解度 + Lipinski/Veber)
        sol_norm = clamp((sol + 6) / 7)
        vio = sum([mw > 500, lp > 5, Lipinski.NumHDonors(mol) > 5, hba_val > 10])
        oral = clamp(1.0 - vio * 0.15 - 0.1 * (rot > 10) - 0.15 * (tp > 140))

        # 综合ADMET (毒性安全权重提升至40%)
        parts = {
            "toxicity_safety": 1.0 - total_tox,
            "bbb_safety": bbb_score,
            "organ_safety": safety,
            "oral_absorption": oral,
            "solubility": sol_norm,
        }
        ws = {"toxicity_safety": 0.40, "organ_safety": 0.25,
              "oral_absorption": 0.15, "solubility": 0.10, "bbb_safety": 0.10}
        admet_raw = sum(parts[k] * ws.get(k, 0) for k in parts)
        admet_raw = clamp(admet_raw)

        # ── 红牌机制: 严重毒性基团 → L3强制天花板 ──
        if has_heavy_metal:
            admet = min(admet_raw, 0.25)
        elif n_severe >= 2:
            admet = min(admet_raw, 0.25)
        elif n_severe >= 1:
            admet = min(admet_raw, 0.35)
        elif n_moderate >= 2:
            admet = min(admet_raw, 0.50)
        else:
            admet = admet_raw

        # 构建毒性标志描述
        tox_flags = []
        if has_heavy_metal:
            tox_flags.append("HEAVY_METAL(高毒)")
        if n_severe > 0:
            severe_names = [k for k in tox_smarts if k in TOXIC_SEVERE]
            tox_flags.append(f"SEVERE: {';'.join(severe_names)}")
        if n_moderate > 0:
            moderate_names = [k for k in tox_smarts if k not in TOXIC_SEVERE]
            tox_flags.append(f"SMARTS: {';'.join(moderate_names)}")
        if sider:
            tox_flags.append(f"SIDER: {';'.join(sider.keys())}")
        if ct_tox > 0.5:
            tox_flags.append(f"ClinTox={ct_tox:.2f}")

        return {
            "admet_score": round(admet, 4),
            "toxicity_count": len(tox_smarts) + len(sider) + (1 if has_heavy_metal else 0),
            "toxicity_severe": n_severe,
            "toxicity_flags": " | ".join(tox_flags) if tox_flags else "low",
            "solubility_logS": round(sol, 2),
            "bbb_prob": bbb_p,
            "cyp3a4_risk": round(cyp3a4, 2),
            "cyp2d6_risk": round(cyp2d6, 2),
            "herg_risk": round(herg, 2),
            "hepatotoxicity_risk": round(ct_tox, 2),
            "oral_bioavailability": round(oral, 2),
            "tox21_risk": round(tox21_risk, 2),
            "sider_risk": round(sider_risk, 2),
        }

    def score_many(self, smiles_list: List[str]) -> List[Dict]:
        """Batch the four sklearn model forwards, then reuse the scalar rules exactly."""
        smiles = list(smiles_list)
        if not smiles:
            return []
        molecules = [Chem.MolFromSmiles(smi) for smi in smiles]
        valid_indices = [index for index, mol in enumerate(molecules) if mol is not None]
        result = [
            {"admet_score": 0.0, "toxicity_count": 0, "toxicity_flags": "invalid"}
            if mol is None else None
            for mol in molecules
        ]
        if not valid_indices:
            return result

        fps = np.vstack([self._fp(smiles[index]) for index in valid_indices])
        models = self.models

        tox21_rows = [{} for _ in valid_indices]
        tox_model = models.get("tox21")
        if tox_model is not None:
            outputs = tox_model["model"].predict_proba(fps)
            for task_index, task in enumerate(tox_model["meta"]["tasks"]):
                values = outputs[task_index]
                probabilities = values[:, 1] if values.shape[1] > 1 else values[:, 0]
                for row_index, probability in enumerate(probabilities):
                    tox21_rows[row_index][task] = round(float(probability), 4)

        bbb_rows = [0.5] * len(valid_indices)
        bbb_model = models.get("bbbp")
        if bbb_model is not None:
            values = bbb_model["model"].predict_proba(fps)
            probabilities = values[:, 1] if values.shape[1] > 1 else values[:, 0]
            bbb_rows = [round(float(value), 4) for value in probabilities]

        clintox_rows = [{} for _ in valid_indices]
        clintox_model = models.get("clintox")
        if clintox_model is not None:
            outputs = clintox_model["model"].predict_proba(fps)
            for task_index, task in enumerate(clintox_model["meta"]["tasks"]):
                values = outputs[task_index]
                probabilities = values[:, 1] if values.shape[1] > 1 else values[:, 0]
                for row_index, probability in enumerate(probabilities):
                    clintox_rows[row_index][task] = round(float(probability), 4) if values.shape[1] > 1 else 0.0

        sider_rows = [{} for _ in valid_indices]
        sider_model = models.get("sider")
        if sider_model is not None:
            outputs = sider_model["model"].predict_proba(fps)
            for task_index, task in enumerate(sider_model["meta"]["tasks"]):
                values = outputs[task_index]
                probabilities = values[:, 1] if values.shape[1] > 1 else values[:, 0]
                for row_index, probability in enumerate(probabilities):
                    if probability > 0.8:
                        sider_rows[row_index][task] = round(float(probability), 4)

        for row_index, original_index in enumerate(valid_indices):
            result[original_index] = self.score(
                smiles[original_index],
                _predictions=(tox21_rows[row_index], bbb_rows[row_index], clintox_rows[row_index], sider_rows[row_index]),
            )
        return result


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Layer 4: 综合评分 + 虚拟筛选验证
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

class Layer4Aggregator:
    """综合打分 + EF/AUC/BEDROC"""

    WEIGHTS = {"layer1": 0.25, "docking": 0.40, "admet": 0.35}

    @staticmethod
    def combine(l1, l2, l3, l4):
        values = (0.20 * np.asarray(l1) + 0.50 * np.asarray(l2)
                  + 0.20 * np.asarray(l3) + 0.10 * np.asarray(l4))
        rounded = np.round(values, 4)
        return float(rounded) if rounded.ndim == 0 else rounded

    @classmethod
    def final_score(cls, l1: float, dock_norm: float, admet: float) -> float:
        return round(l1 * cls.WEIGHTS["layer1"]
                     + dock_norm * cls.WEIGHTS["docking"]
                     + admet * cls.WEIGHTS["admet"], 4)

    @staticmethod
    def quality_gate(row: Dict) -> Tuple[str, str]:
        reasons = []
        if float(row.get("docking_normalized", 0) or 0) < 0.2:
            reasons.append("L2结合弱")
        if float(row.get("admet_score", 0) or 0) < 0.5:
            reasons.append("ADMET风险高")
        if int(row.get("toxicity_count", 0) or 0) >= 2:
            reasons.append("多个毒性基团")
        if int(row.get("toxicity_severe", 0) or 0) >= 1:
            reasons.append("严重毒性基团")
        if float(row.get("oral_bioavailability", 0) or 0) < 0.4:
            reasons.append("口服利用度低")
        if float(row.get("herg_risk", 0) or 0) > 0.6:
            reasons.append("hERG风险高")
        return ("PASS", "") if len(reasons) < 2 else ("FAIL", "; ".join(reasons))

    @staticmethod
    def roc_auc(pos: List[float], neg: List[float], *, higher_is_better: bool = True) -> float:
        if not pos or not neg:
            return 0.0
        wins = sum(
            1.0 if (p > n if higher_is_better else p < n) else 0.5 if p == n else 0.0
            for p in pos
            for n in neg
        )
        return round(wins / (len(pos) * len(neg)), 4)

    @staticmethod
    def enrichment_factor(
        pos: List[float],
        neg: List[float],
        frac: float = 0.05,
        *,
        higher_is_better: bool = True,
    ) -> float:
        ranked = sorted(
            [(float(score), 1) for score in pos] + [(float(score), 0) for score in neg],
            key=lambda item: item[0],
            reverse=higher_is_better,
        )
        n, n_pos = len(ranked), len(pos)
        if n == 0 or n_pos == 0:
            return 0.0
        k = max(1, int(n * frac))
        boundary = ranked[k - 1][0]
        if higher_is_better:
            better = [item for item in ranked if item[0] > boundary]
        else:
            better = [item for item in ranked if item[0] < boundary]
        tied = [item for item in ranked if item[0] == boundary]
        remaining = k - len(better)
        hits = sum(label for _, label in better)
        if tied and remaining > 0:
            hits += remaining * sum(label for _, label in tied) / len(tied)
        return round((hits / k) / (n_pos / n), 4) if n_pos / n > 0 else 0.0

    @staticmethod
    def bedroc(
        pos: List[float],
        neg: List[float],
        alpha: float = 20.0,
        *,
        higher_is_better: bool = True,
    ) -> float:
        ranked = sorted(
            [(float(score), 1) for score in pos] + [(float(score), 0) for score in neg],
            key=lambda item: item[0],
            reverse=higher_is_better,
        )
        n, n_pos = len(ranked), len(pos)
        if n == 0 or n_pos == 0:
            return 0.0
        if alpha <= 0:
            raise ValueError("alpha must be greater than zero")

        # BEDROC is rank based. Average the active contribution inside score-tie
        # groups so input order cannot make tied positives look artificially early.
        sum_exp = 0.0
        rank = 1
        index = 0
        while index < n:
            score = ranked[index][0]
            group = []
            while index < n and ranked[index][0] == score:
                group.append(ranked[index])
                index += 1
            active_fraction = sum(label for _, label in group) / len(group)
            sum_exp += active_fraction * sum(
                math.exp(-alpha * position / n)
                for position in range(rank, rank + len(group))
            )
            rank += len(group)

        denominator = ((1.0 - math.exp(-alpha))
                       / (n * (math.exp(alpha / n) - 1.0)))
        rie = sum_exp / (n_pos * denominator)
        ratio = n_pos / n
        rie_max = (1.0 - math.exp(-alpha * ratio)) / (ratio * (1.0 - math.exp(-alpha)))
        rie_min = (1.0 - math.exp(alpha * ratio)) / (ratio * (1.0 - math.exp(alpha)))
        if rie_max == rie_min:
            return 1.0
        value = (rie - rie_min) / (rie_max - rie_min)
        return round(clamp(value), 4)

    @classmethod
    def vs_benchmark(cls, results: List[Dict],
                     positive_ids: Optional[set] = None) -> Dict:
        if positive_ids is None:
            return {
                "status": "not_evaluated",
                "reason": "explicit positive_ids are required",
                "total": len(results),
                "positive": 0,
                "negative": 0,
                "roc_auc": None,
                "ef_1pct": None,
                "ef_5pct": None,
                "ef_10pct": None,
                "bedroc_alpha20": None,
            }
        def measured_score(row: Dict) -> Optional[float]:
            value = row.get("docking_score_kcal_mol")
            if value in (None, ""):
                return None
            try:
                score = float(value)
            except (TypeError, ValueError):
                return None
            return score if math.isfinite(score) else None

        pos = [score for row in results
               if row["id"] in positive_ids and (score := measured_score(row)) is not None]
        neg = [score for row in results
               if row["id"] not in positive_ids and (score := measured_score(row)) is not None]
        if not pos or not neg:
            return {
                "status": "not_evaluated",
                "reason": "positive and negative docking scores are both required",
                "total": len(results),
                "positive": len(pos),
                "negative": len(neg),
                "roc_auc": None,
                "ef_1pct": None,
                "ef_5pct": None,
                "ef_10pct": None,
                "bedroc_alpha20": None,
            }
        return {
            "status": "evaluated",
            "total": len(results), "positive": len(pos), "negative": len(neg),
            "roc_auc": cls.roc_auc(pos, neg, higher_is_better=True),
            "ef_1pct": cls.enrichment_factor(pos, neg, 0.01, higher_is_better=True),
            "ef_5pct": cls.enrichment_factor(pos, neg, 0.05, higher_is_better=True),
            "ef_10pct": cls.enrichment_factor(pos, neg, 0.10, higher_is_better=True),
            "bedroc_alpha20": cls.bedroc(pos, neg, 20.0, higher_is_better=True),
        }


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Layer 2 (新): DeepPurpose 结合亲和力预测
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

class Layer2DeepPurpose:
    """DeepPurpose BindingDB 结合亲和力预测 (替代分子对接)

    基于序列的深度学习模型, 输入SMILES+蛋白序列, 直接预测结合强度。
    使用BindingDB预训练CNN_CNN模型, 覆盖92万条蛋白-配体亲和力数据。

    用法:
        l2 = Layer2DeepPurpose(target_name="HIV-1_protease")
        result = l2.score("CC(=O)Oc1ccccc1C(=O)O")
    """

    TARGET_SEQUENCES = {
        "HIV-1_protease": "PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF",
    }

    def __init__(self, target_name: str = "HIV-1_protease",
                 target_sequence: Optional[str] = None,
                 model_name: str = "cnn_cnn_bindingdb"):
        self.target_name = target_name
        self.target_seq = target_sequence or self.TARGET_SEQUENCES.get(target_name)
        if self.target_seq is None:
            raise ValueError(f"未知靶点 '{target_name}', 请提供 target_sequence 参数")

        self.model_name = model_name
        self._model = None
        self._utils = None
        self._cache = {}

    @property
    def model(self):
        if self._model is None:
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                from DeepPurpose import DTI
                path_dir = str(Path(__file__).parent / "models" / "deeppurpose"
                               / "pretrained_models" / "model_cnn_cnn_bindingdb")
                self._model = DTI.model_pretrained(path_dir=path_dir)
        return self._model

    def _get_utils(self):
        if self._utils is None:
            from DeepPurpose import utils as dp_utils
            self._utils = dp_utils
        return self._utils

    @staticmethod
    def normalize(raw: float) -> float:
        """BindingDB原始值 → 0-1归一化: ~4(不结合)→0.0, ~9(强结合)→1.0"""
        return round(clamp((raw - 4.0) / 5.0), 4)

    def score(self, smiles: str, mol_id: str = "mol") -> Dict:
        if smiles in self._cache:
            raw = self._cache[smiles]
        else:
            try:
                utils = self._get_utils()
                X = utils.data_process_repurpose_virtual_screening(
                    X_repurpose=[smiles], target=self.target_seq,
                    drug_encoding='CNN', target_encoding='CNN', mode='repurposing')
                y = self.model.predict(X)
                raw = float(y[0])
                self._cache[smiles] = raw
            except Exception as e:
                return {"docking_score_kcal_mol": "", "docking_normalized": 0.0,
                        "docking_status": f"failed:{e}"}

        return {"docking_score_kcal_mol": round(raw, 3),
                "docking_normalized": self.normalize(raw),
                "docking_status": "success",
                "docking_method": "DeepPurpose-BindingDB"}


class Layer2Unavailable:
    """Explicit non-strict state for a requested BindingDB backend that did not load."""

    model_path = Path("bindingdb-unavailable")

    def __init__(self, error: Exception):
        self.error = error

    def score(self, smiles: str, target_text: str = "", mol_id: str = "mol") -> Dict:
        del smiles, target_text, mol_id
        return {
            "docking_score_kcal_mol": "",
            "docking_normalized": 0.0,
            "docking_status": f"failed:backend_unavailable:{type(self.error).__name__}",
            "docking_method": "BindingDB-L2-unavailable",
        }


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 统一接口
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

class MoleculeScorer:
    """
    一体化打分器 (四层漏斗)

    用法:
        scorer = MoleculeScorer()
        result = scorer.score_one("CC(=O)Oc1ccccc1C(=O)O")
        results = scorer.score_batch([("mol1", "CCO"), ("mol2", "c1ccccc1")])

        # 自定义靶点
        scorer = MoleculeScorer(deeppurpose_target="HIV-1_protease")
        scorer = MoleculeScorer(deeppurpose_target_seq="PQITLWQR...")
    """

    def __init__(self, deeppurpose_target: str = "HIV-1_protease",
                 deeppurpose_target_seq: Optional[str] = None,
                 use_unimol: bool = True,
                 l2_method: str = "bindingdb",
                 default_target_text: Optional[str] = None,
                 l2_model_path: Optional[str] = None,
                 strict_backends: bool = False,
                 unimol_device: str = "cpu",
                 asset_root: Optional[str | Path] = None):
        asset_paths = resolve_asset_paths(asset_root)
        if asset_paths is not None:
            os.environ["FOUR_LEVEL_ASSET_ROOT"] = str(asset_paths.root)
            if asset_paths.manifest.is_file():
                os.environ["FOUR_LEVEL_ASSET_MANIFEST"] = str(asset_paths.manifest)
        l2_model_path = l2_model_path or (str(asset_paths.l2_model) if asset_paths is not None else None)
        self.l1 = Layer1Scorer()
        self.l3 = Layer3Scorer(
            strict_backends=strict_backends,
            model_dir=asset_paths.admet_model_dir if asset_paths is not None else None,
        )
        self.l4 = Layer4Aggregator()
        self.unimol = None
        self.l2_method = l2_method
        self.default_target_text = default_target_text
        self.strict_backends = strict_backends

        if use_unimol:
            try:
                try:
                    from .scripts.unimol_scorer import UniMolScorer
                except ImportError:  # Direct script/PYTHONPATH=scoring compatibility.
                    from scripts.unimol_scorer import UniMolScorer
                self.unimol = UniMolScorer(
                    device=unimol_device,
                    model_dir=asset_paths.unimol_model_dir if asset_paths is not None else None,
                )
            except Exception as e:
                if strict_backends:
                    raise RuntimeError(f"L4 UniMol 加载失败: {e}") from e
                print(f"Uni-Mol加载失败: {e}")

        # ── Layer 2: 靶点感知结合亲和力 (取代失效的 DeepPurpose) ──
        self.l2 = None
        if l2_method in ("bindingdb", "bindingdb_seq"):
            try:
                if l2_method == "bindingdb":
                    try:
                        from .l2_bindingdb import Layer2BindingDB
                    except ImportError:  # Direct script/PYTHONPATH=scoring compatibility.
                        from l2_bindingdb import Layer2BindingDB
                    self.l2 = Layer2BindingDB(
                        prefer="mlp",
                        model_path=l2_model_path,
                        params_path=(str(asset_paths.l2_params) if asset_paths is not None else None),
                    )
                    print(f"L2: BindingDB 靶点感知模型(靶点哈希, 权重0.50, 模型={self.l2.model_kind})"
                          + (f" [覆盖={os.path.basename(l2_model_path)}]" if l2_model_path else ""))
                else:
                    try:
                        from .l2_bindingdb import Layer2BindingDBSeq
                    except ImportError:  # Direct script/PYTHONPATH=scoring compatibility.
                        from l2_bindingdb import Layer2BindingDBSeq
                    self.l2 = Layer2BindingDBSeq()
                    print(f"L2: BindingDB 序列嵌入靶点模型 (BLOSUM62+CNN, 权重0.50)")
            except Exception as e:
                if strict_backends:
                    raise RuntimeError(f"L2 BindingDB 加载失败: {e}") from e
                print(f"[WARN] BindingDB L2 加载失败: {e}; 标记 L2 后端不可用")
                self.l2 = Layer2Unavailable(e)
        if l2_method == "deeppurpose":
            self.l2_method = "deeppurpose"
            self.l2 = Layer2DeepPurpose(
                target_name=deeppurpose_target,
                target_sequence=deeppurpose_target_seq,
            )
            print(f"L2: DeepPurpose BindingDB (target={deeppurpose_target})")
        elif self.l2 is None:
            raise ValueError(f"不支持的 L2 后端: {l2_method}")

    def score_one(self, smiles: str, mol_id: str = "mol",
                  all_smiles: Optional[List[str]] = None,
                  target_text: Optional[str] = None,
                  _unimol_result: Optional[Dict] = None) -> Dict:
        """对单个SMILES打分, 返回所有层的结果

        target_text: 靶点文本 (BindingDB 命名空间), 供 L2 靶点感知使用。
                     省略时使用 self.default_target_text。
        """
        # Layer 1
        r1 = self.l1.score(smiles, all_smiles)
        l1_score = r1["layer1_score"]
        l1_status = "ok" if r1.get("valid", l1_score is not None) and l1_score is not None else "error"

        # Layer 2: 靶点感知结合亲和力预测
        tt = target_text or self.default_target_text
        if self.l2_method in ("bindingdb", "bindingdb_seq"):
            r2 = self.l2.score(smiles, target_text=tt or "", mol_id=mol_id)
        else:
            r2 = self.l2.score(smiles, mol_id)

        layer2_status = "ok" if r2.get("docking_status") == "success" else str(r2.get("docking_status", "error"))

        # Layer 3
        r3 = self.l3.score(smiles)
        layer3_status = (
            "error"
            if r3.get("admet_score") is None or r3.get("toxicity_flags") == "invalid"
            else "ok"
        )

        # Uni-Mol zero-shot
        unimol_score = 0.0
        unimol_pos = 0.0
        unimol_neg = 0.0
        if self.unimol:
            um = _unimol_result if _unimol_result is not None else self.unimol.score(smiles)
            if um.get("error") and self.strict_backends:
                raise RuntimeError(f"L4 UniMol 分子评分失败: {um['error']}")
            unimol_score = um.get("unimol_score", 0)
            unimol_pos = um.get("pos_similarity", 0)
            unimol_neg = um.get("neg_similarity", 0)
            l4_status = "error" if um.get("error") else "ok"
        else:
            l4_status = "disabled"

        if self.strict_backends:
            failures = {
                "L1": l1_status,
                "L2": layer2_status,
                "L3": layer3_status,
                "L4": l4_status if self.unimol else "disabled",
            }
            failed = {layer: status for layer, status in failures.items() if status not in {"ok", "disabled"}}
            if failed:
                details = "; ".join(f"{layer}={status}" for layer, status in failed.items())
                raise RuntimeError(f"strict backend scoring failed: {details}")

        # Layer 4: 最终分 = Layer1*0.20 + Layer2*0.50 + Layer3*0.20 + UniMol*0.10
        dock_norm = float(r2.get("docking_normalized", 0) or 0)
        admet = float(r3.get("admet_score", 0) or 0)
        final = self.l4.combine(l1_score, dock_norm, admet, unimol_score)
        gate, reason = self.l4.quality_gate({**r2, **r3})

        return {
            "id": mol_id, "smiles": smiles,
            "layer1_status": l1_status,
            "layer1_backend": "RDKit",
            "layer1_model_asset_id": "rdkit-runtime",
            "layer2_status": layer2_status,
            "layer2_backend": str(r2.get("docking_method", "unknown")),
            "layer2_model_asset_id": Path(getattr(self.l2, "model_path", "runtime-model")).name,
            "layer3_status": layer3_status,
            "layer3_backend": "ADMET-RF+RDKit",
            "layer3_model_asset_id": "tox21.pkl+bbbp.pkl+clintox.pkl+sider.pkl",
            "layer4_status": l4_status,
            "layer4_backend": "UniMolRepr-mol_pre_all_h_220816" if self.unimol else "disabled",
            "layer4_model_asset_id": "mol_pre_all_h_220816.pt+ref_embeddings.npz" if self.unimol else "none",
            # Layer 1
            "mw": r1.get("mw"), "logp": r1.get("logp"), "tpsa": r1.get("tpsa"),
            "hbd": r1.get("hbd"), "hba": r1.get("hba"),
            "qed": r1.get("qed"), "sa": r1.get("sa"),
            "lipinski_violations": r1.get("lipinski_violations"),
            "layer1_score": l1_score, "diversity": r1.get("diversity"),
            # Uni-Mol zero-shot
            "unimol_score": unimol_score,
            "unimol_pos_sim": unimol_pos,
            "unimol_neg_sim": unimol_neg,
            # Layer 2
            "docking_score_kcal_mol": r2.get("docking_score_kcal_mol"),
            "docking_normalized": r2.get("docking_normalized"),
            "docking_status": r2.get("docking_status"),
            "docking_method": r2.get("docking_method"),
            # Layer 3
            "admet_score": r3.get("admet_score"),
            "toxicity_count": r3.get("toxicity_count"),
            "toxicity_flags": r3.get("toxicity_flags"),
            "solubility_logS": r3.get("solubility_logS"),
            "bbb_prob": r3.get("bbb_prob"),
            "cyp3a4_risk": r3.get("cyp3a4_risk"),
            "cyp2d6_risk": r3.get("cyp2d6_risk"),
            "herg_risk": r3.get("herg_risk"),
            "hepatotoxicity_risk": r3.get("hepatotoxicity_risk"),
            "oral_bioavailability": r3.get("oral_bioavailability"),
            # Layer 4
            "final_score": final, "gate_status": gate, "gate_reason": reason,
        }

    def score_batch(self, molecules: List[Tuple[str, str]],
                     target_text: Optional[str] = None) -> List[Dict]:
        """批量打分: [(id, smiles), ...]
        target_text: 若提供, 所有分子按同一靶点打分 (L2 靶点感知)。"""
        all_smiles = [s for _, s in molecules]
        results = []
        unimol_rows = None
        if self.unimol:
            try:
                # Keep the production CLI's batch path numerically aligned with
                # the frozen target-independent UniMol cache.
                unimol_rows = self.unimol.score_many(all_smiles, batch_size=128)
            except Exception as exc:
                if self.strict_backends:
                    raise RuntimeError(f"L4 UniMol 批量评分失败: {exc}") from exc
        for i, (mid, smi) in enumerate(molecules, 1):
            print(f"[{i}/{len(molecules)}] {mid}", end=" ... ", flush=True)
            unimol_result = unimol_rows[i - 1] if unimol_rows is not None else None
            r = self.score_one(
                smi,
                mid,
                all_smiles,
                target_text=target_text,
                _unimol_result=unimol_result,
            )
            print(f"L1={r['layer1_score']:.3f}  L2={r['docking_normalized']}  "
                  f"L3={r['admet_score']:.3f}  final={r['final_score']:.4f}  {r['gate_status']}")
            results.append(r)
        results.sort(key=lambda x: x["final_score"], reverse=True)
        return results

    @staticmethod
    def save_csv(results: List[Dict], path: str):
        if not results:
            return
        fieldnames = list(results[0].keys())
        destination = Path(path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        with destination.open("w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(results)
        print(f"结果已写入: {path}")

    @staticmethod
    def print_ranking(results: List[Dict], top_n: int = 15):
        if not results:
            print("没有可排名的分子")
            return
        score_key = "final_score_dock" if any("final_score_dock" in row for row in results) else "final_score"
        score_label = "级联分" if score_key == "final_score_dock" else "最终分"
        print(f"\n{'排名':>4}  {'ID':<16} {score_label:>7} {'L1':>6} {'L2':>6} {'L3':>6} {'L4':>6} {'门槛':>6}")
        print(f"{'-' * 80}")
        for i, r in enumerate(results[:top_n], 1):
            l4 = r.get("unimol_score", 0)
            print(f"{i:>4}  {r['id']:<16} {r.get(score_key, r['final_score']):>7} "
                  f"{r['layer1_score']:>6} {r['docking_normalized']:>6} "
                  f"{r['admet_score']:>6} {l4:>6} {r['gate_status']:>6}")


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# CLI入口
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def read_molecules_csv(path: str | Path) -> List[Tuple[str, str]]:
    with Path(path).open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = set(reader.fieldnames or [])
        if not {"id", "smiles"}.issubset(fields):
            raise ValueError("输入 CSV 必须包含 id,smiles 表头")
        molecules = []
        seen_ids = set()
        for line_number, row in enumerate(reader, start=2):
            mol_id = str(row.get("id") or "").strip()
            smiles = str(row.get("smiles") or "").strip()
            if not mol_id or not smiles:
                raise ValueError(f"输入 CSV 第 {line_number} 行的 id/smiles 不能为空")
            if mol_id in seen_ids:
                raise ValueError(f"输入 CSV 包含重复 id: {mol_id}")
            seen_ids.add(mol_id)
            molecules.append((mol_id, smiles))
    return molecules


def _parse_docking_triplet(value: str, *, name: str, positive: bool = False) -> Tuple[float, float, float]:
    try:
        parsed = tuple(float(item.strip()) for item in value.split(","))
    except (AttributeError, TypeError, ValueError) as exc:
        raise ValueError(f"{name} 必须是三个逗号分隔的数值") from exc
    if len(parsed) != 3 or not all(math.isfinite(item) for item in parsed):
        raise ValueError(f"{name} 必须是三个有限数值")
    if positive and not all(item > 0 for item in parsed):
        raise ValueError(f"{name} 的三个值都必须大于 0")
    return parsed


def _parse_unit_interval(value: str) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise argparse.ArgumentTypeError("数值必须位于 0 到 1 之间") from exc
    if not math.isfinite(parsed) or not 0.0 <= parsed <= 1.0:
        raise argparse.ArgumentTypeError("数值必须是 0 到 1 之间的有限值")
    return parsed


def _parse_positive_int(value: str) -> int:
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise argparse.ArgumentTypeError("数值必须是整数") from exc
    if parsed <= 0:
        raise argparse.ArgumentTypeError("数值必须大于 0")
    return parsed


def validate_docking_request(
    *,
    receptor: Optional[str],
    box_center: Optional[str],
    box_size: Optional[str],
    dock_rerank: bool,
    manual_receptor: bool,
) -> Tuple[Optional[Tuple[float, float, float]], Optional[Tuple[float, float, float]]]:
    if dock_rerank and not receptor:
        raise ValueError("--dock-rerank 需要 --receptor 或有效的注册受体")
    if (box_center or box_size) and not receptor:
        raise ValueError("对接盒参数需要有效受体")
    if not receptor:
        return None, None
    if not Path(receptor).expanduser().is_file():
        raise ValueError(f"受体文件不存在: {receptor}")
    if not box_center or not box_size:
        source = "手工受体" if manual_receptor else "注册受体"
        raise ValueError(f"{source}需要同时提供 box-center 和 box-size 对接盒参数")
    center = _parse_docking_triplet(box_center, name="box-center")
    size = _parse_docking_triplet(box_size, name="box-size", positive=True)
    return center, size


def main():
    parser = argparse.ArgumentParser(description="AI药物分子打分系统 (四层漏斗 + 可选对接精排[LE校正])")
    parser.add_argument("--input", "-i", required=True, help="候选分子CSV (需含id,smiles列)")
    parser.add_argument("--output", "-o", default="scores.csv", help="输出CSV路径")
    parser.add_argument("--target", "-t", default="HIV-1_protease",
                        help="靶点: ChEMBL ID / 靶点名 / 蛋白名 (bindingdb 模式自动解析为靶点文本)")
    parser.add_argument("--target-seq", help=argparse.SUPPRESS)
    parser.add_argument("--l2", default="bindingdb",
                        choices=["bindingdb"],
                        help="L2 结合亲和力模型（正式 CLI 仅发布已验证的 bindingdb 后端）")
    parser.add_argument("--target-text", default=None,
                        help="显式指定 BindingDB 命名空间的靶点文本; 省略时尝试用 -t 经同义词表解析")
    # ---- 选择树: 按可用证据(受体?)自动路由 级联/文库 分支 ----
    parser.add_argument("--mode", default="auto", choices=["auto", "cascade", "library"],
                        help="管线选择树: auto=按受体注册表自动路由(默认) / cascade=强制级联对接精排 / library=强制文库分支(L2+ADMET+UniMol)")
    # ---- 可选 M2 对接腿 ----
    parser.add_argument("--receptor", default=None,
                        help="受体 PDB/PDBQT 路径; 提供后启用 smina 对接腿 (并输出亲和力列)")
    parser.add_argument("--box-center", default=None,
                        help="对接盒中心 x,y,z (逗号分隔), 与 --receptor 搭配")
    parser.add_argument("--box-size", default=None,
                        help="对接盒尺寸 x,y,z (逗号分隔), 与 --receptor 搭配")
    parser.add_argument("--flex", default=None,
                        help="柔性残基 (chain:resid 逗号分隔); 提供后启用 smina --flexres 柔性对接")
    parser.add_argument("--mcce", action="store_true",
                        help="对接前用 Open Babel -p 7.4 做生理 pH 质子化 (兼容旧 --mcce 名称)")
    # ---- 对接精排 (② 验证: smina + 配体效率LE校正) ----
    parser.add_argument("--dock-rerank", action="store_true",
                        help="启用对接精排: 对所有分子跑 smina 刚性对接, 并用配体效率(LE)校正重排 "
                             "(需配合 --receptor; 默认 LE 模式, 校正 Vina 尺寸/亲脂偏差)")
    parser.add_argument("--dock-mode", default=None, choices=["le", "raw"],
                        help="对接重排信号: le=配体效率(affinity/重原子数, 校正尺寸偏差, 默认) / raw=原始 affinity")
    parser.add_argument("--dock-fusion", type=_parse_unit_interval, default=None,
                        help="修正级联融合权重 W_LE(0-1): 化学基线分 + W_LE·σ(base)·z(−LE); "
                             "自动级联默认0.30 (0=仅作独立列与重排CSV, 不混分)")
    parser.add_argument("--cascade-top-n", type=_parse_positive_int, default=None,
                        help="级联分支仅对 L2 粗筛前 N 个候选运行结构对接 (默认300)")
    parser.add_argument("--strict-backends", action="store_true",
                        help="任一四级后端加载或分子评分失败时立即退出，不允许静默降级")
    parser.add_argument("--asset-root", default=None,
                        help="外部离线资产根目录 (含 ASSET_MANIFEST.json、scoring/models 和 scoring/receptors)")
    args = parser.parse_args()

    if args.asset_root:
        configured_assets = resolve_asset_paths(args.asset_root)
        os.environ["FOUR_LEVEL_ASSET_ROOT"] = str(configured_assets.root)
        os.environ["FOUR_LEVEL_ASSET_MANIFEST"] = str(configured_assets.manifest)

    if args.mode == "library":
        conflicts = [
            option
            for option, active in (
                ("--receptor", args.receptor is not None),
                ("--box-center", args.box_center is not None),
                ("--box-size", args.box_size is not None),
                ("--flex", args.flex is not None),
                ("--mcce", args.mcce),
                ("--dock-rerank", args.dock_rerank),
                ("--dock-fusion", args.dock_fusion is not None),
                ("--dock-mode", args.dock_mode is not None),
                ("--cascade-top-n", args.cascade_top_n is not None),
            )
            if active
        ]
        if conflicts:
            parser.error(f"--mode library 不允许对接参数: {', '.join(conflicts)}")
    manual_receptor = args.receptor is not None

    # 读取输入
    try:
        rows = read_molecules_csv(args.input)
    except ValueError as exc:
        parser.error(str(exc))

    print(f"读入 {len(rows)} 个分子")
    if not rows:
        parser.error("输入文件不包含可评分分子（需要至少一行 id,smiles）")

    # ---- 选择树: 按可用证据(受体?)路由到 级联 / 文库 分支 ----
    cid = args.target if args.target.startswith("CHEMBL") else None
    l2_model_path = os.environ.get("L2_MODEL_PATH")
    decision = None
    try:
        try:
            from .pipeline_router import route
        except ImportError:  # Direct script/PYTHONPATH=scoring compatibility.
            from pipeline_router import route
        route_mode = "auto" if manual_receptor and args.mode == "cascade" else args.mode
        decision = route(args.target, chembl_id=cid, mode=route_mode, l2_model_path=l2_model_path)
        print(f"[选择树] branch={decision['branch']}  reason={decision['rationale']}")
    except Exception as e:
        if args.mode == "cascade":
            raise
        print(f"[WARN] 选择树路由失败: {e}; 回退文库分支")
    if decision is None:
        decision = {"branch": "library", "target_text": args.target,
                    "receptor": None, "box_center": None, "box_size": None,
                    "l2_model_path": l2_model_path,
                    "rationale": "回退文库分支"}
    if manual_receptor and args.mode in {"auto", "cascade"}:
        decision["branch"] = "cascade"
        decision["rationale"] = "显式受体与对接盒 -> 强制级联分支"

    if decision["branch"] == "library":
        conflicts = [
            option
            for option, active in (
                ("--box-center", args.box_center is not None),
                ("--box-size", args.box_size is not None),
                ("--flex", args.flex is not None),
                ("--mcce", args.mcce),
                ("--dock-rerank", args.dock_rerank),
                ("--dock-fusion", args.dock_fusion is not None),
                ("--dock-mode", args.dock_mode is not None),
                ("--cascade-top-n", args.cascade_top_n is not None),
            )
            if active
        ]
        if conflicts:
            parser.error(f"有效文库分支不允许对接参数: {', '.join(conflicts)}")

    args.dock_mode = args.dock_mode or "le"
    args.cascade_top_n = args.cascade_top_n or 300

    # 靶点文本 (L2 靶点感知)
    default_tt = args.target_text or decision["target_text"]
    if default_tt and len(default_tt) > 70:
        print(f"靶点解析: '{args.target}' -> {default_tt[:70]}...")

    # 自动级联: 命中受体注册表 -> 注入对接腿参数 + 默认融合权重
    if decision["branch"] == "cascade" and (args.receptor or decision.get("receptor")):
        if not args.receptor:
            args.receptor = decision["receptor"]
        if not args.box_center and decision.get("box_center"):
            args.box_center = ",".join(map(str, decision["box_center"]))
        if not args.box_size and decision.get("box_size"):
            args.box_size = ",".join(map(str, decision["box_size"]))
        args.dock_rerank = True
        if args.dock_fusion is None:
            args.dock_fusion = 0.30
        print(f"[选择树] 自动启用级联对接精排 (受体={os.path.basename(args.receptor)}, 融合权重={args.dock_fusion})")

    try:
        dock_center, dock_size = validate_docking_request(
            receptor=args.receptor,
            box_center=args.box_center,
            box_size=args.box_size,
            dock_rerank=args.dock_rerank,
            manual_receptor=manual_receptor,
        )
    except ValueError as exc:
        parser.error(str(exc))

    scorer = MoleculeScorer(
        deeppurpose_target=args.target,
        deeppurpose_target_seq=args.target_seq,
        l2_method=args.l2,
        default_target_text=default_tt,
        l2_model_path=decision.get("l2_model_path"),
        strict_backends=args.strict_backends,
        asset_root=args.asset_root,
    )

    # 先完成四级粗筛，级联只对 L2 top-N 运行结构对接。
    results = scorer.score_batch(rows)
    l2_scores = {r["id"]: float(r.get("docking_normalized", 0.0) or 0.0) for r in results}

    # ---- 可选 对接腿 (smina + obabel, 配体效率LE校正) ----
    dock_map = {}   # id -> {affinity, heavy_atoms, ligand_efficiency, dock_rerank_rank, status}
    if args.receptor:
        try:
            try:
                from .dock_rerank import (
                    DockingReranker,
                    find_binary,
                    cascade_corrected_fusion,
                    protonated_receptor,
                )
            except ImportError:  # Direct script/PYTHONPATH=scoring compatibility.
                from dock_rerank import (
                    DockingReranker,
                    find_binary,
                    cascade_corrected_fusion,
                    protonated_receptor,
                )
            ob = find_binary("obabel") if args.mcce else None
            if args.mcce and not ob:
                raise RuntimeError("Open Babel binary is unavailable for --mcce")
            receptor_scope = (
                protonated_receptor(args.receptor, ob)
                if args.mcce
                else nullcontext(args.receptor)
            )
            with receptor_scope as rec:
                rr = DockingReranker(receptor=rec, center=dock_center, size=dock_size, flexres=args.flex or "")
                dock_rows = (select_cascade_candidates(rows, l2_scores, args.cascade_top_n)
                             if args.dock_rerank else rows)
                recs = rr.dock_all(dock_rows, mode=args.dock_mode)
                dock_map = {r["id"]: r for r in recs}
                n_ok = sum(1 for r in recs if r["affinity"] is not None)
                if n_ok == 0 and (args.mode == "cascade" or args.strict_backends):
                    raise RuntimeError("级联对接零成功结果，拒绝生成伪级联输出")
                print(f"对接腿完成: {n_ok}/{len(dock_rows)} 成功; 总候选={len(rows)} (mode={args.dock_mode})")
        except Exception as e:
            if args.mode == "cascade" or args.strict_backends:
                raise
            import traceback; traceback.print_exc()
            print(f"[WARN] 对接腿失败, 仅输出 L1/L2/L3: {e}")

    # 合并对接结果 (affinity / 重原子数 / 配体效率 / 重排名次)
    for r in results:
        d = dock_map.get(r["id"])
        if d:
            r["docking_affinity_kcal_mol"] = d["affinity"]
            r["heavy_atoms"] = d["heavy_atoms"]
            r["ligand_efficiency"] = d["ligand_efficiency"]
            r["dock_rerank_rank"] = d.get("dock_rerank_rank")
            r["structure_docking_status"] = d["status"]
        elif args.receptor:
            r["docking_affinity_kcal_mol"] = None
            r["heavy_atoms"] = None
            r["ligand_efficiency"] = None
            r["dock_rerank_rank"] = None
            r["structure_docking_status"] = "not_selected"

    # 对接精排: 归一化 LE 列 + 可选融合 + 重排 CSV
    if args.dock_rerank and dock_map:
        if args.receptor is None:
            print("[WARN] --dock-rerank 需要 --receptor, 已跳过精排")
        else:
            les = [r["ligand_efficiency"] for r in results if r.get("ligand_efficiency") is not None]
            if les:
                lo, hi = min(les), max(les)
                for r in results:
                    le = r.get("ligand_efficiency")
                    r["dock_le_norm"] = round((hi - le) / (hi - lo), 4) if (le is not None and hi > lo) else 0.0
            if args.dock_fusion and args.dock_fusion > 0:
                w = args.dock_fusion
                # 历史展示采用的融合公式；0.935 是人工平衡集结果，不是通用性能声明。
                base_arr = np.array([r["final_score"] for r in results], dtype=float)
                le_arr = [r.get("ligand_efficiency") for r in results]
                fused = cascade_corrected_fusion(base_arr, le_arr, w)
                for r, f in zip(results, fused):
                    r["final_score_dock"] = round(float(f), 4)
                results.sort(key=lambda x: x.get("final_score_dock", x["final_score"]), reverse=True)
                n_docked = sum(1 for r in results if r.get("ligand_efficiency") is not None)
                print(f"[对接精排] 修正融合: base + {w:.2f}·σ(base)·z(−LE); 成功对接 {n_docked} 个")
            # 重排 CSV (纯物理精排视图, 按 dock_rerank_rank)
            rerank_sorted = sorted([r for r in results if r.get("dock_rerank_rank") is not None],
                                   key=lambda x: x["dock_rerank_rank"])
            rerank_out = os.path.splitext(args.output)[0] + "_reranked.csv"
            scorer.save_csv(rerank_sorted, rerank_out)

    # 输出主 CSV；启用融合时保留 base 分并按 final_score_dock 排序。
    scorer.save_csv(results, args.output)
    scorer.print_ranking(results)

    # 对接精排 Top 榜单
    if args.dock_rerank and dock_map and args.receptor:
        mode_label = "配体效率LE(affinity/重原子)" if args.dock_mode == "le" else "原始affinity"
        final_key = "final_score_dock" if any("final_score_dock" in row for row in results) else "final_score"
        final_label = "fused" if final_key == "final_score_dock" else "final"
        print(f"\n[对接精排 Top] 按 {mode_label} 重排 (rank 1=最佳):")
        print(f"{'排名':>4}  {'ID':<16} {'affinity':>9} {'HAC':>4} {'LE':>7} {'L2-norm':>7} {final_label:>7}")
        print("-" * 70)
        top = sorted([r for r in results if r.get("dock_rerank_rank") is not None],
                     key=lambda x: x["dock_rerank_rank"])[:15]
        for r in top:
            print(f"{r['dock_rerank_rank']:>4}  {r['id']:<16} "
                  f"{str(r.get('docking_affinity_kcal_mol')):>9} "
                  f"{str(r.get('heavy_atoms')):>4} {str(r.get('ligand_efficiency')):>7} "
                  f"{str(r.get('docking_normalized')):>7} {r.get(final_key, r.get('final_score')):>7}")

    # 虚拟筛选基准
    bench = Layer4Aggregator.vs_benchmark(results)
    if bench["status"] == "evaluated":
        print(f"\n虚拟筛选基准: AUC={bench['roc_auc']}  EF5%={bench['ef_5pct']}  "
              f"BEDROC={bench['bedroc_alpha20']}")
    else:
        print(f"\n虚拟筛选基准: not_evaluated ({bench['reason']})")


if __name__ == "__main__":
    main()

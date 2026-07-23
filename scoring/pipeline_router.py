# -*- coding: utf-8 -*-
from __future__ import annotations
"""
pipeline_router —— 药物分子打分"选择树" (按可用证据类型路由, 非按模型版本挑)

为什么这样设计 (见 investigate 结论):
  - 在 L2 的多个版本间挑 (哈希 0.689 / CE重训 0.571 / BPR 0.506) 是【无价值】的:
      逐靶完美路由(oracle)天花板仅 0.734 (+0.045), 且无法用廉价特征(n_pos等)恢复;
      CE/BPR 重训在聚合域反而退化 —— 单模型重训修不好这个"bug"。
  - 真正有价值的分支在【管线层级】: 该靶点"有没有可用 3D 受体"?
      * 有受体 -> 级联分支 (L2 粗筛 + 对接精排 LE + 融合)。历史 0.935
        来自人工平衡/特定展示集，不代表真实大池或跨靶点泛化。
      * 无受体 -> 文库分支 (L2 + ADMET + UniMol), 通用中位 0.689
    这是基于"可用证据"的规则分支, 干净、可解释、已被数据证明。

route() 返回决策:
  {
    "branch": "cascade" | "library",
    "target_text": str,            # 供 L2 靶点感知
    "method": str,                 # 靶点解析方式
    "receptor": str|None,          # pdbqt 路径 (cascade 时非空)
    "box_center": [x,y,z]|None,
    "box_size":   [x,y,z]|None,
    "l2_model_path": str|None,     # 若设了 L2_MODEL_PATH 则用其覆盖哈希默认
    "rationale": str               # 人类可读的路由理由
  }
"""
import os, sys, json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REGISTRY_PATH = Path(os.environ.get("FOUR_LEVEL_RECEPTOR_REGISTRY", HERE / "receptor_registry.json")).resolve()


def _resolved_entry(entry: dict) -> dict:
    """Resolve registry-relative asset paths without embedding a host path."""
    resolved = dict(entry)
    receptor = Path(str(resolved.get("pdbqt", "")))
    if not receptor.is_absolute():
        external_root = os.environ.get("FOUR_LEVEL_ASSET_ROOT", "").strip()
        external = (Path(external_root).expanduser() / "scoring" / receptor).resolve() if external_root else None
        receptor = external if external is not None and external.is_file() else (REGISTRY_PATH.parent / receptor).resolve()
    resolved["pdbqt"] = str(receptor)
    return resolved


def load_registry() -> dict:
    if REGISTRY_PATH.exists():
        try:
            return json.load(open(REGISTRY_PATH))
        except Exception:
            pass
    return {"entries": {}}


def _normalize(name: str) -> str:
    return (name or "").strip().lower().replace("_", " ")


def lookup_receptor(target: str, chembl_id: str = None) -> dict | None:
    reg = load_registry().get("entries", {})
    # 优先 chembl_id 精确匹配
    if chembl_id and chembl_id in reg:
        return _resolved_entry(reg[chembl_id])
    # 其次规范化名匹配
    nt = _normalize(target)
    for k, v in reg.items():
        if _normalize(k) == nt:
            return _resolved_entry(v)
        if _normalize(v.get("target_name", "")) == nt:
            return _resolved_entry(v)
    return None


def validate_receptor(entry: dict | None) -> tuple[bool, str]:
    if not entry:
        return False, "未注册受体"
    receptor = Path(entry.get("pdbqt", ""))
    if not receptor.is_file():
        return False, f"受体文件不存在: {receptor}"
    center = entry.get("box_center")
    size = entry.get("box_size")
    if not isinstance(center, list) or len(center) != 3:
        return False, "box_center 必须是三个数值"
    if not isinstance(size, list) or len(size) != 3:
        return False, "box_size 必须是三个数值"
    try:
        if any(float(value) <= 0 for value in size):
            return False, "box_size 必须为正数"
        [float(value) for value in center]
    except (TypeError, ValueError):
        return False, "对接盒参数不是数值"
    return True, "ok"


def resolve_target_text(target: str, chembl_id: str = None):
    """返回 (target_text, method); 失败则优雅回退。"""
    try:
        sys.path.insert(0, str(HERE))
        try:
            from .target_resolver import TargetResolver
        except ImportError:  # Direct script/PYTHONPATH=scoring compatibility.
            from target_resolver import TargetResolver
        tt, method = TargetResolver().resolve(target, chembl_id=chembl_id)
        if tt:
            return tt, method
    except Exception as e:
        print(f"[router] 靶点解析失败: {e}")
    return target, "name_fallback"


def route(target: str, chembl_id: str = None, mode: str = "auto",
          l2_model_path: str = None) -> dict:
    """
    mode:
      "auto"    -> 按受体可用性自动选 (默认)
      "cascade" -> 强制级联 (无受体则报错提示)
      "library" -> 强制文库分支
    """
    target_text, method = resolve_target_text(target, chembl_id)
    l2_model_path = l2_model_path or os.environ.get("L2_MODEL_PATH") or None

    if mode == "library":
        return {
            "branch": "library", "target_text": target_text, "method": method,
            "receptor": None, "box_center": None, "box_size": None,
            "l2_model_path": l2_model_path,
            "rationale": "强制文库分支 (L2 + ADMET + UniMol)",
        }

    rec = lookup_receptor(target, chembl_id)
    receptor_ok, receptor_reason = validate_receptor(rec)
    if mode == "cascade" and not receptor_ok:
        raise RuntimeError(f"强制级联不可用: {receptor_reason}")

    if receptor_ok:
        return {
            "branch": "cascade", "target_text": target_text, "method": method,
            "receptor": rec.get("pdbqt"),
            "box_center": rec.get("box_center"),
            "box_size": rec.get("box_size"),
            "l2_model_path": l2_model_path,
            "rationale": (f"命中有效受体资产 -> 级联分支 "
                          f"(L2 粗筛 + 对接精排 LE + 修正融合); 受体={Path(rec.get('pdbqt','')).name}"),
        }

    if rec is not None:
        return {
            "branch": "library", "target_text": target_text, "method": method,
            "receptor": None, "box_center": None, "box_size": None,
            "l2_model_path": l2_model_path,
            "rationale": f"受体资产无效 -> 文库分支: {receptor_reason}",
        }

    # 默认: 无受体 -> 文库分支
    return {
        "branch": "library", "target_text": target_text, "method": method,
        "receptor": None, "box_center": None, "box_size": None,
        "l2_model_path": l2_model_path,
        "rationale": "无可用受体 -> 文库分支 (L2 + ADMET + UniMol)",
    }


if __name__ == "__main__":
    import pprint
    tests = [("CHEMBL2051", "CHEMBL2051"), ("neuraminidase", None),
             ("EGFR", None), ("HIV-1_protease", None)]
    for t, c in tests:
        print(f"\n=== route(target={t!r}, chembl={c}) ===")
        pprint.pprint(route(t, c, mode="auto"))

# -*- coding: utf-8 -*-
"""
靶点解析器 (同义词表): 把任意查询 (靶点名 / ChEMBL ID / 基因名) 解析为
MoleculeScorer L2 训练时见过的 target_text 字符串。

解决的问题:
  上一轮 L2 只训在 BindingDB, 164 个面板靶点里 48 个无法子串匹配到
  BindingDB target_text -> 推理时只能拿面板 name 当弱信号文本。
  实际上这 48 个靶点都来自 ChEMBL, 而 ChEMBL examples 自带 target_text
  与 target_chembl_id。用 target_chembl_id 精确映射即可补上这 48 个。

解析优先级:
  1. target_chembl_id 精确命中 ChEMBL target_text  (最准, 覆盖全部 ChEMBL 源靶点)
  2. 蛋白名子串命中 BindingDB target_text            (覆盖纯 BindingDB 靶点)
  3. 蛋白名子串命中 ChEMBL target_name -> 取其 target_text
  4. 基因符号 (大写, 含在名中) 兜底
  5. 都不中 -> 返回 None (调用方降级为 name 弱信号)

训练侧 (train_l2_combined.py) 用同一套 target_text 字符串, 保证推理/训练一致。
"""
from __future__ import annotations

import csv
import json
import os
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

REPO = Path(os.environ.get("FOUR_LEVEL_ROOT", Path(__file__).resolve().parents[1])).resolve()
BINDINGDB_EX = Path(os.environ.get(
    "FOUR_LEVEL_BINDINGDB_EXAMPLES",
    REPO / "data" / "bindingdb_202606_target_match_examples.csv",
)).resolve()
CHEMBL_EX = Path(os.environ.get(
    "FOUR_LEVEL_CHEMBL_EXAMPLES",
    REPO / "data" / "chembl_37_target_match_examples.csv",
)).resolve()
PANEL = Path(os.environ.get(
    "FOUR_LEVEL_TARGET_PANEL",
    REPO / "data" / "target_panel.json",
)).resolve()

UNIPROT = re.compile(r"\b[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}\b")
GENE = re.compile(r"\b[A-Z][A-Z0-9]{2,5}\b")  # 粗糙基因符号


class TargetResolver:
    def __init__(self):
        self.chembl_id_to_text: Dict[str, str] = {}
        self.chembl_name_to_text: Dict[str, str] = {}
        self.bindingdb_texts: List[str] = []
        self._bd_counts: Dict[str, int] = {}
        self._chembl_text_counts: Dict[str, int] = {}
        self._load()

    def _load(self):
        # BindingDB
        if BINDINGDB_EX.exists():
            with open(BINDINGDB_EX) as f:
                for row in csv.DictReader(f):
                    t = row.get("target_text", "").strip()
                    if t:
                        self.bindingdb_texts.append(t)
                        self._bd_counts[t] = self._bd_counts.get(t, 0) + 1
        # ChEMBL
        if CHEMBL_EX.exists():
            with open(CHEMBL_EX) as f:
                for row in csv.DictReader(f):
                    cid = row.get("target_chembl_id", "").strip()
                    t = row.get("target_text", "").strip()
                    name = row.get("target_name", "").strip()
                    if cid and t:
                        # 同一 chembl id 可能有多行 target_text, 取出现最多的
                        self._chembl_text_counts[t] = self._chembl_text_counts.get(t, 0) + 1
                        if cid not in self.chembl_id_to_text or \
                           self._chembl_text_counts[t] > self._chembl_text_counts.get(self.chembl_id_to_text.get(cid, ""), 0):
                            self.chembl_id_to_text[cid] = t
                    if name and t:
                        if name not in self.chembl_name_to_text or \
                           self._chembl_text_counts[t] > self._chembl_text_counts.get(self.chembl_name_to_text.get(name, ""), 0):
                            self.chembl_name_to_text[name] = t

    def resolve(self, name: str, chembl_id: Optional[str] = None) -> Tuple[Optional[str], str]:
        """返回 (target_text, method)。method 用于记录解析方式。"""
        # 1. chembl id 精确
        if chembl_id and chembl_id in self.chembl_id_to_text:
            return self.chembl_id_to_text[chembl_id], "chembl_id_exact"
        # 2. 蛋白名子串命中 BindingDB
        ln = (name or "").lower()
        cands = [(t, self._bd_counts[t]) for t in self.bindingdb_texts if ln and ln in t.lower()]
        if cands:
            cands.sort(key=lambda x: -x[1])
            return cands[0][0], "bindingdb_substr"
        # 3. 蛋白名子串命中 ChEMBL target_name
        cands = [(self.chembl_name_to_text[n], 1) for n in self.chembl_name_to_text
                 if ln and ln in n.lower()]
        if cands:
            return cands[0][0], "chembl_name_substr"
        # 4. 基因符号
        genes = GENE.findall(name or "")
        for g in genes:
            for n, t in self.chembl_name_to_text.items():
                if g.lower() in n.lower():
                    return t, "gene_symbol"
        return None, "none"


def build_panel_synonyms(out_json: Optional[Path] = None) -> Dict:
    """为 164 个面板靶点构建同义词解析表, 写出 JSON。

    返回 dict: { chembl_id: {"name":.., "resolved_text":.., "method":.., "mapped": bool} }
    """
    out_json = out_json or (REPO / "scientific_validation/multitarget_benchmark/target_synonyms.json")
    panel = json.load(open(PANEL))
    res = TargetResolver()
    out = {}
    n_mapped = 0
    for cid, v in panel.items():
        name = v.get("name", "")
        text, method = res.resolve(name, chembl_id=cid)
        mapped = method != "none"
        n_mapped += int(mapped)
        out[cid] = {"name": name, "resolved_text": text, "method": method, "mapped": mapped}
    json.dump(out, open(out_json, "w"), indent=1, ensure_ascii=False)
    print(f"面板同义词表: {n_mapped}/{len(panel)} 已映射 -> {out_json}")
    # 统计方法分布
    from collections import Counter
    cnt = Counter(d["method"] for d in out.values())
    print("解析方法分布:", dict(cnt))
    return out


if __name__ == "__main__":
    build_panel_synonyms()

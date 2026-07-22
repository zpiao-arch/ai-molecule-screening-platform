# -*- coding: utf-8 -*-
"""深层调查(单靶点模式): 每个靶点独立进程, 互不牵连; 缓存训练集。
用法: python investigate_auc2.py <chembl_id>   (无参数=列出样本cid)
输出: 追加一行到 `AUC2_ROWS` 环境变量指定的 TSV（默认当前目录）。
"""
import csv, json, sys, os, random, pickle
import numpy as np
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
sys.path.insert(0, str(Path("<validated-workspace>/评分_work_package/评分")))
from target_resolver import TargetResolver
from l2_bindingdb import Layer2BindingDB

REPO = Path("<validated-workspace>")
BD = REPO / "data_lake/bindingdb/aligned_model_input/bindingdb_202606_target_match_examples.csv"
CH = REPO / "data_lake/chembl/latest/aligned_model_input/chembl_37_target_match_examples.csv"
PANEL = json.load(open(REPO / "scientific_validation/multitarget_benchmark/target_panel.json"))
RES = {r["target"]: r for r in json.load(open(REPO / "scientific_validation/multitarget_benchmark/hts_closed_loop_results.json"))["per_target"]}
random.seed(0)

def canon(s):
    try: return Chem.MolToSmiles(Chem.MolFromSmiles(s))
    except Exception: return None
def mfp(s):
    try: return AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(s), 2, nBits=512)
    except Exception: return None

OUTPUT_DIR = Path(os.environ.get("AUC2_OUTPUT_DIR", "."))
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
CACHE_LIB = str(OUTPUT_DIR / "_auc2_lib.pkl")
CACHE_TR = str(OUTPUT_DIR / "_auc2_train.pkl")
AUC2_ROWS = Path(os.environ.get("AUC2_ROWS", str(OUTPUT_DIR / "auc2_rows.tsv")))

def load_lib():
    if os.path.exists(CACHE_LIB):
        return pickle.load(open(CACHE_LIB, "rb"))
    lib = set()
    with open(REPO / "scientific_validation/multitarget_benchmark/candidates_10k.csv") as f:
        for r in csv.DictReader(f):
            c = canon(r["smiles"])
            if c: lib.add(c)
    lib = list(lib)
    pickle.dump(lib, open(CACHE_LIB, "wb"))
    return lib

def load_train():
    if os.path.exists(CACHE_TR):
        return pickle.load(open(CACHE_TR, "rb"))
    train_act = {}
    for path in (BD, CH):
        if not path.exists(): continue
        with open(path) as f:
            for row in csv.DictReader(f):
                t = (row.get("target_text", "") or "").strip()
                lab = str(row.get("label", ""))
                smi = row.get("canonical_smiles") or row.get("Ligand SMILES") or ""
                if t and lab in ("1", "1.0") and smi:
                    c = canon(smi)
                    if c: train_act.setdefault(t, set()).add(c)
    pickle.dump(train_act, open(CACHE_TR, "wb"))
    return train_act

def analyze(cid):
    lib = load_lib()
    train_act = load_train()
    res = TargetResolver()
    l2 = Layer2BindingDB()
    r = RES[cid]
    text, method = res.resolve(PANEL[cid].get("name", ""), chembl_id=cid)
    libset = set(lib)
    libact = [canon(s) for s in PANEL[cid].get("pos", []) if canon(s) in libset]
    if not libact:
        libact = random.sample(lib, 30)
    tract = list(train_act.get(text, set()))
    inact = random.sample([s for s in lib if s not in set(libact)], 800)
    def l2s(s):
        try: return float(l2.score(s, text).get("docking_normalized") or 0.0)
        except Exception: return 0.0
    sa = np.array([l2s(s) for s in libact])
    si = np.array([l2s(s) for s in inact])
    sep = float(sa.mean() - si.mean())
    from sklearn.metrics import roc_auc_score
    y = np.concatenate([np.ones(len(sa)), np.zeros(len(si))])
    try: l2auc = float(roc_auc_score(y, np.concatenate([sa, si])))
    except Exception: l2auc = float("nan")
    medmaxT = float("nan")
    if tract:
        tfp = [mfp(s) for s in tract if mfp(s) is not None]
        if tfp:
            medmax = []
            for s in libact:
                fp = mfp(s)
                if fp is None: continue
                try:
                    ts = [AllChem.DataStructs.TanimotoSimilarity(fp, t) for t in tfp]
                    medmax.append(max(ts))
                except Exception: pass
            medmaxT = float(np.median(medmax)) if medmax else float("nan")
    grp = "REVERSED" if r["auc_full"] < 0.5 else "STRONG"
    return (cid, grp, round(r["auc_full"], 3), len(libact), len(tract),
            round(float(sa.mean()), 3), round(float(si.mean()), 3),
            round(sep, 3), round(medmaxT, 3), PANEL[cid].get("name", "")[:20])

if __name__ == "__main__":
    if len(sys.argv) < 2:
        rev = sorted([c for c, r in RES.items() if r["auc_full"] < 0.5], key=lambda c: RES[c]["auc_full"])[:4]
        hi = sorted([c for c, r in RES.items() if r["auc_full"] > 0.7], key=lambda c: -RES[c]["auc_full"])[:4]
        print("SAMPLE_CIDS=" + " ".join(rev + hi))
    else:
        cid = sys.argv[1]
        try:
            row = analyze(cid)
            line = "\t".join(str(x) for x in row)
            with AUC2_ROWS.open("a", encoding="utf-8") as f:
                f.write(line + "\n")
            print("OK", line)
        except Exception as e:
            print(f"FAIL {cid}: {e}")

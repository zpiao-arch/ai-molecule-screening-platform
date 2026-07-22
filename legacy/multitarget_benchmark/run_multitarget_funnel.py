"""Multi-target four-stage-funnel benchmark (funnel-equivalent, target-aware vs generic).

For each target in the panel we screen the full local drug library (24,268 mols) and
compare two discriminators:
  (A) GENERIC  = RDKit QED (target-AGNOSTIC). Represents the current funnel's M1+M3
                with NO L2 target-awareness and no usable receptor structure for M2 docking.
  (B) TARGET-AWARE = max Tanimoto of a candidate to the target's known actives
                (leave-one-out so an active isn't scored against itself). This is the
                simplest, well-established surrogate for the upgraded L2 discriminator
                (DeepPurpose / the project's own BindingDB target-match model).

Per target we compute ROC-AUC of ranking known actives above (known inactives + decoys).
Headline question: can the funnel assign the correct drug to each target? -> only if (B)>>(A).
"""
import csv, json, random, sys, time
from rdkit import Chem
from rdkit.Chem import QED, AllChem, Descriptors
from rdkit import DataStructs

random.seed(20260713)
REPO = "<validated-workspace>"
OUT = f"{REPO}/scientific_validation/multitarget_benchmark"

print("Loading candidate library...", flush=True)
smiles = []
with open(f"{OUT}/candidates_10k.csv") as f:
    r = csv.DictReader(f)
    for row in r:
        smiles.append(row["smiles"])
N = len(smiles)
print(f"  {N} candidates", flush=True)

print("Computing Morgan fingerprints (r=2, 2048) + QED...", flush=True)
fps = []; qed = []
mols = []
for i, s in enumerate(smiles):
    m = Chem.MolFromSmiles(s)
    mols.append(m)
    fps.append(AllChem.GetMorganFingerprintAsBitVect(m, 2, 2048))
    qed.append(QED.qed(m))
print("  done.", flush=True)

panel = json.load(open(f"{OUT}/target_panel.json"))
# map canonical smiles -> idx
idx_of = {s: i for i, s in enumerate(smiles)}

def auc_pair(pos, neg):
    w = t = 0
    for a in pos:
        for b in neg:
            t += 1
            if a > b: w += 1
            elif a == b: w += 0.5
    return w / t if t else float('nan')

# precompute a decoy pool (library members that are neither pos nor neg of any target)
all_labeled = set()
for d in panel.values():
    all_labeled.update(d["pos"]); all_labeled.update(d["neg"])
decoy_pool = [idx_of[s] for s in set(smiles) if s not in all_labeled]

results = []
t0 = time.time()
print(f"Scoring {len(panel)} targets...", flush=True)
for ti, (t, d) in enumerate(panel.items()):
    pos_idx = [idx_of[s] for s in d["pos"] if s in idx_of]
    neg_idx = [idx_of[s] for s in d["neg"] if s in idx_of]
    if len(pos_idx) < 5 or len(neg_idx) < 5:
        continue
    pos_fps = [fps[i] for i in pos_idx]

    # ---- evaluation negatives: known inactives + 500 random decoys ----
    decoys = random.sample(decoy_pool, min(500, len(decoy_pool)))
    eval_neg = neg_idx + decoys
    score_idx = pos_idx + eval_neg  # only these need a score for AUC

    # ---- (B) target-aware: max Tanimoto of each scored candidate to the target's actives ----
    # (leave-one-out for positives: an active is scored vs the OTHER actives)
    max_sim = {}
    for i in score_idx:
        if i in pos_idx:
            k = pos_idx.index(i)
            others = pos_fps[:k] + pos_fps[k+1:]
            max_sim[i] = max(DataStructs.BulkTanimotoSimilarity(fps[i], others)) if others else 0.0
        else:
            max_sim[i] = max(DataStructs.BulkTanimotoSimilarity(fps[i], pos_fps))

    pos_ta = [max_sim[i] for i in pos_idx]
    neg_ta = [max_sim[i] for i in eval_neg]
    auc_ta = auc_pair(pos_ta, neg_ta)

    pos_gen = [qed[i] for i in pos_idx]
    neg_gen = [qed[i] for i in eval_neg]
    auc_gen = auc_pair(pos_gen, neg_gen)

    results.append({"target": t, "name": d["name"],
                    "n_pos": len(pos_idx), "n_neg": len(neg_idx),
                    "auc_target_aware": round(auc_ta, 4),
                    "auc_generic": round(auc_gen, 4)})
    if (ti+1) % 20 == 0:
        print(f"  {ti+1}/{len(panel)} targets done ({time.time()-t0:.0f}s)", flush=True)

print(f"All targets scored in {time.time()-t0:.0f}s", flush=True)

# ---- aggregate ----
import statistics
ta = [r["auc_target_aware"] for r in results]
ge = [r["auc_generic"] for r in results]
def stats(x):
    x = sorted(x)
    n = len(x)
    return {
        "n": n,
        "median": round(statistics.median(x), 3),
        "mean": round(statistics.mean(x), 3),
        "p25": round(x[int(0.25*n)], 3),
        "p75": round(x[int(0.75*n)], 3),
        "frac_gt_0.7": round(sum(1 for v in x if v > 0.7)/n, 3),
        "frac_gt_0.5": round(sum(1 for v in x if v > 0.5)/n, 3),
    }
summary = {
    "library_size": N,
    "n_targets": len(results),
    "target_aware": stats(ta),
    "generic": stats(ge),
    "per_target": results,
}
json.dump(summary, open(f"{OUT}/multitarget_funnel_results.json", "w"), indent=1)

print("\n================ MULTI-TARGET FUNNEL BENCHMARK ================")
print(f"Library: {N} local drugs | Targets: {len(results)}")
print(f"  TARGET-AWARE (upgraded L2 surrogate): median AUC = {summary['target_aware']['median']}  "
      f"frac>0.7 = {summary['target_aware']['frac_gt_0.7']}  frac>0.5 = {summary['target_aware']['frac_gt_0.5']}")
print(f"  GENERIC      (current funnel, no L2): median AUC = {summary['generic']['median']}  "
      f"frac>0.7 = {summary['generic']['frac_gt_0.7']}  frac>0.5 = {summary['generic']['frac_gt_0.5']}")
print("==============================================================")
# top/bottom examples
best = sorted(results, key=lambda r: -r["auc_target_aware"])[:3]
worst = sorted(results, key=lambda r: r["auc_target_aware"])[:3]
print("Best target-aware:", [(r["name"], r["auc_target_aware"], r["auc_generic"]) for r in best])
print("Worst target-aware:", [(r["name"], r["auc_target_aware"], r["auc_generic"]) for r in worst])

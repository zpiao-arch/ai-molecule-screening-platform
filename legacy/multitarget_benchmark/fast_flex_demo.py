# 精简版 NA 柔性对接演示: 少柔性残基 + 低 exhaustiveness, 快速出数
import os, json, tempfile, subprocess
from pathlib import Path
import flexible_docking as fdmod

root = Path("<validated-workspace>")
rec_pdbqt = str(root / "ai_mol_loop/influenza_na_2000_project_20260625/stage4/receptors/3TI6_protein_only_obabel.pdbqt")
center = (-28.914, 14.334, 20.794); size = (23.585, 20.45, 24.18)
flexres = "A:118,A:119,A:151,A:152"   # 4 个活性腔核心残基 (精简)
ligands = {
    "oseltamivir_carboxylate(NA活性)": "CC(C)O[C@@H]1C=C[C@H](OC(=O)C)C(=O)N1C",
    "caffeine(非NA对照)": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
}
fd = fdmod.FlexibleDocking()
print(f"受体存在={os.path.exists(rec_pdbqt)}")
fd.boltz_multiconf("x", "x"); fd.mcce_protonate("x", "x")
results = []
tmp = tempfile.mkdtemp(prefix="flexfast_")
for name, smi in ligands.items():
    lig = os.path.join(tmp, name.split("(")[0] + ".pdbqt")
    if not fd.smiles_to_pdbqt(smi, lig):
        print(f"  跳过 {name}"); continue
    rr = fd.dock(rec_pdbqt, lig, center, size, flexres="", exhaustiveness=4)
    rf = fd.dock(rec_pdbqt, lig, center, size, flexres=flexres, exhaustiveness=4)
    print(f"  {name:32s} rigid={rr['affinity_kcal_mol']}  flex={rf['affinity_kcal_mol']}")
    results.append({"ligand": name, "rigid_aff": rr["affinity_kcal_mol"], "flex_aff": rf["affinity_kcal_mol"]})
out = root / "scientific_validation/multitarget_benchmark/na_flexible_docking_demo.json"
json.dump(results, open(out, "w"), indent=1, ensure_ascii=False)
print("saved", out)

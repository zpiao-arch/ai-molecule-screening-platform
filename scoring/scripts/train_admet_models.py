"""训练DeepChem ADMET模型并保存为pickle文件"""
import pickle
import numpy as np
from pathlib import Path
from sklearn.ensemble import RandomForestClassifier
from deepchem.molnet import load_tox21, load_bbbp, load_clintox, load_sider
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

MODEL_DIR = Path(__file__).parent / "models" / "admet"
MODEL_DIR.mkdir(exist_ok=True)


def featurize_smiles(smiles_list):
    """SMILES → Morgan指纹 (1024-bit)"""
    fpg = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
    fps = []
    for s in smiles_list:
        mol = Chem.MolFromSmiles(s)
        if mol is None:
            fps.append(np.zeros(1024, dtype=np.float32))
        else:
            fps.append(np.array(fpg.GetFingerprintAsNumPy(mol), dtype=np.float32))
    return np.array(fps)


def train_and_save(name, loader_fn):
    print(f"\n{'='*50}")
    print(f"训练 {name} 模型...")

    tasks, datasets, transformers = loader_fn()
    train_data = datasets[0]
    smiles_list = train_data.ids
    X = featurize_smiles(smiles_list)
    y = train_data.y.astype(np.float32)

    print(f"  样本数: {len(smiles_list)}, 任务数: {len(tasks)}")
    print(f"  任务: {tasks[:5]}...")

    model = RandomForestClassifier(
        n_estimators=100, max_depth=15, random_state=42,
        n_jobs=-1, class_weight='balanced'
    )
    model.fit(X, y)

    # 计算训练准确率
    train_acc = model.score(X, y)
    print(f"  训练准确率: {train_acc:.4f}")

    # 保存
    model_path = MODEL_DIR / f"{name}.pkl"
    meta = {"name": name, "tasks": tasks, "n_samples": len(smiles_list)}
    with open(model_path, "wb") as f:
        pickle.dump({"model": model, "meta": meta}, f)
    print(f"  已保存: {model_path} ({model_path.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    train_and_save("tox21", load_tox21)
    train_and_save("bbbp", load_bbbp)
    train_and_save("clintox", load_clintox)
    train_and_save("sider", load_sider)
    print(f"\n全部模型已保存到: {MODEL_DIR}")

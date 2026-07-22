#!/usr/bin/env python3
# 探究: 是否需要"多模型 + 选择树"?
#   1) 算 oracle 天花板 = 逐靶 max(hash_auc, ce_auc) 的中位  -> 路由的理论上限
#   2) 看 CE 到底在哪些靶上赢了 hash (赢了多少个 / 是哪类靶)
#   3) 测试"按 n_pos 阈值路由"能恢复多少 oracle  (n_pos 是推理期可用的廉价信号)
# 用法: .venv_mlx_qwen35/bin/python3 analyze_model_selection.py
import json, numpy as np

H = json.load(open("bigrun_results_hash.json"))["per_target_metrics"]
C = json.load(open("bigrun_results_retrain471.json"))["per_target_metrics"]

rows = []  # (target, hash_auc, ce_auc, n_pos)
for t, v in H.items():
    ha = v.get("auc"); np_ = v.get("n_pos")
    if ha is None: continue
    ca = C.get(t, {}).get("auc")
    rows.append((t, float(ha), (float(ca) if ca is not None else np.nan), np_ or 0))
rows = [r for r in rows if not (np.isnan(r[1]) or (not np.isnan(r[2]) and np.isnan(r[2])))]
print(f"共同靶数: {len(rows)}\n")

ha = np.array([r[1] for r in rows])
ca = np.array([r[2] for r in rows])
np_ = np.array([r[3] for r in rows])

def med(a): return float(np.nanmedian(a))
def rev(a): return float(np.mean(a < 0.5) * 100)

print("=== 单模型基线 (471 标注靶) ===")
print(f"  哈希基线  中位AUC={med(ha):.4f}  反向率={rev(ha):.1f}%  AUC>=0.7={int((ha>=0.7).sum())}")
print(f"  CE重训    中位AUC={med(ca):.4f}  反向率={rev(ca):.1f}%  AUC>=0.7={int((ca>=0.7).sum())}")

# 1) oracle 天花板
oracle = np.maximum(ha, ca)
print(f"\n=== Oracle 天花板 (逐靶 max, 假设有完美路由器) ===")
print(f"  Oracle 中位AUC={med(oracle):.4f}  反向率={rev(oracle):.1f}%  AUC>=0.7={int((oracle>=0.7).sum())}")
print(f"  -> 路由相对哈希基线增益: {med(oracle)-med(ha):+.4f}")

# 2) CE 在哪些靶赢了
ce_win = ca > ha + 0.02
hash_win = ha > ca + 0.02
print(f"\n=== CE vs 哈希 胜负 ===")
print(f"  CE 显著更优 (Δ>0.02): {int(ce_win.sum())} 个靶")
print(f"  哈希显著更优 (Δ>0.02): {int(hash_win.sum())} 个靶")
print(f"  大致持平: {len(rows)-int(ce_win.sum())-int(hash_win.sum())} 个靶")

# CE 赢的靶的 n_pos 分布 vs 哈希赢的
print(f"\n  CE赢的靶  n_pos中位={np.median(np_[ce_win]):.0f} 均值={np.mean(np_[ce_win]):.1f}")
print(f"  哈希赢的靶 n_pos中位={np.median(np_[hash_win]):.0f} 均值={np.mean(np_[hash_win]):.1f}")

# 3) n_pos 阈值路由测试
print(f"\n=== 按 n_pos 阈值路由: 'n_pos<T 用CE, 否则用哈希' 的恢复率 ===")
best_single = med(ha)
for T in [5, 10, 15, 20, 30]:
    chosen = np.where(np_ < T, ca, ha)
    m = med(chosen)
    rec = (m - best_single) / (med(oracle) - best_single) * 100 if med(oracle) != best_single else 0
    print(f"  T={T:>2}: 中位AUC={m:.4f}  恢复oracle比例={rec:.0f}%  (路由错配靶={int(((np_<T)&hash_win).sum())+int(((np_>=T)&ce_win).sum())})")

# 4) 关键反问: 哈希在"良采样"靶已很好, CE 是否只救了"本就差的靶"却坑了"本就好的"?
good_hash = ha >= 0.7
print(f"\n=== CE 对好靶(哈希>=0.7)的伤害 ===")
print(f"  哈希好靶中, CE 反而<0.5 的有 {int((good_hash & (ca<0.5)).sum())} 个 (被坑)")
print(f"  哈希差靶(哈希<0.5)中, CE 救回>0.7 的有 {int(((ha<0.5) & (ca>=0.7)).sum())} 个")

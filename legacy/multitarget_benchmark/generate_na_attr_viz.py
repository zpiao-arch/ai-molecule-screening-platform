#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""从 diagnose_na_recall.json 生成两张归因证据图(SVG)."""
import json, os

D = json.load(open("<external-library-cache>/diagnose_na_recall.json"))
OUT = "<validated-workspace>/scientific_validation/multitarget_benchmark"

# ── 图1: L2 分数全库分布(塌缩) ──
W, H, L, T, R, B = 680, 340, 55, 30, 15, 60
pw, ph = W - L - R, H - T - B
edges = D["hist_edges"]; counts = D["hist_counts"]
maxc = max(counts)
def col(i):
    c = (i + 0.5) / 20
    if c >= 0.9: return "#E24B4A"
    if c >= 0.8: return "#EF9F27"
    if c >= 0.5: return "#85B7EB"
    return "#378ADD"
bw = pw / 20
bars = ""
for i, c in enumerate(counts):
    h = ph * c / maxc
    x = L + i * bw
    bars += f'<rect x="{x:.1f}" y="{T+ph-h:.1f}" width="{bw-2:.1f}" height="{h:.1f}" fill="{col(i)}"/>'
med = D["l2_median"]; mx = L + med * pw
medline = f'<line x1="{mx:.1f}" y1="{T-6}" x2="{mx:.1f}" y2="{T+ph}" stroke="#3B6D11" stroke-width="1.5" stroke-dasharray="5 4"/>'
xt = ""
for v in [0, 0.2, 0.4, 0.6, 0.8, 1.0]:
    x = L + v * pw
    xt += f'<text x="{x:.1f}" y="{T+ph+18}" font-size="11" fill="#5F5E5A" text-anchor="middle">{v:.1f}</text>'
svg1 = f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" role="img">
<title>NA靶点 L2 分数全库分布坍缩</title>
<desc>100k分子L2分数直方图: {D["frac_gt_0_90"]*100:.0f}%分子分数&gt;0.9, {D["frac_gt_0_95"]*100:.0f}%在0.95-1.0; 18活性中位0.888与背景中位{med:.3f}几乎重合</desc>
<line x1="{L}" y1="{T+ph}" x2="{L+pw}" y2="{T+ph}" stroke="#888780" stroke-width="1"/>
{bars}
{medline}
{xt}
<text x="{L-8}" y="{T+10}" font-size="11" fill="#5F5E5A" text-anchor="end">频数</text>
<text x="{mx:.1f}" y="{T-10}" font-size="11" fill="#3B6D11" text-anchor="middle">背景中位 {med:.3f}</text>
<text x="{L+pw}" y="{T+ph+40}" font-size="10" fill="#A32D2D" text-anchor="end">红区=分数&gt;0.9: {D["frac_gt_0_90"]*100:.0f}%分子(47k) — 活性被淹没其中</text>
</svg>'''
open(f"{OUT}/diagnose_na_l2_hist.svg", "w").write(svg1)

# ── 图2: 18 活性在全库的 L2 排名(散布=无排序能力) ──
W2, H2, L2_, T2, R2, B2 = 680, 300, 55, 30, 15, 55
pw2, ph2 = W2 - L2_ - R2, H2 - T2 - B2
ranks = D["active_ranks"]; N = D["library_size"]
# y: 排名(线性, 0=榜首, N=榜尾)
dots = ""
for rk in ranks:
    x = L2_ + (rk / N) * pw2
    dots += f'<circle cx="{x:.1f}" cy="{T2+ph2/2:.1f}" r="4" fill="#A32D2D" stroke="#fff" stroke-width="0.7"/>'
# top-N 窗口阴影 (250/500/1000/5000)
shade = ""
for topN, col2 in [(250, "#F3C9C9"), (500, "#E9D6A8"), (1000, "#CFE0F2"), (5000, "#BFD8EE")]:
    w = (topN / N) * pw2
    shade += f'<rect x="{L2_:.1f}" y="{T2:.1f}" width="{w:.1f}" height="{ph2:.1f}" fill="{col2}" opacity="0.55"/>'
xt2 = ""
for p in [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
    x = L2_ + p * pw2
    xt2 += f'<text x="{x:.1f}" y="{T2+ph2+18}" font-size="10" fill="#5F5E5A" text-anchor="middle">{int(p*N/1000)}k</text>'
labels = (f'若 L2 top-250 接对接: 落入窗口活性 = {D["topN_window_recall"]["250"]}/18   '
          f'top-1000 = {D["topN_window_recall"]["1000"]}/18   top-5000 = {D["topN_window_recall"]["5000"]}/18')
svg2 = f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W2} {H2}" role="img">
<title>NA 18个活性在100k库中的L2排名散布</title>
<desc>18活性散布于5665~92602名, 几乎贯穿全库, 说明L2对NA检索近乎随机排序</desc>
{shade}
<line x1="{L2_}" y1="{T2}" x2="{L2_}" y2="{T2+ph2}" stroke="#888780" stroke-width="1"/>
<line x1="{L2_}" y1="{T2+ph2}" x2="{L2_+pw2}" y2="{T2+ph2}" stroke="#888780" stroke-width="1"/>
{dots}
{xt2}
<text x="{L2_}" y="{T2-10}" font-size="11" fill="#5F5E5A">L2排名轴(榜首→榜尾)</text>
<text x="{L2_+pw2}" y="{T2+ph2+40}" font-size="10" fill="#A32D2D" text-anchor="end">{labels}</text>
<text x="{L2_}" y="{T2+ph2+40}" font-size="10" fill="#3B6D11">浅蓝=top-5000窗口, 红点=18活性(中位排名{D["active_rank_median"]})</text>
</svg>'''
open(f"{OUT}/diagnose_na_active_ranks.svg", "w").write(svg2)
print("SVG1/2 written")

#!/usr/bin/env python3
"""
Rearranged Figure 6 (main) + Supplementary Figure 6 for the PPARG/PNPLA3 manuscript.

Inputs (same folder or UP/OUT below):
  figure_6.pdf (original; panels a-d are the authors' luciferase/EMSA artwork and are re-used as-is)
  PNPLA3_signal_profile.csv, FABP4_signal_profile.csv
  PNPLA3_/FABP4_quantitative_signal_metrics.csv, combined_8gene_heatmap_data.csv
  all_peaks_annotated.csv, ASCRosi_DE_with_ChIP_annotation.csv, ChIP_bound_genes_Hallmark_ORA.csv
Outputs:
  MAIN_Figure6_rearranged.pdf, Supplementary_Figure6_PPARG_ChIP.pdf
"""
import subprocess, os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch
from scipy.stats import mannwhitneyu, fisher_exact

import tempfile
# All input files (see header) and outputs live in ONE folder: set PPARG_DIR or run from that folder.
UP = OUT = os.environ.get("PPARG_DIR", ".")
TMP = tempfile.mkdtemp(prefix="fig6_")
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 6.5, "axes.linewidth": 0.6,
                     "pdf.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False})
HASC, TIS = "#5B4A94", "#E0A526"
SAMP = ["GSM3177950", "GSM3177951", "GSM3177952", "GSM3177953", "GSM3177954",
        "GSM3416379", "GSM3416380", "GSM3416381", "GSM3416382", "GSM3416383"]
col = lambda s: HASC if s.startswith("GSM31779") else TIS
TSS = {"PNPLA3": 43923805, "FABP4": 81483350}

def label(ax, t, x=-0.12, y=1.04):
    ax.text(x, y, t, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")

# ---------- panels a-d: crop from the authors' original artwork (300 dpi) ----------
subprocess.run(["pdftoppm", "-r", "300", "-png", "-f", "1", "-l", "1", f"{UP}/figure_6.pdf", f"{TMP}/orig"], check=True)
from PIL import Image
Image.MAX_IMAGE_PIXELS = None
im = Image.open(f"{TMP}/orig-1.png")
W, H = im.size
top = im.crop((0, int(H*0.005), W, int(H * 0.1895)))          # a-d occupy the top ~19% of the page
bb = Image.eval(top.convert("L"), lambda v: 255 if v < 245 else 0).getbbox()
top = top.crop((max(bb[0]-20,0), max(bb[1]-20,0), min(bb[2]+20,top.size[0]), min(bb[3]+20,top.size[1])))
ASPECT = top.size[0] / top.size[1]
top.save(f"{TMP}/panels_a_d.png")

# ---------- data ----------
def load_prof(g):
    d = pd.read_csv(f"{UP}/{g}_signal_profile.csv"); d["rel"] = (d.position - TSS[g]) / 1000.0; return d
prof = {g: load_prof(g) for g in TSS}
mp = pd.read_csv(f"{OUT}/PNPLA3_quantitative_signal_metrics.csv").assign(gene="PNPLA3")
mf = pd.read_csv(f"{OUT}/FABP4_quantitative_signal_metrics.csv").assign(gene="FABP4")
met = pd.concat([mf, mp])
de = pd.read_csv(f"{OUT}/ASCRosi_DE_with_ChIP_annotation.csv")
heat = pd.read_csv(f"{OUT}/combined_8gene_heatmap_data.csv")
peaks = pd.read_csv(f"{UP}/all_peaks_annotated.csv")

# PPRE (EMSA/luciferase) positions relative to PNPLA3 TSS
PPRE = [(-0.790, -0.780), (-0.713, -0.707)]
PEAK = (43925319, 43925915)

def tracks(fig, spec, g, ttl, tag, annotate):
    sub = spec.subgridspec(10, 1, hspace=0.0)
    d = prof[g]; ymax = d[SAMP].values.max() * 1.05
    axs = []
    for i, s in enumerate(SAMP):
        ax = fig.add_subplot(sub[i]); axs.append(ax)
        ax.fill_between(d.rel, 0, d[s], color=col(s), lw=0, step="mid")
        ax.axvline(0, color="firebrick", ls="--", lw=0.6)
        if g == "PNPLA3":
            for a, b in PPRE: ax.axvspan(a - 0.04, b + 0.04, color="#2E86C1", alpha=0.45, lw=0)
            ax.axvspan((PEAK[0]-TSS[g])/1000, (PEAK[1]-TSS[g])/1000, color="#C0392B", alpha=0.12, lw=0)
        ax.set_ylim(0, ymax); ax.set_xlim(d.rel.min(), d.rel.max())
        ax.set_yticks([]); ax.patch.set_alpha(0)
        for sp in ["left"]: ax.spines[sp].set_visible(False)
        if i < 9: ax.spines["bottom"].set_linewidth(0.4); ax.set_xticks([])
        ax.text(0.005, 0.55, s.replace("GSM", ""), transform=ax.transAxes, fontsize=5.5, va="center")
        ax.text(0.995, 0.15, "hASC" if i < 5 else "tissue", transform=ax.transAxes, ha="right", fontsize=5.5, color=col(s))
    axs[0].set_title(ttl, fontsize=7, pad=52)
    axs[-1].set_xlabel("distance from TSS (kb)")
    axs[0].text(0, 1.05, "TSS", transform=axs[0].get_xaxis_transform(), ha="center", color="firebrick", fontsize=5.5)
    if annotate:
        a0 = axs[0]; tr = a0.get_xaxis_transform()
        a0.annotate("EMSA/luciferase PPREs\n(−0.79/−0.71 kb)", xy=(-0.75, 1.0), xycoords=tr, xytext=(-1.9, 1.3), textcoords=tr,
                    fontsize=5.3, color="#1F618D", ha="left", va="bottom", annotation_clip=False,
                    arrowprops=dict(arrowstyle="-", color="#1F618D", lw=0.5))
        a0.annotate("ChIP peak (intron 1)\n+1.5 to +2.1 kb", xy=(1.8, 1.0), xycoords=tr, xytext=(2.3, 1.3), textcoords=tr,
                    fontsize=5.3, color="#922B21", ha="left", va="bottom", annotation_clip=False,
                    arrowprops=dict(arrowstyle="-", color="#922B21", lw=0.5))
    return axs

# ================= MAIN FIGURE =================
fig = plt.figure(figsize=(7.2, (6.6/ASPECT + 4.3 + 2.3) / 0.84))
gs = GridSpec(3, 6, figure=fig, height_ratios=[6.6/ASPECT, 4.3, 2.3], hspace=0.32, wspace=1.0,
              left=0.07, right=0.985, top=0.985, bottom=0.05)
axtop = fig.add_subplot(gs[0, :]); axtop.imshow(np.asarray(Image.open(f"{TMP}/panels_a_d.png"))); axtop.axis("off")

axs_e = tracks(fig, gs[1, 0:3], "PNPLA3", "PPARG ChIP-seq signal at PNPLA3 (10 samples)", "e", True)
axs_f = tracks(fig, gs[1, 3:6], "FABP4", "PPARG ChIP-seq signal at FABP4 (10 samples)", "f", False)
axs_e[0].text(-0.05, 3.3, "e", transform=axs_e[0].transAxes, fontsize=10, fontweight="bold")
axs_f[0].text(-0.05, 3.3, "f", transform=axs_f[0].transAxes, fontsize=10, fontweight="bold")

# g quantitative occupancy
ag = fig.add_subplot(gs[2, 0:2]); label(ag, "g", -0.28)
rng = np.random.default_rng(1)
for xi, gene in enumerate(["FABP4", "PNPLA3"]):
    for k, (grp, c) in enumerate([("hASC", HASC), ("tissue", TIS)]):
        v = met[(met.gene == gene) & (met.group == grp)].mean_signal.values
        x = xi + (k - .5) * 0.5 + rng.uniform(-.06, .06, len(v))
        ag.scatter(x, v, s=14, color=c, edgecolor="k", lw=0.4, zorder=3)
        m, se = v.mean(), v.std(ddof=1) / np.sqrt(len(v))
        ag.errorbar(xi + (k - .5) * 0.5, m, se, color="k", capsize=2, lw=0.7, zorder=4)
        ag.plot([xi + (k - .5) * 0.5 - .12, xi + (k - .5) * 0.5 + .12], [m, m], color="k", lw=0.9, zorder=4)
ag.set_xticks([0, 1]); ag.set_xticklabels(["FABP4", "PNPLA3"]); ag.set_xlim(-.6, 1.6); ag.set_ylim(0)
ag.set_ylabel("Mean signal, TSS ± 2 kb (RPGC)"); ag.set_title("Quantitative occupancy", fontsize=7)
ag.legend(handles=[Patch(color=HASC, label="hASC (n=5)"), Patch(color=TIS, label="tissue (n=5)")], frameon=False, fontsize=5.5, loc="upper right")

# h ECDF of rosiglitazone t-stat
ah = fig.add_subplot(gs[2, 2:4]); label(ah, "h", -0.28)
b, u = de[de.bound].t_stat.values, de[~de.bound].t_stat.values
mw = mannwhitneyu(b, u, alternative="greater")
for v, c, n in [(u, "#8C9A9A", f"not bound (n={len(u):,})"), (b, "#B05FA6", f"PPARG-bound (n={len(b):,})")]:
    x = np.sort(v); ah.plot(x, np.arange(1, len(x) + 1) / len(x), color=c, lw=1.3, label=n)
pn = de[de.Gene == "PNPLA3"].t_stat.iloc[0]
pct = (b < pn).mean() * 100
ah.axvline(pn, color="k", ls=":", lw=0.7); ah.text(pn + .3, 0.12, f"PNPLA3\n(t = {pn:.1f};\ntop {100-pct:.0f}% of bound genes)", fontsize=5.3)
ah.set_xlim(-12, 12); ah.set_xlabel("Rosiglitazone response (paired t, Rosi vs DMSO)"); ah.set_ylabel("Cumulative fraction")
ah.set_title("Response by binding status", fontsize=7); ah.legend(frameon=False, fontsize=5.5, loc="upper left")
ah.text(0.97, 0.05, f"Mann–Whitney P = {mw.pvalue:.1e}", transform=ah.transAxes, ha="right", fontsize=5.5)

# i enrichment bars
ai = fig.add_subplot(gs[2, 4:6]); label(ai, "i", -0.3)
nb, ni = len(b), int(((de.bound) & (de.induced)).sum()); nu, nui = len(u), int(((~de.bound) & (de.induced)).sum())
orr, p = fisher_exact([[ni, nb - ni], [nui, nu - nui]], alternative="greater")
ai.bar([0, 1], [100 * ni / nb, 100 * nui / nu], color=["#B05FA6", "#8C9A9A"], width=0.6)
for x, v, n, N in [(0, 100 * ni / nb, ni, nb), (1, 100 * nui / nu, nui, nu)]:
    ai.text(x, v + 0.08, f"{v:.1f}%\n({n}/{N:,})", ha="center", fontsize=5.8)
ai.set_xticks([0, 1]); ai.set_xticklabels(["PPARG-bound", "Not bound"]); ai.set_ylim(0, 4.8)
ai.set_ylabel("% induced by rosiglitazone (padj < 0.05)")
ai.set_title(f"OR = {orr:.2f}; Fisher P = {p:.1e}", fontsize=7)
fig.savefig(f"{OUT}/MAIN_Figure6_rearranged.pdf"); plt.close(fig)
print("main:", f"OR={orr:.2f} P={p:.2e} MW P={mw.pvalue:.2e} PNPLA3 t={pn:.2f} pct={pct:.0f}")

# ================= SUPPLEMENTARY FIGURE =================
fig = plt.figure(figsize=(7.2, 8.6))
gs = GridSpec(3, 2, figure=fig, height_ratios=[1.0, 1.25, 1.1], hspace=0.55, wspace=0.55,
              left=0.1, right=0.97, top=0.965, bottom=0.07)
# a reproducibility
genes = ["ADIPOQ", "CD36", "FABP4", "PCK1", "PLIN1", "LPL", "SLC2A4", "PPARG", "PNPLA3"]
q = peaks[(peaks.SYMBOL.isin(genes)) & (peaks.distanceToTSS.abs() <= 5000)]
cnt = q.groupby(["SYMBOL", "dataset"])["sample"].nunique().unstack(fill_value=0).reindex(genes)
a = fig.add_subplot(gs[0, 0]); label(a, "a", -0.2)
x = np.arange(len(genes))
a.bar(x - .2, cnt["GSE115421_hASC"], .4, color=HASC, label="hASC"); a.bar(x + .2, cnt["GSE115421_tissue"], .4, color=TIS, label="tissue")
a.set_xticks(x); a.set_xticklabels(genes, rotation=45, ha="right"); a.set_ylabel("Samples with peak (±5 kb TSS)"); a.set_ylim(0, 5.6)
a.axhline(3, color="grey", ls=":", lw=.6); a.legend(frameon=False, fontsize=5.5, ncol=2, loc="upper center")
a.set_title("Binding reproducibility", fontsize=7)
# d ORA
o = pd.read_csv(f"{OUT}/ChIP_bound_genes_Hallmark_ORA.csv").sort_values("padj").head(10)
d_ = fig.add_subplot(gs[0, 1]); label(d_, "b", -0.5)
d_.barh(range(len(o))[::-1], -np.log10(o.padj), color="#B05FA6")
d_.set_yticks(range(len(o))[::-1]); d_.set_yticklabels([t.replace("HALLMARK_", "").replace("_", " ").title() for t in o.Term], fontsize=5.5)
d_.set_xlabel("−log10(adjusted P)"); d_.set_title("Hallmark enrichment of PPARG-bound genes", fontsize=7)
# heatmaps
order = ["ADIPOQ", "CD36", "FABP4", "PLIN1", "PNPLA3", "MYOD1", "RHO", "INS"]
M = heat.pivot(index="gene", columns="sample", values="mean_signal").loc[order, SAMP]
ax = fig.add_subplot(gs[1, 0]); label(ax, "c", -0.35)
imv = ax.imshow(np.log2(M.values + 1), aspect="auto", cmap="viridis")
ax.set_yticks(range(8)); ax.set_yticklabels(order); ax.set_xticks(range(10)); ax.set_xticklabels([s.replace("GSM", "") for s in SAMP], rotation=90, fontsize=5.5)
ax.axhline(4.5, color="w", lw=1); ax.axhline(4 - .5, color="w", lw=1)
ax.set_title("Raw signal, log2(mean RPGC + 1), TSS ± 2 kb", fontsize=7); plt.colorbar(imv, ax=ax, fraction=.04, pad=.02)
ax.axvline(4.5, color="w", lw=1)
Z = (M.values - M.values.mean(1, keepdims=True)) / M.values.std(1, keepdims=True)
ax2 = fig.add_subplot(gs[1, 1]); label(ax2, "d", -0.35)
im2 = ax2.imshow(Z, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
ax2.set_yticks(range(8)); ax2.set_yticklabels(order); ax2.set_xticks(range(10)); ax2.set_xticklabels([s.replace("GSM", "") for s in SAMP], rotation=90, fontsize=5.5)
ax2.axvline(4.5, color="k", lw=.8); ax2.set_title("Within-gene z-score (pattern only)", fontsize=7); plt.colorbar(im2, ax=ax2, fraction=.04, pad=.02)
# e FABP4/PNPLA3 window sanity: per-sample dot of peak vs upstream window
up = prof["PNPLA3"]
def win(lo, hi): return up[(up.position >= lo) & (up.position <= hi)][SAMP].mean()
r = pd.DataFrame({"Upstream PPRE window": win(43922950, 43923150), "Upstream background": win(43921800, 43922600),
                  "Intron-1 ChIP peak": win(PEAK[0], PEAK[1])})
e = fig.add_subplot(gs[2, :]); label(e, "e", -0.05)
for i, s in enumerate(SAMP):
    e.plot(range(3), r.loc[s].values, "-o", color=col(s), ms=3, lw=.6, alpha=.8)
e.set_xticks(range(3)); e.set_xticklabels(r.columns); e.set_ylabel("Mean signal (RPGC)")
e.set_title("PPARG signal at the reporter-validated upstream PPREs vs the intron-1 peak (PNPLA3; each line = one sample)", fontsize=7)
e.text(0.01, 0.9, "PPRE window coordinates are estimated (chr22:43,922,950–43,923,150); confirm against construct.",
       transform=e.transAxes, fontsize=5.5, color="#555")
e.legend(handles=[Patch(color=HASC, label="hASC"), Patch(color=TIS, label="tissue")], frameon=False, fontsize=5.5, loc="upper left", bbox_to_anchor=(0, .8))
fig.savefig(f"{OUT}/Supplementary_Figure6_PPARG_ChIP.pdf"); plt.close(fig)
print(r.mean().round(2).to_dict())

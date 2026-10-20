#!/usr/bin/env python3
"""Manuscript figures for the HepG2 PPARG ChIP-seq analysis: one MAIN figure and one SUPPLEMENTARY figure.

Reads the outputs of steps 1-5 (tracks, peaks, full_analysis/*.csv) - nothing is recomputed except the TZD-response
summary (from the DE tables and the saved bound-gene lists) - and writes, in $OUT/full_analysis/manuscript_figures/:
  HepG2_main_figure.pdf / .png          (183 mm wide)
  HepG2_supplementary_figure.pdf / .png (183 mm wide)
  source_data_*.csv                     (plotted values, one file per panel group)

Environment variables: OUT, GTF, DE_CSV, HASC_CSV  (same as step 5)
Needs: pyBigWig pandas numpy scipy matplotlib
"""
import os, re, gzip, glob, sys
import numpy as np, pandas as pd, pyBigWig
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.patches import Rectangle, FancyArrowPatch

OUT = os.environ.get("OUT", "results/hepg2_ppar"); GTF = os.environ["GTF"]
DE_CSV, HASC_CSV = os.environ.get("DE_CSV", ""), os.environ.get("HASC_CSV", "")
FA = f"{OUT}/full_analysis"; FD = f"{FA}/manuscript_figures"; os.makedirs(FD, exist_ok=True)
MM = 1 / 25.4
ML = os.environ.get("MAIN_LETTERS", "abcd")   # panel letters of the MAIN figure, e.g. MAIN_LETTERS=fghi when it extends an existing figure
RS_POS = int(os.environ.get("RS_POS", 43928847))   # rs738409 (PNPLA3 I148M), GRCh38 chr22; checked against the exon model below
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"], "font.size": 6.5,
                     "axes.labelsize": 6.5, "axes.titlesize": 7, "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 5.5,
                     "axes.linewidth": .6, "xtick.major.width": .6, "ytick.major.width": .6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
                     "pdf.fonttype": 42, "ps.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False})
C = {"rep1": "#0072B2", "rep2": "#D55E00", "pooled": "#222222", "input": "#8c8c8c", "bound": "#7b3f98", "unbound": "#b0b0b0",
     "site_up": "#d62728", "site_in": "#2ca02c", "site_48": "#ff9f1c"}
SITE_COL = {"upstream_PPRE": C["site_up"], "intron1_site": C["site_in"], "site_plus4.8kb": C["site_48"]}
SITE_LAB = {"upstream_PPRE": "upstream PPREs\n(−850/−650)", "intron1_site": "intron-1 site\n(+1.0/+2.3 kb)", "site_plus4.8kb": "+4.8 kb peak"}

# ---------------------------------------------------------------- inputs
def need(f):
    if not os.path.exists(f): sys.exit(f"Missing {f}. Run steps 1-5 first.")
    return f
reps = sorted(os.path.basename(p).replace(".RPGC.bw", "") for p in glob.glob(f"{OUT}/tracks/*.RPGC.bw")
              if not os.path.basename(p).startswith(("control", "ChIP_pooled"))
              and os.path.exists(f"{OUT}/peaks/" + os.path.basename(p).replace(".RPGC.bw", "") + "_peaks.narrowPeak"))
if len(reps) < 2: sys.exit("need two replicate tracks")
S = reps + ["ChIP_pooled"]; SC = dict(zip(S, [C["rep1"], C["rep2"], C["pooled"]]))
SL = {reps[0]: reps[0], reps[1]: reps[1], "ChIP_pooled": "pooled"}
trk = {s: pyBigWig.open(need(f"{OUT}/tracks/{s}.log2ChIPoverInput.bw")) for s in S}
ctrl = pyBigWig.open(need(f"{OUT}/tracks/control.RPGC.bw"))
def rp(f):
    try: d = pd.read_csv(f, sep="\t", header=None, usecols=[0, 1, 2, 6, 7, 8, 9], names=["chr", "s", "e", "fc", "nlp", "nlq", "summit"])
    except pd.errors.EmptyDataError: d = pd.DataFrame(columns=["chr", "s", "e", "fc", "nlp", "nlq", "summit"])
    return d
P = {(s, "q01"): rp(need(f"{OUT}/peaks/{s}_peaks.narrowPeak")) for s in S}
P.update({(s, "q05"): rp(need(f"{OUT}/peaks/relaxed/{s}_q05_peaks.narrowPeak")) for s in S})
W = pd.read_csv(need(f"{FA}/window_enrichment_tests.csv"))
QC = pd.read_csv(need(f"{FA}/peak_QC_annotation.csv")); CONC = pd.read_csv(need(f"{FA}/replicate_peak_concordance.csv"))
HALL = pd.read_csv(f"{FA}/Hallmark_ORA_HepG2_bound_genes.csv") if os.path.exists(f"{FA}/Hallmark_ORA_HepG2_bound_genes.csv") else None
OVL = pd.read_csv(f"{FA}/overlap_with_hASC_bound_genes.csv") if os.path.exists(f"{FA}/overlap_with_hASC_bound_genes.csv") else None
MOT = pd.read_csv(f"{FA}/motif_scan_named_sites.csv") if os.path.exists(f"{FA}/motif_scan_named_sites.csv") else None
REV = f"{FA}/revalidation"
def rd(f): return pd.read_csv(f) if os.path.exists(f) else None
PI, PX, UNI = rd(f"{REV}/panel_i_statistics.csv"), rd(f"{REV}/panel_i_interaction.csv"), rd(f"{REV}/common_universe_table.csv")
if PI is None: print("NOTE: step 7 (revalidation) outputs not found - panel d/i falls back to the symbol-matched analysis. Run scripts/chipseq_hepg2/07_hepg2_revalidation.py first.")
def macs_hit(gene, win):
    """best MACS3 peak (pooled; q<0.01 else q<0.05) overlapping a named window -> (-log10 q, fold change, level) or None"""
    r = W[(W.gene == gene) & (W.window == win)].iloc[0]
    for lvl in ("q01", "q05"):
        d = P[("ChIP_pooled", lvl)]; d = d[(d.chr == r.chr) & (d.s < r.end) & (d.e > r.start)]
        if len(d): k = d.nlq.astype(float).idxmax(); return float(d.loc[k, "nlq"]), float(d.loc[k, "fc"]), lvl
    return None
INT = pd.read_csv(f"{FA}/integration_with_RNAseq.csv") if os.path.exists(f"{FA}/integration_with_RNAseq.csv") else None

# ---------------------------------------------------------------- annotation (TSS, strand, exon model)
GENES = ["PNPLA3", "PLIN1", "PCK1", "FABP4"]
op = gzip.open if GTF.endswith(".gz") else open
tx, ex = {}, {}
with op(GTF, "rt") as f:
    for ln in f:
        if ln[0] == "#": continue
        c = ln.rstrip("\n").split("\t")
        if len(c) < 9 or c[2] not in ("transcript", "exon"): continue
        m = re.search(r'gene_name "([^"]+)"', c[8])
        if not m or m.group(1) not in GENES: continue
        t = re.search(r'transcript_id "([^"]+)"', c[8]).group(1); g = m.group(1)
        if c[2] == "exon": ex.setdefault((g, t), []).append((int(c[3]) - 1, int(c[4])))
        else: tx[(g, t)] = (c[0], c[6], int(c[3]) - 1, int(c[4]))
MODEL = {}
for g in GENES:
    cand = [(len(v), max(e for _, e in v) - min(s for s, _ in v), t) for (gg, t), v in ex.items() if gg == g]
    if not cand: continue
    t = sorted(cand)[-1][2]; ch, st, a, b = tx[(g, t)]
    ex_s = sorted(ex[(g, t)]); MODEL[g] = dict(chr=ch, strand=st, tss=a if st == "+" else b - 1, exons=ex_s, tid=t)
def describe(g, pos):
    m = MODEL[g]; ex_s = m["exons"] if m["strand"] == "+" else m["exons"][::-1]
    for i, (s, e) in enumerate(ex_s, 1):
        if s <= pos < e: return f"exon {i}"
    for i in range(len(ex_s) - 1):
        lo, hi = sorted((ex_s[i][1] if m["strand"] == "+" else ex_s[i][0], ex_s[i + 1][0] if m["strand"] == "+" else ex_s[i + 1][1]))
        if lo <= pos <= hi: return f"intron {i + 1}"
    return "outside the transcript"

# ---------------------------------------------------------------- helpers
def flet(fig, ax, l, dx=.06, dy=.02):
    pos = ax.get_position(); fig.text(pos.x0 - dx, pos.y1 + dy, l, fontsize=9, fontweight="bold", va="bottom", ha="left")
def letter(ax, l, x=-.18, y=1.08): ax.text(x, y, l, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom", ha="left")
def binned(bw, ch, a, b, n=400):
    v = bw.stats(ch, int(a), int(b), nBins=n, type="mean"); return np.linspace(a, b, n), np.array([0 if x is None else x for x in v])
def stars(p): return "***" if p < 1e-6 else "**" if p < 1e-3 else "*" if p < .05 else "n.s."
def locus(fig, spec, gene, up, down, samples, label_sites=False, show_input=True, title=None, ylim=None, compact=False):
    m = MODEL[gene]; ch, sg, tss = m["chr"], 1 if m["strand"] == "+" else -1, m["tss"]
    a, b = (tss - up, tss + down) if sg == 1 else (tss - down, tss + up)
    nr = 1 + len(samples) + (1 if show_input else 0) + 1
    hr = [.55] + [1] * len(samples) + ([1] if show_input else []) + [.35 * len(samples) + .3]
    inner = GridSpecFromSubplotSpec(nr, 1, subplot_spec=spec, hspace=.12, height_ratios=hr); axs = []
    # gene model
    ax = fig.add_subplot(inner[0]); axs.append(ax); ax.set_ylim(0, 1); ax.axis("off")
    ax.plot([max(a, m["exons"][0][0]), min(b, m["exons"][-1][1])], [.5, .5], color="k", lw=.6)
    for i, (s, e) in enumerate(m["exons"] if sg == 1 else m["exons"][::-1], 1):
        if e > a and s < b: ax.add_patch(Rectangle((max(s, a), .2), min(e, b) - max(s, a), .6, color="k")); ax.text((max(s, a) + min(e, b)) / 2, .95, str(i), ha="center", va="bottom", fontsize=5)
    if not compact: ax.text(.005 if sg == 1 else .995, .5, f"{gene}", transform=ax.transAxes, ha="left" if sg == 1 else "right", va="center", fontsize=6, style="italic", fontweight="bold")
    ax.set_xlim(a, b)
    # signal tracks
    for i, s in enumerate(samples):
        ax = fig.add_subplot(inner[1 + i]); axs.append(ax); x, v = binned(trk[s], ch, a, b)
        ax.fill_between(x, 0, v, color=SC[s], lw=0, step="mid"); ax.axhline(0, color="k", lw=.3)
        ax.text(.005 if sg == 1 else .995, .93, SL[s], transform=ax.transAxes, ha="left" if sg == 1 else "right", va="top", fontsize=5.2, color=SC[s], fontweight="bold")
        if i == len(samples) // 2: ax.set_ylabel("log$_2$ ChIP/input", fontsize=5.5, labelpad=1)
        if ylim: ax.set_ylim(*ylim)
        ax.set_xlim(a, b); ax.tick_params(labelbottom=False, length=2)
    if show_input:
        ax = fig.add_subplot(inner[1 + len(samples)]); axs.append(ax); x, v = binned(ctrl, ch, a, b)
        ax.fill_between(x, 0, v, color=C["input"], lw=0, step="mid"); ax.text(.005 if sg == 1 else .995, .93, "input", transform=ax.transAxes, ha="left" if sg == 1 else "right", va="top", fontsize=5.2, color="#555", fontweight="bold"); ax.set_ylabel("RPGC", fontsize=5.5, labelpad=1)
        ax.set_xlim(a, b); ax.tick_params(labelbottom=False, length=2)
    # peaks
    ax = fig.add_subplot(inner[-1]); axs.append(ax); ax.set_ylim(0, len(samples)); ax.set_xlim(a, b)
    for j, s in enumerate(samples):
        for q, al in (("q05", .35), ("q01", 1)):
            d = P[(s, q)]; d = d[(d.chr == ch) & (d.s < b) & (d.e > a)]
            for r in d.itertuples(): ax.add_patch(Rectangle((r.s, len(samples) - 1 - j + .15), r.e - r.s, .7, color=SC[s], alpha=al, lw=0))
    ax.set_yticks([])
    for j, s_ in enumerate(samples): ax.text(.005 if sg == 1 else .995, (len(samples) - 1 - j + .5) / len(samples), f"{SL[s_]} peaks", transform=ax.transAxes, ha="left" if sg == 1 else "right", va="center", fontsize=4.6, color=SC[s_], bbox=dict(fc="w", ec="none", pad=.4, alpha=.8))
    ax.spines["left"].set_visible(False)
    kb = np.arange(np.ceil(-up / 1000), np.floor(down / 1000) + 1) if True else []
    ax.set_xticks([tss + sg * k * 1000 for k in kb if a <= tss + sg * k * 1000 <= b]); ax.set_xticklabels([f"{int(k)}" for k in kb if a <= tss + sg * k * 1000 <= b])
    ax.set_xlabel(f"distance from TSS (kb), {ch}", labelpad=1)
    if sg == -1:
        for x_ in axs: x_.set_xlim(b, a)
    for x_ in axs[1:]: x_.spines["left"].set_linewidth(.5)
    return axs, (ch, a, b, sg, tss)

def shade(axs, gene_sites):
    for ax in axs[1:]:
        for w, (lo, hi) in gene_sites.items(): ax.axvspan(lo, hi, color=SITE_COL[w], alpha=.13, lw=0)

# window coordinates for PNPLA3 (from the window table)
wn = W[W.gene == "PNPLA3"].drop_duplicates("window").set_index("window")
sites_xy = {w: (int(wn.loc[w, "start"]), int(wn.loc[w, "end"])) for w in SITE_COL if w in wn.index}

# ================================================================ MAIN FIGURE
fig = plt.figure(figsize=(183 * MM, 150 * MM))
gs = GridSpec(2, 12, figure=fig, height_ratios=[1.5, 1], hspace=.8, wspace=3.1, left=.075, right=.985, top=.885, bottom=.17)
axA, (ch, a, b, sg, tss) = locus(fig, gs[0, :], "PNPLA3", 2500, 7500, S, ylim=(-1.2, 3.5))
shade(axA, sites_xy)
fig.text(.012, .985, ML[0], fontsize=9, fontweight="bold", va="top")
fig.text(.04, .985, "PPARG ChIP-seq in HepG2 (ENCODE ENCSR130VQL, hg38): PNPLA3 locus", fontsize=7, va="top")
fig.text(.985, .985, "peak rows: dark = q<0.01, light = q<0.05 only", fontsize=5, va="top", ha="right", color="#555")
snp_desc = describe("PNPLA3", RS_POS)
if snp_desc == "exon 3" and a <= RS_POS <= b:
    for x_ in axA[1:]: x_.axvline(RS_POS, color="k", lw=.5, ls=":", zorder=5)
    axA[0].text(RS_POS + 130, .5, "rs738409", fontsize=4.6, va="center", ha="left", bbox=dict(fc="w", ec="none", pad=.5))
else: print(f"WARNING: rs738409 position {RS_POS:,} falls in '{snp_desc}' of the exon model, not exon 3 - marker not drawn; check RS_POS / genome build")
for w, (lo, hi) in sites_xy.items():
    axA[0].text((lo + hi) / 2, 1.55, SITE_LAB[w].replace("\n", " "), color=SITE_COL[w], ha="center", va="bottom", fontsize=5.3, transform=axA[0].get_xaxis_transform(), clip_on=False)
# annotation of peaks relative to exons (printed for the legend)
pk = P[("ChIP_pooled", "q01")]; pk = pk[(pk.chr == MODEL["PNPLA3"]["chr"]) & (pk.s < b) & (pk.e > a)]
anno = []
for r in pk.itertuples():
    s_ = r.s + r.summit; dsnp = 0 if r.s <= RS_POS < r.e else min(abs(r.s - RS_POS), abs(r.e - RS_POS))
    anno.append(f"pooled q<0.01 peak chr{r.chr.replace('chr','')}:{r.s:,}-{r.e:,}, summit {s_ - tss:+,} bp from TSS, in {describe('PNPLA3', s_)}; summit {abs(s_ - RS_POS):,} bp and nearest edge {dsnp:,} bp from rs738409 (chr22:{RS_POS:,}, {snp_desc})")

# b: enrichment over input
cats = [("PNPLA3", "upstream_PPRE", "PNPLA3 upstream PPREs"), ("PNPLA3", "intron1_site", "PNPLA3 intron-1 site"), ("PNPLA3", "site_plus4.8kb", "PNPLA3 +4.8 kb"),
        ("PLIN1", "promoter", "PLIN1"), ("PCK1", "promoter", "PCK1"), ("PLIN2", "promoter", "PLIN2"),
        ("FABP4", "promoter", "FABP4"), ("ADIPOQ", "promoter", "ADIPOQ"), ("MYOD1", "promoter", "MYOD1"), ("RHO", "promoter", "RHO"), ("INS", "promoter", "INS")]
cats = [c for c in cats if ((W.gene == c[0]) & (W.window == c[1])).any()]
axB = fig.add_subplot(gs[1, 0:6]); flet(fig, axB, ML[1], dx=.062, dy=.05)
rowsb = []
for i, (g, w, lab) in enumerate(cats):
    d = W[(W.gene == g) & (W.window == w)].set_index("sample")
    pooled = d.loc["ChIP_pooled"]; zero = (pooled.chip_reads + pooled.input_reads) < 5   # no usable coverage in either library
    if zero:   # a pseudocount would draw a spurious positive bar: show the window as uninformative instead
        axB.text(i, .12, "no\ncoverage", ha="center", va="bottom", fontsize=4.6, color="#777", linespacing=1.1)
    else:
        axB.bar(i, pooled.log2_enrichment, .62, color="#d9d9d9", edgecolor="none", zorder=1)
        for s, dx in zip(reps, (-.14, .14)):
            axB.scatter(i + dx, d.loc[s, "log2_enrichment"], s=14, color=SC[s], zorder=3, linewidths=.4, edgecolors="w", marker="o")
        top = max(pooled.log2_enrichment, d.loc[reps[0], "log2_enrichment"], d.loc[reps[1], "log2_enrichment"], 0) + .12
        mh = macs_hit(g, w)
        axB.text(i, top, "no peak" if mh is None else (f"q=10$^{{-{int(round(mh[0]))}}}$" if mh[0] < 100 else "q<10$^{-100}$") + ("" if mh[2] == "q01" else "†"), ha="center", va="bottom", fontsize=4.4, color="#777" if mh is None else "k")
    rowsb.append(dict(gene=g, window=w, informative=not zero, **{f"log2_enrichment_{SL[s]}": (np.nan if zero else d.loc[s, "log2_enrichment"]) for s in S}, pooled_binomial_p_descriptive=pooled.binom_p, macs3_pooled_neg_log10_q=(np.nan if macs_hit(g, w) is None else macs_hit(g, w)[0]), macs3_pooled_fold_change_at_summit=(np.nan if macs_hit(g, w) is None else macs_hit(g, w)[1]), pooled_chip_reads=pooled.chip_reads, input_reads=pooled.input_reads))
axB.axhline(0, color="k", lw=.5); axB.set_xticks(range(len(cats))); axB.set_xticklabels([c[2] for c in cats], fontsize=5.2, rotation=40, ha="right", rotation_mode="anchor")
axB.set_ylabel("log$_2$ (ChIP / input read density)"); axB.set_xlim(-.6, len(cats) - .4)
_v = np.concatenate([W[(W.gene == g) & (W.window == w)].log2_enrichment.values for g, w, _ in cats]); axB.set_ylim(min(-.9, np.nanmin(_v) - .3), np.nanmax(_v) + .95)
for x0, x1, t in ((-.5, 2.5, "PNPLA3 sites"), (2.5, 5.5, "bound genes"), (5.5, 7.5, "adipocyte"), (7.5, len(cats) - .5, "negative")):
    if x1 > x0 and x0 < len(cats):
        axB.annotate(t, xy=((x0 + x1) / 2, 1.0), xycoords=("data", "axes fraction"), ha="center", va="bottom", fontsize=5.3, color="#555")
        if x0 > 0: axB.axvline(x0, color="#ddd", lw=.5, zorder=0)
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
axB.legend(handles=[Line2D([], [], marker="o", ls="", color=C["rep1"], ms=3.5, label=reps[0]), Line2D([], [], marker="o", ls="", color=C["rep2"], ms=3.5, label=reps[1]), Patch(color="#d9d9d9", label="pooled (bar)")],
           loc="upper right", frameon=False, bbox_to_anchor=(1.0, .9), handletextpad=.2, borderaxespad=0)
axB.text(.995, .56, "labels: MACS3 q-value of the pooled\npeak overlapping the window\n(† = q<0.05 peak only)", transform=axB.transAxes, ha="right", va="top", fontsize=4.1, color="#555", linespacing=1.15)
pd.DataFrame(rowsb).to_csv(f"{FD}/source_data_main_panel_b.csv", index=False)

# c: genome-wide percentile (pooled) + input percentile
axC = fig.add_subplot(gs[1, 6:9]); flet(fig, axC, ML[2], dx=.075, dy=.06)
cc = [(g, w, l) for g, w, l in cats if not np.isnan(W[(W.gene == g) & (W.window == w) & (W["sample"] == "ChIP_pooled")].percentile_log2ratio.iloc[0])][:6]
for i, (g, w, lab) in enumerate(cc[::-1]):
    r = W[(W.gene == g) & (W.window == w) & (W["sample"] == "ChIP_pooled")].iloc[0]
    axC.plot(sorted([r.input_percentile, r.percentile_log2ratio]), [i, i], color="#bbb", lw=.8, zorder=1)
    axC.scatter(r.percentile_log2ratio, i, s=18, color=C["pooled"], zorder=3, label="ChIP/input" if i == 0 else None)
    axC.scatter(r.input_percentile, i, s=18, marker="D", facecolor="w", edgecolor=C["input"], linewidth=.9, zorder=3, label="input only" if i == 0 else None)
axC.axvline(99, color="r", lw=.6, ls="--"); axC.set_yticks(range(len(cc))); axC.set_yticklabels([{"PNPLA3 intron-1 site": "PNPLA3\nintron-1", "PNPLA3 +4.8 kb": "PNPLA3\n+4.8 kb", "PNPLA3 upstream PPREs": "PNPLA3\nupstream"}.get(c[2], c[2]) for c in cc[::-1]], fontsize=4.8)
axC.set_xlim(0, 101); axC.set_xticks([0, 50, 100]); axC.set_ylim(-.9, len(cc) - .4); axC.set_xlabel("genome-wide percentile"); axC.legend(frameon=False, loc="lower left", bbox_to_anchor=(-.02, 1.0), ncol=2, handletextpad=.1, columnspacing=.8)
axC.text(98.6, -.62, "99th", color="r", fontsize=5, ha="right", va="center")

# d: response to TZD for HepG2-bound vs unbound genes
def response(de, sym, lfc, padj, bset):
    de = de.dropna(subset=[padj]).copy(); de["bound"] = de[sym].isin(bset); de["ind"] = (de[padj] < .05) & (de[lfc] > 0)
    a_, b_ = int((de.bound & de.ind).sum()), int((de.bound & ~de.ind).sum()); c_, d_ = int((~de.bound & de.ind).sum()), int((~de.bound & ~de.ind).sum())
    o, p = stats.fisher_exact([[a_, b_], [c_, d_]], alternative="greater")
    return dict(pct_bound=100 * a_ / max(a_ + b_, 1), pct_unbound=100 * c_ / max(c_ + d_, 1), n_bound=a_ + b_, n_unbound=c_ + d_, OR=o, p=p), de
primary = next((d for d in ("pooled_supported_by_both_reps_q05", "pooled_q01", "consensus_q01_both_replicates") if os.path.exists(f"{FA}/bound_genes_{d}.csv")), None)
RESP, DEs = {}, {}
def load_resp():
    out = {}
    if DE_CSV and os.path.exists(DE_CSV):
        de = pd.read_csv(DE_CSV); de["sym"] = de.GeneName.str.upper(); de = de.sort_values("baseMean", ascending=False).drop_duplicates("sym"); out["Mouse liver, TZD\n(GSE162276)"] = (de, "sym", "log2FoldChange", "padj", "stat")
    if HASC_CSV and os.path.exists(HASC_CSV):
        h = pd.read_csv(HASC_CSV).drop_duplicates("Gene"); h["sym"] = h.Gene.str.upper(); out["hASC adipocytes, rosiglitazone\n(GSE115421)"] = (h, "sym", "log2FC_Rosi_vs_DMSO", "padj", "t_stat")
    return out
DATA = load_resp()
axD = fig.add_subplot(gs[1, 9:12]); flet(fig, axD, ML[3], dx=.075, dy=.05)
if PI is not None and len(PI):
    rows = []; topmax = float(max(PI.pct_bound_induced.max(), PI.pct_unbound_induced.max()))
    for i, r in enumerate(PI.itertuples()):
        axD.bar(i - .2, r.pct_bound_induced, .36, color=C["bound"], label="HepG2 PPARG-bound" if i == 0 else None); axD.bar(i + .2, r.pct_unbound_induced, .36, color=C["unbound"], label="unbound" if i == 0 else None)
        nu = r.unbound_induced + r.unbound_not_induced
        axD.text(i, topmax * 1.22, f"OR {r.OR:.2f} ({r.OR_ci_lo:.1f}–{r.OR_ci_hi:.1f})\nP = {r.fisher_p_two_sided:.2g}", ha="center", va="bottom", fontsize=4.6, linespacing=1.15)
        rows.append(r._asdict())
    if PX is not None and len(PX):
        x_ = PX.iloc[0]; pp = x_.get("perm_p_difference_expression_stratified", np.nan); pp = x_["perm_p_difference"] if pd.isna(pp) else pp
        axD.plot([0, 1], [topmax * 1.92] * 2, color="k", lw=.5)
        axD.text(.5, topmax * 1.96, f"OR ratio hASC/liver {x_['ratio_of_ORs_hASC_over_liver']:.1f}\nP(difference) = {pp:.2g}", ha="center", va="bottom", fontsize=4.2, linespacing=1.15)
    pd.DataFrame(rows).drop(columns=["Index"], errors="ignore").to_csv(f"{FD}/source_data_main_panel_d.csv", index=False)
    axD.set_xticks(range(len(PI))); axD.set_xticklabels([("Mouse liver,\nTZD\n(GSE162276)" if str(r.dataset).startswith("liver") else "hASC adipocytes,\nrosiglitazone\n(GSE115421)") + f"\nbound {r.bound_induced:,}/{r.n_bound:,}\nunbound {r.unbound_induced:,}/{r.unbound_induced + r.unbound_not_induced:,}" for r in PI.itertuples()], fontsize=4.6)
    axD.set_ylabel("genes induced (%, padj < 0.05)"); axD.set_ylim(0, topmax * 3.05); axD.legend(frameon=False, loc="upper left", bbox_to_anchor=(0, 1.04), fontsize=5)
elif DATA and primary:
    bs = {x.upper() for x in pd.read_csv(f"{FA}/bound_genes_{primary}.csv").gene}; rows = []
    for i, (lab, (de, sym, lfc, padj, st)) in enumerate(DATA.items()):
        r, de2 = response(de, sym, lfc, padj, bs); RESP[lab] = r; DEs[lab] = de2; rows.append(dict(dataset=lab.replace("\n", " "), bound_definition=primary, **r))
        axD.bar(i - .2, r["pct_bound"], .36, color=C["bound"], label="HepG2 PPARG-bound" if i == 0 else None); axD.bar(i + .2, r["pct_unbound"], .36, color=C["unbound"], label="unbound" if i == 0 else None)
        top = max(r["pct_bound"], r["pct_unbound"]); axD.text(i, top * 1.08 + .1, f"OR {r['OR']:.2f}\nP = {r['p']:.2g}", ha="center", va="bottom", fontsize=5, linespacing=1.1)
        axD.text(i - .2, .1, f"n={r['n_bound']:,}", ha="center", va="bottom", fontsize=4.6, color="w", rotation=90); axD.text(i + .2, .1, f"n={r['n_unbound']:,}", ha="center", va="bottom", fontsize=4.6, color="#333", rotation=90)
    pd.DataFrame(rows).to_csv(f"{FD}/source_data_main_panel_d.csv", index=False)
    axD.set_xticks(range(len(DATA))); axD.set_xticklabels([l.replace(", ", ",\n") for l in DATA], fontsize=5); axD.set_ylabel("genes induced (%, padj < 0.05)")
    axD.set_ylim(0, max(max(v["pct_bound"], v["pct_unbound"]) for v in RESP.values()) * 2.1); axD.legend(frameon=False, loc="upper left", bbox_to_anchor=(0, 1.05))
else: axD.text(.5, .5, "response tables not found", transform=axD.transAxes, ha="center")
fig.savefig(f"{FD}/HepG2_main_figure.pdf"); fig.savefig(f"{FD}/HepG2_main_figure.png", dpi=300); plt.close(fig)

# ================================================================ SUPPLEMENTARY FIGURE
fig = plt.figure(figsize=(183 * MM, 250 * MM))
gs = GridSpec(4, 3, figure=fig, hspace=.85, wspace=.6, left=.1, right=.985, top=.965, bottom=.05, height_ratios=[1, 1.5, 1.9, 1.1])
# a: peak numbers
ax = fig.add_subplot(gs[0, 0]); flet(fig, ax, "a", dx=.075, dy=.02)
t = QC.drop_duplicates(["sample", "q"]); x = np.arange(len(S)); w_ = .36
for j, (q, al) in enumerate((("q01", 1), ("q05", .4))):
    v = [int(t[(t["sample"] == s) & (t.q == q)].n_peaks.iloc[0]) for s in S]; ax.bar(x + (j - .5) * w_, v, w_, color=[SC[s] for s in S], alpha=al, label=f"q < {0.01 if q=='q01' else 0.05}")
    for xi, vi in zip(x + (j - .5) * w_, v): ax.text(xi, vi, f"{vi/1000:.1f}k", ha="center", va="bottom", fontsize=4.8)
ax.set_xticks(x); ax.set_xticklabels([SL[s].replace("ENCFF", "ENCFF\n") for s in S], fontsize=5.2); ax.set_ylabel("peaks (MACS3)"); ax.set_title("Peak numbers", fontsize=6.5, loc="left")
# b: genomic annotation
ax = fig.add_subplot(gs[0, 1]); flet(fig, ax, "b", dx=.075, dy=.02)
cats_ = ["promoter (TSS+/-3 kb)", "exon", "intron", "intergenic"]; qq = QC[QC.q == "q01"].pivot(index="sample", columns="feature", values="frac")[cats_].loc[S]
bot = np.zeros(len(S))
for cname, col in zip(cats_, ["#7b3f98", "#e69f00", "#56b4e9", "#b0b0b0"]):
    ax.bar(range(len(S)), qq[cname].values, bottom=bot, color=col, label=cname.replace(" (TSS+/-3 kb)", " (±3 kb)"), width=.65); bot += qq[cname].values
ax.set_xticks(range(len(S))); ax.set_xticklabels([SL[s].replace("ENCFF", "ENCFF\n") for s in S], fontsize=5.2); ax.set_ylabel("fraction of peak summits"); ax.set_title("Genomic annotation (q < 0.01)", fontsize=6.5, loc="left")
ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(.5, -.2), ncol=2, fontsize=5, columnspacing=.8, handlelength=1)
# c: replicate concordance
ax = fig.add_subplot(gs[0, 2]); flet(fig, ax, "c", dx=.12, dy=.02)
lab_c = [f"{r.split(' peaks overlapping ')[0].replace('ENCFF','')}\n→ {r.split(' peaks overlapping ')[1].replace('ENCFF','')}" for r in CONC["set"]]
ax.barh(range(len(CONC))[::-1], 100 * CONC.fraction, color="#7b3f98", height=.6); ax.set_yticks(range(len(CONC))[::-1]); ax.set_yticklabels(lab_c, fontsize=4.6)
ax.set_xlim(0, 118); ax.set_xticks([0, 50, 100]); ax.set_xlabel("% of peaks with an overlapping peak"); ax.set_title("Replicate concordance", fontsize=6.5, loc="left")
for yi, v in zip(range(len(CONC))[::-1], 100 * CONC.fraction): ax.text(v + 1, yi, f"{v:.0f}", va="center", fontsize=5)

# d-f: control loci
sub = gs[1, :].subgridspec(1, 3, wspace=.75)
for k, (g, lab, up, dn) in enumerate((("PLIN1", "d", 4000, 4000), ("PCK1", "e", 4000, 4000), ("FABP4", "f", 4000, 4000))):
    if g not in MODEL: continue
    axs, _ = locus(fig, sub[0, k], g, up, dn, ["ChIP_pooled"], compact=True, ylim=(-1.2, 3.5))
    flet(fig, axs[0], lab, dx=.075, dy=.03)
    axs[0].set_title(f"{g}: {'bound control' if g != 'FABP4' else 'adipocyte gene, unbound'}", fontsize=6, loc="left", pad=9, style="italic")

# g: PNPLA3 window table
ax = fig.add_subplot(gs[2, 0:2]); flet(fig, ax, "g", dx=.075, dy=.03); ax.axis("off")
pn = W[W.gene == "PNPLA3"]; hdr = ["site", "chr22 window (hg38)", "sample", "ChIP\nreads", "input\nreads", "norm.\nChIP/input", "MACS3 FC\n(summit)", "MACS3\n−log10 q", "peak\nq<.01", "peak\nq<.05"]
GN = {"upstream_PPRE": "upstream PPRE", "intron1_site": "intron-1 site", "site_plus4.8kb": "+4.8 kb site", "promoter": "promoter ±2 kb"}
def peak_stat(sample, r):
    for lvl in ("q01", "q05"):
        d = P[(sample, lvl)]; d = d[(d.chr == r.chr) & (d.s < r.end) & (d.e > r.start)]
        if len(d): k = d.nlq.astype(float).idxmax(); return f"{float(d.loc[k, 'fc']):.1f}", f"{float(d.loc[k, 'nlq']):.1f}" + ("" if lvl == "q01" else "†")
    return "–", "–"
cells = []
for w in ("upstream_PPRE", "intron1_site", "site_plus4.8kb", "promoter"):
    for s in S:
        r = pn[(pn.window == w) & (pn["sample"] == s)]
        if r.empty: continue
        r = r.iloc[0]; f0 = s == S[0]; cov = (r.chip_reads + r.input_reads) >= 5; fc_, q_ = peak_stat(s, r)
        cells.append([GN[w] if f0 else "", f"{int(r.start):,}-{int(r.end):,}" if f0 else "", SL[s].replace("ENCFF", ""), f"{int(r.chip_reads)}", f"{int(r.input_reads)}",
                      f"{r.fold_enrichment:.1f}" if cov else "n/a", fc_, q_, "yes" if r.peak_q01 is True or r.peak_q01 == True else "no", "yes" if r.peak_q05 is True or r.peak_q05 == True else "no"])
cw = np.array([.12, .215, .07, .055, .055, .09, .08, .08, .06, .06]); cw = list(cw / cw.sum())
tb = ax.table(cellText=cells, colLabels=hdr, bbox=[0, 0, .94, .96], cellLoc="center", colWidths=cw); tb.auto_set_font_size(False); tb.set_fontsize(4.1)
for (i, j), c_ in tb.get_celld().items():
    c_.set_linewidth(.2); c_.set_edgecolor("#ccc")
    if i == 0: c_.set_facecolor("#eee"); c_.set_text_props(fontweight="bold")
    if j in (0, 1) and i > 0: c_.set_text_props(ha="left"); c_._loc = "left"
ax.text(0, -.04, "n/a: < 5 reads in ChIP + input combined (no usable coverage). Windows overlap (promoter contains the start of the intron-1 window).\nDepth-norm. ChIP/input: (ChIP reads / library) / (input reads / library), 0.5 pseudocount - descriptive only. MACS3 columns: fold change at summit and −log10 q of the overlapping\npeak (q<0.01; † = q<0.05 peak only), from the MACS3 output. Per-sample peak statistics and a local-background re-test are in revalidation/peak_level_statistics.csv and local_background_test.csv.", transform=ax.transAxes, fontsize=3.9, va="top", color="#555")
ax.set_title("PNPLA3 windows: read counts against input", fontsize=6.5, loc="left")
# h: motif
ax = fig.add_subplot(gs[2, 2]); flet(fig, ax, "h", dx=.12, dy=.03)
if MOT is not None and len(MOT):
    m_ = MOT.copy(); m_["lab"] = m_.site.str.replace("_", " ").str.replace(r"\(.*\)", "", regex=True).str.strip(); yy = np.arange(len(m_))[::-1]
    ax.barh(yy, m_.percentile_vs_random, color=["#7b3f98" if "PLIN1" in s or "PCK1" in s else "#2ca02c" if "intron1" in s else "#ff9f1c" if "4.8" in s else "#d62728" for s in m_.site], height=.6)
    ax.axvline(95, color="k", ls="--", lw=.6); ax.set_yticks(yy); ax.set_yticklabels(m_.lab, fontsize=4.6); ax.set_xlim(0, 100); ax.set_xlabel("percentile of motif score\nvs GC-matched random windows")
    ax.text(95, len(m_) - .35, "95th", fontsize=5, ha="right"); ax.set_title("PPARG::RXRA motif score", fontsize=6.5, loc="left", pad=8)

# i: Hallmark
ax = fig.add_subplot(gs[3, 0]); flet(fig, ax, "i", dx=.075, dy=.02)
if HALL is not None and len(HALL):
    h_ = HALL.sort_values("pval").head(8)[::-1]; import textwrap
    ax.barh([textwrap.fill(t.replace("HALLMARK_", "").replace("_", " ").capitalize(), 15) for t in h_.Term], -np.log10(h_.padj.clip(lower=1e-300)), color="#7b3f98", height=.65)
    ax.axvline(-np.log10(.05), color="r", lw=.6, ls="--"); ax.set_xlabel("−log$_{10}$ FDR"); ax.tick_params(axis="y", labelsize=4.6); ax.set_title("Hallmark, HepG2-bound genes", fontsize=6.5, loc="left", pad=8)
# j: overlap with hASC
ax = fig.add_subplot(gs[3, 1]); flet(fig, ax, "j", dx=.075, dy=.02)
if OVL is not None and len(OVL):
    names = {"consensus_q01_both_replicates": "both replicates\nq<0.01", "pooled_q01": "pooled\nq<0.01", "pooled_supported_by_both_reps_q05": "pooled, both\nreps q<0.05"}
    x = np.arange(len(OVL))
    ax.bar(x - .2, OVL.expected, .36, color=C["unbound"], label="expected"); ax.bar(x + .2, OVL.overlap, .36, color=C["bound"], label="observed")
    for xi, r in zip(x, OVL.itertuples()): ax.text(xi + .2, r.overlap, f"P={r.hypergeom_p:.0e}", ha="center", va="bottom", fontsize=4.6, rotation=90)
    ax.set_xticks(x); ax.set_xticklabels([names.get(d, d) for d in OVL.bound_definition], fontsize=4.8); ax.set_ylabel("genes bound in HepG2 and hASC"); ax.legend(frameon=False, loc="upper left")
    ax.set_ylim(0, max(OVL.overlap.max(), OVL.expected.max()) * 1.7); ax.set_title("Overlap with hASC PPARG-bound genes", fontsize=6.5, loc="left", pad=8)
# k: ECDFs
ax = fig.add_subplot(gs[3, 2]); flet(fig, ax, "k", dx=.075, dy=.02)
if UNI is not None and len(UNI):
    for col_, ls in (("stat_liver", "-"), ("stat_hasc", "--")):
        for flag, colr in ((False, C["unbound"]), (True, C["bound"])):
            v = np.sort(UNI.loc[UNI.bound_A == flag, col_].dropna().values); ax.plot(v, np.arange(1, len(v) + 1) / len(v), color=colr, ls=ls, lw=.9)
    ax.plot([], [], color=C["bound"], label="bound"); ax.plot([], [], color=C["unbound"], label="unbound"); ax.plot([], [], color="k", ls="-", lw=.9, label="liver TZD"); ax.plot([], [], color="k", ls="--", lw=.9, label="hASC rosi")
    ax.set_xlim(-6, 6); ax.set_xlabel("response statistic (t or Wald)"); ax.set_ylabel("cumulative fraction"); ax.legend(frameon=False, ncol=1, loc="upper left", fontsize=5); ax.set_title("Response, bound vs unbound", fontsize=6.5, loc="left", pad=8)
elif DEs:
    for (lab, de2), ls in zip(DEs.items(), ("-", "--")):
        st = DATA[lab][4]; de2 = de2.dropna(subset=[st])
        for flag, col in ((False, C["unbound"]), (True, C["bound"])):
            v = np.sort(de2.loc[de2.bound == flag, st].values); ax.plot(v, np.arange(1, len(v) + 1) / len(v), color=col, ls=ls, lw=.9)
    ax.plot([], [], color=C["bound"], label="bound"); ax.plot([], [], color=C["unbound"], label="unbound"); ax.plot([], [], color="k", ls="-", lw=.9, label="liver TZD"); ax.plot([], [], color="k", ls="--", lw=.9, label="hASC rosi")
    ax.set_xlim(-6, 6); ax.set_xlabel("response statistic (t or Wald)"); ax.set_ylabel("cumulative fraction"); ax.legend(frameon=False, ncol=1, loc="upper left", fontsize=5); ax.set_title("Response, bound vs unbound", fontsize=6.5, loc="left", pad=8)
fig.savefig(f"{FD}/HepG2_supplementary_figure.pdf"); fig.savefig(f"{FD}/HepG2_supplementary_figure.png", dpi=300); plt.close(fig)

# ---------------------------------------------------------------- legends helper (printed)
print("\nPNPLA3 peak positions (pooled q<0.01, within the plotted window):")
for s_ in anno: print("  " + s_)
if RESP:
    for k, v in RESP.items(): print(f"  {k.replace(chr(10), ' ')}: bound {v['pct_bound']:.1f}% (n={v['n_bound']}) vs unbound {v['pct_unbound']:.1f}% (n={v['n_unbound']}); OR {v['OR']:.2f}, one-sided Fisher P={v['p']:.2g}")
extra = []
if PI is not None:
    for r in PI.itertuples(): extra.append(f"{r.dataset}: bound {r.bound_induced}/{r.n_bound} induced ({r.pct_bound_induced:.1f}%) vs unbound {r.unbound_induced}/{r.unbound_induced + r.unbound_not_induced} ({r.pct_unbound_induced:.1f}%); OR {r.OR:.2f} (95% CI {r.OR_ci_lo:.2f}-{r.OR_ci_hi:.2f}), two-sided Fisher P={r.fisher_p_two_sided:.3g}; expression-adjusted OR {r.adj_OR_expression:.2f} ({r.adj_OR_ci_lo:.2f}-{r.adj_OR_ci_hi:.2f}), P={r.adj_p_expression:.3g}")
if PX is not None:
    for r in PX.itertuples(): extra.append(f"between-dataset difference: OR ratio hASC/liver {r.ratio_of_ORs_hASC_over_liver:.2f}; Woolf P={r.woolf_p:.3g}; GEE interaction P={getattr(r, 'gee_interaction_p', float('nan')):.3g}; permutation P={r.perm_p_difference:.3g}; expression-stratified permutation P={r.perm_p_difference_expression_stratified:.3g}")
for e_ in extra: print("  " + e_)
open(f"{FD}/figure_legend_numbers.txt", "w").write("\n".join(anno + extra) + "\n")
print(f"\nDone. Figures in {FD}/")

#!/usr/bin/env python3
"""HepG2 PPARG ChIP-seq (ENCODE ENCSR130VQL): the same analysis as GSE115421, in one figure.

Panels (analogous to Fig. 6e-l for GSE115421):
  e  PNPLA3 locus: per-replicate RPGC, control RPGC, per-replicate log2(ChIP/input)
  f  FABP4 locus: same
  g  mean signal at TSS +/-2 kb (RPGC and log2 ChIP/input), positive vs negative control genes
  h  binding reproducibility: replicates with a peak within TSS +/-2 kb (2 replicates -> 0/1/2)
  i  response to TZD in liver (GSE162276 mouse, human symbol = upper-case mouse symbol), bound vs unbound
  j  % induced among bound vs unbound genes (Fisher) + expression-adjusted logistic regression
  k  raw signal heatmap, genes x replicates
  l  PNPLA3 windows (promoter, upstream PPRE, intron-1 site): log2 ChIP/input

Environment variables:
  OUT   results/hepg2_ppar   (same OUT as step 1)
  GTF   path to GENCODE GTF(.gz) (hg38)
  CTRL_BAM  results/hepg2_ppar/bam/control_sub20.bam   (for the control RPGC track)
  DE_CSV    GSE162276_TZDvsControl_DE_results.csv (columns: GeneName, baseMean, log2FoldChange, stat, padj)
Needs: pyBigWig pandas numpy scipy matplotlib statsmodels(optional); deepTools bamCoverage on PATH
(only if the control RPGC track does not exist yet).
"""
import os, re, gzip, glob, subprocess, sys
import numpy as np, pandas as pd, pyBigWig
from scipy.stats import fisher_exact, mannwhitneyu
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

OUT = os.environ.get("OUT", "results/hepg2_ppar")
GTF = os.environ["GTF"]
CTRL_BAM = os.environ.get("CTRL_BAM", f"{OUT}/bam/control_sub20.bam")
DE_CSV = os.environ.get("DE_CSV", "")
if not DE_CSV:      # look for the file in the usual places
    for d in ("ref", "data", ".", "..", "../.."):
        c = os.path.join(d, "GSE162276_TZDvsControl_DE_results.csv")
        if os.path.exists(c): DE_CSV = c; print("found DE table:", c); break
if not DE_CSV or not os.path.exists(DE_CSV):
    sys.exit(f"STOP: the TZD DE table was not found (DE_CSV='{DE_CSV}'). Upload GSE162276_TZDvsControl_DE_results.csv to ref/ and rerun, "
             "or set DE_CSV=<real path>. (Running without it silently drops panels i and j.)")
PROMOTER = 3000          # bound = consensus peak within TSS +/- 3 kb (same as GSE115421 analysis)
POS = ["FABP4", "ADIPOQ", "CD36", "PLIN1", "PCK1", "LPL", "ANGPTL4", "PLIN2", "PNPLA2"]
NEG = ["MYOD1", "RHO", "INS"]
SHOW = ["PNPLA3"] + POS + NEG
WIN = {"promoter": (-2000, 2000), "upstream_PPRE": (-850, -650), "intron1_site": (1000, 2300)}
GSIZE = 2913022398

# ---------------------------------------------------------------- inputs
reps = sorted(os.path.basename(p).replace(".RPGC.bw", "") for p in glob.glob(f"{OUT}/tracks/*.RPGC.bw")
              if not os.path.basename(p).startswith(("control", "ChIP_pooled")) and os.path.exists(f"{OUT}/peaks/" + os.path.basename(p).replace(".RPGC.bw", "") + "_peaks.narrowPeak"))
if len(reps) < 2: sys.exit(f"Need >=2 replicate RPGC tracks in {OUT}/tracks (found {reps})")
for r in reps:
    for f in (f"{OUT}/tracks/{r}.log2ChIPoverInput.bw", f"{OUT}/peaks/{r}_peaks.narrowPeak"):
        if not os.path.exists(f): sys.exit(f"Missing {f} -- step 1 did not finish for {r}")
ctrl_bw = f"{OUT}/tracks/control.RPGC.bw"
if not os.path.exists(ctrl_bw):
    subprocess.check_call(["bamCoverage", "-b", CTRL_BAM, "-o", ctrl_bw, "--normalizeUsing", "RPGC",
                           "--effectiveGenomeSize", str(GSIZE), "--binSize", "10", "--extendReads", "200", "-p", "8"])
bw = {r: pyBigWig.open(f"{OUT}/tracks/{r}.RPGC.bw") for r in reps}
lr = {r: pyBigWig.open(f"{OUT}/tracks/{r}.log2ChIPoverInput.bw") for r in reps}
cb = pyBigWig.open(ctrl_bw)
peaks = {r: pd.read_csv(f"{OUT}/peaks/{r}_peaks.narrowPeak", sep="\t", header=None, usecols=[0, 1, 2, 8],
                        names=["chr", "s", "e", "nlq"]) for r in reps}

# ---------------------------------------------------------------- annotation (TSS per gene, all transcripts)
op = gzip.open if GTF.endswith(".gz") else open
tss_all, gene_tss = [], {}
with op(GTF, "rt") as f:
    for ln in f:
        if ln.startswith("#"): continue
        c = ln.rstrip("\n").split("\t")
        if len(c) < 9 or c[2] != "transcript": continue
        m = re.search(r'gene_name "([^"]+)"', c[8])
        if not m: continue
        t = int(c[3]) - 1 if c[6] == "+" else int(c[4]) - 1
        tss_all.append((m.group(1), c[0], t, c[6]))
tss_df = pd.DataFrame(tss_all, columns=["gene", "chr", "tss", "strand"]).drop_duplicates()
tss_df["key"] = np.where(tss_df.strand == "+", tss_df.tss, -tss_df.tss)
main = tss_df.sort_values("key").drop_duplicates("gene").set_index("gene")   # most upstream TSS per gene
print("PNPLA3 TSS (expected ~chr22:43,923,805):", main.loc["PNPLA3", ["chr", "tss", "strand"]].tolist())

def merged(df):
    out = {}
    for ch, d in df.groupby("chr"):
        d = d.sort_values("s"); s, e = d.s.values, d.e.values
        ms, me = [s[0]], [e[0]]
        for a, b in zip(s[1:], e[1:]):
            if a <= me[-1]: me[-1] = max(me[-1], b)
            else: ms.append(a); me.append(b)
        out[ch] = (np.array(ms), np.array(me))
    return out

def hit(iv, ch, a, b):
    if ch not in iv: return False
    s, e = iv[ch]; i = np.searchsorted(s, b, "right") - 1
    return i >= 0 and e[i] >= a

# consensus = rep1 peaks overlapping any rep2 peak (both replicates)
r1, r2 = reps[:2]; iv2 = merged(peaks[r2])
cons = peaks[r1][[hit(iv2, c, s, e) for c, s, e in zip(peaks[r1].chr, peaks[r1].s, peaks[r1].e)]]
iv_c = merged(cons); iv_rep = {r: merged(peaks[r]) for r in reps}
bound = {g for g, row in main.iterrows() if hit(iv_c, row.chr, row.tss - PROMOTER, row.tss + PROMOTER)}
pd.DataFrame({"gene": sorted(bound)}).to_csv(f"{OUT}/hepg2_bound_genes_consensus.csv", index=False)
print(f"peaks: {', '.join(f'{r}={len(peaks[r])}' for r in reps)}; consensus={len(cons)}; bound genes (TSS+/-{PROMOTER}) = {len(bound)}")

def wmean(b, ch, a, e):
    v = b.stats(ch, int(max(a, 0)), int(e), type="mean")[0]; return 0.0 if v is None else v

rows = []
for g in SHOW:
    if g not in main.index: print("not in GTF:", g); continue
    ch, t = main.loc[g, "chr"], int(main.loc[g, "tss"]); sg = 1 if main.loc[g, "strand"] == "+" else -1
    for w, (a, b) in WIN.items():
        if g != "PNPLA3" and w != "promoter": continue
        lo, hi = sorted((t + sg * a, t + sg * b))
        for r in reps:
            rows.append(dict(gene=g, window=w, rep=r, RPGC=wmean(bw[r], ch, lo, hi), log2_ChIP_over_input=wmean(lr[r], ch, lo, hi),
                             control_RPGC=wmean(cb, ch, lo, hi), peak_in_rep=hit(iv_rep[r], ch, lo, hi),
                             consensus_peak=hit(iv_c, ch, lo, hi)))
met = pd.DataFrame(rows); met.to_csv(f"{OUT}/hepg2_gene_metrics.csv", index=False)

# ---------------------------------------------------------------- RNA-seq response (liver TZD) integration
stats_txt, rna = [], None
if DE_CSV and os.path.exists(DE_CSV):
    de = pd.read_csv(DE_CSV).dropna(subset=["padj", "stat"]).copy()
    de["sym"] = de.GeneName.str.upper()
    de = de.sort_values("baseMean", ascending=False).drop_duplicates("sym")
    de["bound"] = de.sym.isin(bound); de["induced"] = (de.padj < 0.05) & (de.log2FoldChange > 0)
    a_ = int((de.bound & de.induced).sum()); b_ = int((de.bound & ~de.induced).sum())
    c_ = int((~de.bound & de.induced).sum()); d_ = int((~de.bound & ~de.induced).sum())
    orr, pf = fisher_exact([[a_, b_], [c_, d_]], alternative="greater")
    mw = mannwhitneyu(de[de.bound].stat, de[~de.bound].stat, alternative="greater").pvalue
    stats_txt += [f"Liver TZD (GSE162276, mouse) vs HepG2 PPARG-bound genes (human symbol = upper-case mouse symbol)",
                  f"bound n={int(de.bound.sum())}, unbound n={int((~de.bound).sum())}",
                  f"induced: bound {a_}/{a_+b_} = {100*a_/max(a_+b_,1):.1f}%, unbound {c_}/{c_+d_} = {100*c_/max(c_+d_,1):.1f}%",
                  f"Fisher OR={orr:.2f}, one-sided P={pf:.3g}; Mann-Whitney (t-stat) P={mw:.3g}"]
    try:
        import statsmodels.api as sm
        X = sm.add_constant(pd.DataFrame({"bound": de.bound.astype(float), "log_baseMean": np.log10(de.baseMean + 1)}))
        fit = sm.Logit(de.induced.astype(float), X).fit(disp=0)
        stats_txt.append(f"Expression-adjusted logistic regression: bound OR={np.exp(fit.params['bound']):.2f} "
                         f"(95% CI {np.exp(fit.conf_int().loc['bound',0]):.2f}-{np.exp(fit.conf_int().loc['bound',1]):.2f}), P={fit.pvalues['bound']:.3g}")
    except Exception as ex:
        stats_txt.append(f"(adjusted regression skipped: {ex})")
    if "PNPLA3" in set(de.sym):
        p = de[de.sym == "PNPLA3"].iloc[0]
        stats_txt.append(f"PNPLA3 in liver TZD: log2FC={p.log2FoldChange:.2f}, t={p.stat:.2f}, padj={p.padj:.3g}, bound={bool(p.bound)}")
    rna = de
else:
    stats_txt.append("DE_CSV not found: panels i,j skipped")

# ---------------------------------------------------------------- QC verdict (can the null/positive be interpreted?)
pc = met[(met.window == "promoter") & met.gene.isin(POS)].groupby("gene").consensus_peak.first()
nc = met[(met.window == "promoter") & met.gene.isin(NEG)].groupby("gene").consensus_peak.first()
pn = met[(met.gene == "PNPLA3")]
qc = [f"Positive controls with a consensus peak at TSS+/-{PROMOTER//1000}kb: {int(pc.sum())}/{len(pc)} ({', '.join(pc[pc].index)})",
      f"Negative controls with a consensus peak: {int(nc.sum())}/{len(nc)}",
      "PNPLA3 consensus peak by window: " + "; ".join(f"{w}={bool(pn[pn.window==w].consensus_peak.iloc[0])}" for w in WIN),
      "PNPLA3 mean log2(ChIP/input) by window (rep1,rep2): " + "; ".join(
          f"{w}=" + ",".join(f'{v:.2f}' for v in pn[pn.window==w].log2_ChIP_over_input) for w in WIN)]
open(f"{OUT}/hepg2_integration_and_QC.txt", "w").write("\n".join(stats_txt + [""] + qc) + "\n"); print("\n".join(stats_txt + qc))

# ---------------------------------------------------------------- figure
fig = plt.figure(figsize=(15, 15)); gs = GridSpec(3, 6, figure=fig, height_ratios=[1.15, 1, 1], hspace=.45, wspace=.8)
cols = ["#4C6A9C", "#C0504D"]

def locus(subspec, gene, up, down, letter):
    ch, t = main.loc[gene, "chr"], int(main.loc[gene, "tss"]); sg = 1 if main.loc[gene, "strand"] == "+" else -1
    lo, hi = t - (up if sg > 0 else down), t + (down if sg > 0 else up)
    xs = np.arange(lo, hi, 50); inner = subspec.subgridspec(2 * len(reps) + 1, 1, hspace=.1); axs = []
    items = [(f"{r[-7:]} RPGC", bw[r], cols[i]) for i, r in enumerate(reps)] + [("control RPGC", cb, "gray")] + \
            [(f"{r[-7:]} log2 ChIP/input", lr[r], cols[i]) for i, r in enumerate(reps)]
    for k, (lab, b, col) in enumerate(items):
        ax = fig.add_subplot(inner[k]); ys = [wmean(b, ch, x, x + 50) for x in xs]
        ax.fill_between(xs, ys, step="mid", color=col); ax.axvline(t, color="k", lw=.8); ax.set_ylabel(lab, fontsize=6, rotation=0, ha="right", va="center")
        ax.tick_params(labelsize=6); ax.set_xlim(lo, hi); axs.append(ax)
        if gene == "PNPLA3":
            for (a, b2), c in zip([WIN["upstream_PPRE"], WIN["intron1_site"]], ["tab:red", "tab:green"]):
                ax.axvspan(t + a, t + b2, color=c, alpha=.15)
        if k < len(items) - 1: ax.set_xticklabels([])
    axs[0].set_title(f"{gene} {ch} (hg38; black = TSS)", fontsize=8, loc="left")
    fig.text(subspec.get_position(fig).x0 - .02, subspec.get_position(fig).y1 + .01, letter, fontsize=16, weight="bold")

locus(gs[0, 0:3], "PNPLA3", 6000, 12000, "e"); locus(gs[0, 3:6], "FABP4", 6000, 8000, "f")

genes_show = [g for g in SHOW if g in main.index]
ax = fig.add_subplot(gs[1, 0:2]); w0 = met[met.window == "promoter"]
for i, r in enumerate(reps):
    d = w0[w0.rep == r].set_index("gene").reindex(genes_show); ax.scatter(range(len(genes_show)), d.RPGC, color=cols[i], s=22, label=r)
ax.set_xticks(range(len(genes_show))); ax.set_xticklabels(genes_show, rotation=75, fontsize=7); ax.set_ylabel("mean RPGC TSS$\\pm$2 kb"); ax.legend(fontsize=6)
ax.set_title("g  signal at TSS$\\pm$2 kb", loc="left", fontsize=9)

ax = fig.add_subplot(gs[1, 2:4]); n = [sum(bool(w0[(w0.gene == g) & (w0.rep == r)].peak_in_rep.iloc[0]) for r in reps) for g in genes_show]
ax.bar(range(len(genes_show)), n, color=["#8e44ad" if g in POS or g == "PNPLA3" else "#7f8c8d" for g in genes_show])
ax.set_xticks(range(len(genes_show))); ax.set_xticklabels(genes_show, rotation=75, fontsize=7); ax.set_yticks(range(len(reps) + 1))
ax.set_ylabel("replicates with a peak within TSS$\\pm$2 kb"); ax.set_title("h  reproducibility", loc="left", fontsize=9)
# panel h: peak presence in the TSS+/-2 kb window (WIN['promoter']); consensus 'bound' for the integration uses TSS+/-3 kb

ax = fig.add_subplot(gs[1, 4:6]); H = np.array([[np.log2(w0[(w0.gene == g) & (w0.rep == r)].RPGC.iloc[0] + 1) for r in reps] for g in genes_show])
im = ax.imshow(H, aspect="auto", cmap="viridis"); ax.set_yticks(range(len(genes_show))); ax.set_yticklabels(genes_show, fontsize=7)
ax.set_xticks(range(len(reps))); ax.set_xticklabels([r[-7:] for r in reps], rotation=45, fontsize=7); plt.colorbar(im, ax=ax, label="log2(RPGC+1)")
ax.set_title("k  raw signal (TSS$\\pm$2 kb)", loc="left", fontsize=9)

if rna is not None:
    ax = fig.add_subplot(gs[2, 0:2])
    for flag, col, lab in [(False, "gray", f"not bound (n={int((~rna.bound).sum())})"), (True, "#8e44ad", f"HepG2 PPARG-bound (n={int(rna.bound.sum())})")]:
        v = np.sort(rna[rna.bound == flag].stat.values); ax.plot(v, np.arange(1, len(v) + 1) / len(v), color=col, label=lab)
    if "PNPLA3" in set(rna.sym): ax.axvline(rna[rna.sym == "PNPLA3"].stat.iloc[0], ls="--", color="k", lw=.8)
    ax.set_xlim(-15, 15); ax.set_xlabel("TZD response t-statistic (mouse liver GSE162276)"); ax.set_ylabel("cumulative fraction"); ax.legend(fontsize=7)
    ax.set_title("i  TZD response vs binding (dashed = PNPLA3)", loc="left", fontsize=9)
    ax = fig.add_subplot(gs[2, 2:4]); pb, pu = 100 * a_ / max(a_ + b_, 1), 100 * c_ / max(c_ + d_, 1)
    ax.bar([0, 1], [pb, pu], color=["#8e44ad", "gray"]); ax.set_xticks([0, 1]); ax.set_xticklabels(["bound", "not bound"])
    ax.set_ylabel("% induced (padj<0.05)"); ax.set_title(f"j  OR={orr:.2f}, Fisher P={pf:.2g}", loc="left", fontsize=9)

ax = fig.add_subplot(gs[2, 4:6]); pn = met[met.gene == "PNPLA3"]; wl = list(WIN); wd = .35
for i, r in enumerate(reps):
    ax.bar(np.arange(len(wl)) + (i - .5) * wd, [pn[(pn.window == w) & (pn.rep == r)].log2_ChIP_over_input.iloc[0] for w in wl], wd, color=cols[i], label=r)
ax.axhline(0, color="k", lw=.6); ax.set_xticks(range(len(wl))); ax.set_xticklabels(wl, fontsize=8); ax.set_ylabel("mean log2(ChIP/input)")
ax.set_title("l  PNPLA3 windows", loc="left", fontsize=9); ax.legend(fontsize=6)
fig.suptitle("HepG2 PPARG ChIP-seq (ENCSR130VQL, hg38; control subsampled)", fontsize=11)
fig.savefig(f"{OUT}/hepg2_PPARG_figure.pdf", bbox_inches="tight"); fig.savefig(f"{OUT}/hepg2_PPARG_figure.png", dpi=150, bbox_inches="tight")
print("Wrote", f"{OUT}/hepg2_PPARG_figure.pdf", "and hepg2_gene_metrics.csv, hepg2_bound_genes_consensus.csv, hepg2_integration_and_QC.txt")

#!/usr/bin/env python3
"""HepG2 PPARG ChIP-seq (ENCODE ENCSR130VQL): the full Fig. 6-style pipeline in one script.

Sections (each is wrapped so that a failure is reported and the rest still runs):
  A  peak QC            peak numbers (replicates q<0.01 / q<0.05, pooled), widths, genomic annotation,
                        replicate concordance (peak overlap, genome-wide signal correlation)
  B  input-based tests  at PNPLA3 windows and control promoters: read counts ChIP vs input, conditional binomial
                        test, enrichment over input, genome-wide percentile of the window signal, percentile of the
                        INPUT signal (hyper-ChIPable check), peak calls in replicates and pooled
  C  gene-level binding + integration with RNA-seq responses (liver TZD GSE162276; hASC rosiglitazone) and with the
                        hASC PPARG-bound gene set (GSE115421)
  D  Hallmark over-representation of HepG2-bound genes
  E  PPARG::RXRA motif scan (JASPAR MA0065.1) at named sites, in all pooled peaks, against GC-matched random windows
  F  figure + text report

Environment variables
  OUT          results/hepg2_ppar
  GTF          GENCODE hg38 GTF(.gz)                          (required)
  FASTA        ref/hg38.analysisSet.fa.gz (for the motif scan; section E is skipped if missing)
  CTRL_BAM     results/hepg2_ppar/bam/control_sub20.bam
  DE_CSV       GSE162276_TZDvsControl_DE_results.csv  (GeneName, baseMean, log2FoldChange, stat, padj)
  HASC_CSV     ASCRosi_DE_with_ChIP_annotation.csv     (Gene, log2FC_Rosi_vs_DMSO, t_stat, padj, bound)
  HALLMARK_PKL hallmark_sets.pkl (pickled dict term -> list of human symbols)
Needs: pyBigWig pysam pandas numpy scipy matplotlib (statsmodels optional)
"""
import os, re, gzip, glob, sys, pickle, traceback
import numpy as np, pandas as pd, pyBigWig
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
try: import pysam
except ImportError: pysam = None
try: import statsmodels.api as sm
except ImportError: sm = None

OUT = os.environ.get("OUT", "results/hepg2_ppar")
GTF = os.environ["GTF"]
FASTA = os.environ.get("FASTA", "")
CTRL_BAM = os.environ.get("CTRL_BAM", f"{OUT}/bam/control_sub20.bam")
DE_CSV, HASC_CSV, HALLMARK_PKL = (os.environ.get(k, "") for k in ("DE_CSV", "HASC_CSV", "HALLMARK_PKL"))
OD = f"{OUT}/full_analysis"; os.makedirs(OD, exist_ok=True)
PROM = 3000
WIN = {"promoter": (-2000, 2000), "upstream_PPRE": (-850, -650), "intron1_site": (1000, 2300), "site_plus4.8kb": (4600, 5000)}
POS = ["FABP4", "ADIPOQ", "CD36", "PLIN1", "PCK1", "LPL", "ANGPTL4", "PLIN2", "PNPLA2"]
NEG = ["MYOD1", "RHO", "INS"]
rng = np.random.default_rng(7)
report, status = [], {}
def say(s=""): print(s); report.append(s)

def need(path, what):
    if not path or not os.path.exists(path):
        say(f"!! {what} not found: '{path}'"); return False
    return True

# ================================================================ inputs
reps = sorted(os.path.basename(p).replace(".RPGC.bw", "") for p in glob.glob(f"{OUT}/tracks/*.RPGC.bw")
              if not os.path.basename(p).startswith(("control", "ChIP_pooled"))
              and os.path.exists(f"{OUT}/peaks/" + os.path.basename(p).replace(".RPGC.bw", "") + "_peaks.narrowPeak"))
if len(reps) < 2: sys.exit(f"Need >=2 replicate tracks in {OUT}/tracks (found {reps})")
samples = reps + ["ChIP_pooled"]
for s in samples:
    for f in (f"{OUT}/tracks/{s}.log2ChIPoverInput.bw", f"{OUT}/peaks/{s}_peaks.narrowPeak", f"{OUT}/bam/{s}.bam"):
        if not os.path.exists(f): sys.exit(f"Missing {f}. Run steps 1 and 4 first.")
for f in (f"{OUT}/tracks/control.RPGC.bw", CTRL_BAM):
    if not os.path.exists(f): sys.exit(f"Missing {f}")
if pysam is None: sys.exit("pysam is required for the count-based tests:  pip install pysam   (or conda install -c bioconda pysam)")

def read_peaks(f):
    try: d = pd.read_csv(f, sep="\t", header=None, usecols=[0, 1, 2, 6, 8, 9], names=["chr", "s", "e", "fc", "nlq", "summit"])
    except pd.errors.EmptyDataError: d = pd.DataFrame(columns=["chr", "s", "e", "fc", "nlq", "summit"])
    d["pos"] = d.s + d.summit; return d
P = {}   # peak sets: key -> DataFrame
for s in samples:
    P[(s, "q01")] = read_peaks(f"{OUT}/peaks/{s}_peaks.narrowPeak")
    f05 = f"{OUT}/peaks/relaxed/{s}_q05_peaks.narrowPeak"
    P[(s, "q05")] = read_peaks(f05) if os.path.exists(f05) else None
    if P[(s, "q05")] is None: say(f"note: relaxed peaks missing for {s} ({f05}); run step 4")
trk = {s: pyBigWig.open(f"{OUT}/tracks/{s}.log2ChIPoverInput.bw") for s in samples}
rpgc = {s: pyBigWig.open(f"{OUT}/tracks/{s}.RPGC.bw") for s in samples}
ctrl_bw = pyBigWig.open(f"{OUT}/tracks/control.RPGC.bw")
for f in [f"{OUT}/bam/{s}.bam" for s in samples] + [CTRL_BAM]:
    if not (os.path.exists(f + ".bai") or os.path.exists(f[:-4] + ".bai")): pysam.index(f)
bam = {s: pysam.AlignmentFile(f"{OUT}/bam/{s}.bam") for s in samples}
cbam = pysam.AlignmentFile(CTRL_BAM)
libsize = lambda b: sum(x.mapped for x in b.get_index_statistics())
N = {s: libsize(bam[s]) for s in samples}; NC = libsize(cbam)
say("Library sizes (mapped, filtered, deduplicated): " + ", ".join(f"{s}={N[s]:,}" for s in samples) + f", control={NC:,}")

# ================================================================ annotation
op = gzip.open if GTF.endswith(".gz") else open
tx, ex, tss_rows = [], [], []
with op(GTF, "rt") as f:
    for ln in f:
        if ln[0] == "#": continue
        c = ln.rstrip("\n").split("\t")
        if len(c) < 9 or c[2] not in ("transcript", "exon"): continue
        m = re.search(r'gene_name "([^"]+)"', c[8])
        if not m: continue
        a, b = int(c[3]) - 1, int(c[4])
        if c[2] == "exon": ex.append((c[0], a, b))
        else:
            tx.append((c[0], a, b)); tss_rows.append((m.group(1), c[0], a if c[6] == "+" else b - 1, c[6]))
tss_df = pd.DataFrame(tss_rows, columns=["gene", "chr", "tss", "strand"]).drop_duplicates()
tss_df["key"] = np.where(tss_df.strand == "+", tss_df.tss, -tss_df.tss)
main = tss_df.sort_values("key").drop_duplicates("gene").set_index("gene")   # most upstream TSS per gene
say(f"PNPLA3 TSS (expected chr22:43,923,805, + strand): {main.loc['PNPLA3','chr']}:{main.loc['PNPLA3','tss']:,} {main.loc['PNPLA3','strand']}")

def merge_iv(rows):
    out = {}
    df = pd.DataFrame(rows, columns=["chr", "s", "e"])
    for ch, d in df.groupby("chr"):
        d = d.sort_values("s"); s, e = d.s.values, d.e.values; ms, me = [s[0]], [e[0]]
        for a, b in zip(s[1:], e[1:]):
            if a <= me[-1]: me[-1] = max(me[-1], b)
            else: ms.append(a); me.append(b)
        out[ch] = (np.array(ms), np.array(me))
    return out
def inside(iv, chs, pos):
    r = np.zeros(len(pos), bool); chs = np.asarray(chs); pos = np.asarray(pos)
    for ch in np.unique(chs):
        if ch not in iv: continue
        k = chs == ch; s, e = iv[ch]; i = np.searchsorted(s, pos[k], "right") - 1
        r[k] = (i >= 0) & (pos[k] < e[np.clip(i, 0, None)])
    return r
def overlaps(df, ch, a, b): return df[(df.chr == ch) & (df.s < b) & (df.e > a)]
def merged_of(df): return merge_iv(df[["chr", "s", "e"]].values.tolist()) if len(df) else {}
def ov_mask(A, B):
    """boolean mask: which rows of A overlap any interval of B (vectorised)"""
    iv = merged_of(B); m = np.zeros(len(A), bool); chs = A.chr.values
    for ch in np.unique(chs):
        if ch not in iv: continue
        k = chs == ch; ms, me = iv[ch]; i = np.searchsorted(ms, A.e.values[k], "right") - 1
        m[k] = (i >= 0) & (me[np.clip(i, 0, None)] > A.s.values[k])
    return m
IV_PROM = merge_iv([(r.chr, r.tss - PROM, r.tss + PROM) for r in tss_df.itertuples()])
IV_EXON, IV_TX = merge_iv(ex), merge_iv(tx)

def site(gene, w):
    t = main.loc[gene]; sg = 1 if t.strand == "+" else -1; a, b = WIN[w]
    lo, hi = sorted((int(t.tss) + sg * a, int(t.tss) + sg * b)); return t.chr, max(lo, 0), hi

# genome tiling helpers (percentiles) ------------------------------------------------------------
CHROMS = [c for c in rpgc[samples[0]].chroms() if re.fullmatch(r"chr([0-9]+|X|Y)", c)]
_tile = {}
TILE_FRACTION = float(os.environ.get("TILE_FRACTION", "0.1"))   # fraction of the genome used for the percentile background
SEG = 2_000_000
def tiles(bw_key, bw, width):
    """per-tile mean signal over a random, fixed subset (TILE_FRACTION) of 2-Mb segments of every chromosome.
    The subset is deterministic, so different tracks use the same tiles (needed for the replicate correlation)."""
    import time; k = (bw_key, width)
    if k not in _tile:
        t0 = time.time(); L = bw.chroms(); v = []
        for ci, ch in enumerate(CHROMS):
            nseg = int(np.ceil(L[ch] / SEG)); r = np.random.default_rng([12345, width, ci])
            pick = np.sort(r.choice(nseg, max(1, int(round(TILE_FRACTION * nseg))), replace=False)) if TILE_FRACTION < 1 else np.arange(nseg)
            for sg in pick:
                a0 = sg * SEG; a0 -= a0 % width; b0 = min(a0 + SEG, L[ch]); n = (b0 - a0) // width
                if n < 1: continue
                x = bw.stats(ch, a0, a0 + n * width, type="mean", nBins=n)
                v.append(np.array([np.nan if y is None else y for y in x], float))
        _tile[k] = np.concatenate(v); print(f"   [tiles] {bw_key} width={width}: {len(_tile[k]):,} tiles in {time.time()-t0:.0f}s", flush=True)
    return _tile[k]
def pct(bw_key, bw, width, value):
    t = tiles(bw_key, bw, width); t = t[~np.isnan(t)]; return 100 * float((t < value).mean()) if len(t) else np.nan
def wmean(bw, ch, a, b):
    v = bw.stats(ch, int(a), int(b), type="mean")[0]; return np.nan if v is None else v

# ================================================================ A. peak QC
fig_data = {}
try:
    say("\n=== A. Peak QC ===")
    rows = []
    for (s, q), d in P.items():
        if d is None: continue
        pk = d.copy(); pk["summit_pos"] = pk.pos
        ann = np.where(inside(IV_PROM, pk.chr, pk.summit_pos), "promoter (TSS+/-3 kb)",
              np.where(inside(IV_EXON, pk.chr, pk.summit_pos), "exon",
              np.where(inside(IV_TX, pk.chr, pk.summit_pos), "intron", "intergenic")))
        for cat in ("promoter (TSS+/-3 kb)", "exon", "intron", "intergenic"):
            rows.append(dict(sample=s, q=q, n_peaks=len(d), feature=cat, n=int((ann == cat).sum()), frac=float((ann == cat).mean()) if len(d) else np.nan,
                             median_width=float((d.e - d.s).median()) if len(d) else np.nan))
    qc = pd.DataFrame(rows); qc.to_csv(f"{OD}/peak_QC_annotation.csv", index=False); fig_data["qc"] = qc
    say(qc.groupby(["sample", "q"]).agg(n_peaks=("n_peaks", "first"), median_width=("median_width", "first")).to_string())
    # replicate concordance
    r1, r2 = reps[:2]; conc = []
    for qa in ("q01", "q05"):
        A, B = P[(r1, "q01")], P[(r2, qa)]
        if B is None: continue
        fa = float(ov_mask(A, B).mean()) if len(A) else np.nan
        fb = float(ov_mask(B, A).mean()) if len(B) else np.nan
        conc.append(dict(set=f"{r1} q<0.01 peaks overlapping {r2} {qa}", fraction=fa)); conc.append(dict(set=f"{r2} {qa} peaks overlapping {r1} q01", fraction=fb))
    say(pd.DataFrame(conc).to_string(index=False)); pd.DataFrame(conc).to_csv(f"{OD}/replicate_peak_concordance.csv", index=False)
    a_, b_ = tiles(r1, trk[r1], 1000), tiles(r2, trk[r2], 1000); ok = ~np.isnan(a_) & ~np.isnan(b_)
    pr, sr = stats.pearsonr(a_[ok], b_[ok])[0], stats.spearmanr(a_[ok], b_[ok])[0]
    say(f"log2(ChIP/input) in 1-kb tiles (random {TILE_FRACTION:.0%} of the genome), replicate vs replicate: Pearson r={pr:.3f}, Spearman rho={sr:.3f} (n={ok.sum():,})")
    idx = rng.choice(np.where(ok)[0], min(30000, ok.sum()), replace=False); fig_data["corr"] = (a_[idx], b_[idx], pr, sr); status["A"] = "ok"
except Exception as e: status["A"] = "FAILED"; say("A failed: " + repr(e)); traceback.print_exc()

# ================================================================ B. input-based tests
try:
    say("\n=== B. Input-based enrichment tests ==="); print("   (percentile backgrounds are computed once per track and width; progress is printed)", flush=True)
    tasks = [("PNPLA3", w) for w in WIN] + [(g, "promoter") for g in POS + NEG if g in main.index]
    rows = []
    for g, w in tasks:
        ch, lo, hi = site(g, w); width = hi - lo
        kc = cbam.count(ch, lo, hi)
        for s in samples:
            k = bam[s].count(ch, lo, hi); p0 = N[s] / (N[s] + NC)
            pv = stats.binomtest(k, k + kc, p0, alternative="greater").pvalue if k + kc > 0 else 1.0
            enr = ((k + .5) / N[s]) / ((kc + .5) / NC)
            row = dict(gene=g, window=w, chr=ch, start=lo, end=hi, sample=s, chip_reads=k, input_reads=kc, fold_enrichment=enr,
                       log2_enrichment=np.log2(enr), binom_p=pv, log2_ratio_track=wmean(trk[s], ch, lo, hi), RPGC=wmean(rpgc[s], ch, lo, hi))
            for q in ("q01", "q05"):
                d = P[(s, q)]
                if d is None: row[f"peak_{q}"], row[f"nlq_{q}"] = np.nan, np.nan; continue
                o = overlaps(d, ch, lo, hi); row[f"peak_{q}"] = len(o) > 0; row[f"nlq_{q}"] = float(o.nlq.max()) if len(o) else np.nan
            rows.append(row)
        # tile percentiles use the pooled-track-wide distribution of the same window width
        for r_ in rows[-len(samples):]:
            r_["input_RPGC"] = wmean(ctrl_bw, ch, lo, hi)
            if width < 300:      # 200-bp tiles would be ~15 million per track: skip the percentile
                r_["percentile_log2ratio"] = r_["input_percentile"] = np.nan; continue
            r_["percentile_log2ratio"] = pct(r_["sample"], trk[r_["sample"]], width, r_["log2_ratio_track"])
            r_["input_percentile"] = pct("ctrl", ctrl_bw, width, r_["input_RPGC"])
    W = pd.DataFrame(rows); W.to_csv(f"{OD}/window_enrichment_tests.csv", index=False); fig_data["W"] = W
    show = W[W.gene == "PNPLA3"][["window", "sample", "chip_reads", "input_reads", "fold_enrichment", "binom_p", "percentile_log2ratio",
                                  "input_percentile", "peak_q01", "peak_q05"]]
    say(show.to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    status["B"] = "ok"
except Exception as e: status["B"] = "FAILED"; say("B failed: " + repr(e)); traceback.print_exc()

# ================================================================ C. gene-level binding + integration
def bound_set(df_peaks_list, mode):
    """mode 'all': peak within TSS+/-3 kb in every set; df_peaks_list = list of peak DataFrames"""
    ivs = [merged_of(d) for d in df_peaks_list]
    def hit(iv, ch, a, b):
        if ch not in iv: return False
        s, e = iv[ch]; i = np.searchsorted(s, b, "right") - 1; return i >= 0 and e[i] > a
    return {g for g, r in main.iterrows() if all(hit(iv, r.chr, r.tss - PROM, r.tss + PROM) for iv in ivs)}
BOUND = {}
try:
    say("\n=== C. Gene-level binding and integration ===")
    pool = P[("ChIP_pooled", "q01")]
    # consensus of replicates (q<0.01): rep1 peaks overlapping any rep2 peak
    a_, b_ = P[(reps[0], "q01")], P[(reps[1], "q01")]
    cons = a_[ov_mask(a_, b_)]
    # high-confidence = pooled q<0.01 peaks supported by q<0.05 peaks in BOTH replicates
    rel = [P[(r, "q05")] for r in reps]
    if all(x is not None for x in rel):
        keep = np.ones(len(pool), bool)
        for x in rel: keep &= ov_mask(pool, x)
        hc = pool[keep]
    else: hc = pool.iloc[0:0]
    BOUND = {"consensus_q01_both_replicates": bound_set([cons], "all"), "pooled_q01": bound_set([pool], "all"),
             "pooled_supported_by_both_reps_q05": bound_set([hc], "all")}
    pd.DataFrame([dict(definition=k, n_genes=len(v), PNPLA3=("PNPLA3" in v), **{g: (g in v) for g in POS + NEG}) for k, v in BOUND.items()]).to_csv(f"{OD}/bound_gene_definitions.csv", index=False)
    for k, v in BOUND.items():
        say(f"{k}: {len(v):,} genes; positive controls bound {sum(g in v for g in POS)}/{len(POS)}; negative {sum(g in v for g in NEG)}/{len(NEG)}; PNPLA3={'PNPLA3' in v}")
        pd.DataFrame({"gene": sorted(v)}).to_csv(f"{OD}/bound_genes_{k}.csv", index=False)
    status["C-bound"] = "ok"
except Exception as e: status["C-bound"] = "FAILED"; say("C (bound sets) failed: " + repr(e)); traceback.print_exc()

def integrate(bset, de, label, lfc, stat, padj, base=None):
    de = de.dropna(subset=[padj, stat]).copy(); de["bound"] = de.sym.isin(bset); de["induced"] = (de[padj] < .05) & (de[lfc] > 0)
    a, b = int((de.bound & de.induced).sum()), int((de.bound & ~de.induced).sum()); c, d = int((~de.bound & de.induced).sum()), int((~de.bound & ~de.induced).sum())
    orr, pf = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
    mw = stats.mannwhitneyu(de.loc[de.bound, stat], de.loc[~de.bound, stat], alternative="greater").pvalue if de.bound.sum() > 2 else np.nan
    out = dict(dataset=label, n_tested=len(de), n_bound=int(de.bound.sum()), bound_induced=a, bound_not=b, unbound_induced=c, unbound_not=d,
               pct_bound_induced=100 * a / max(a + b, 1), pct_unbound_induced=100 * c / max(c + d, 1), fisher_OR=orr, fisher_p=pf, mannwhitney_p=mw,
               adj_logit_OR=np.nan, adj_logit_p=np.nan)
    if sm is not None and base and base in de and de.bound.sum() > 5:
        try:
            X = sm.add_constant(pd.DataFrame({"bound": de.bound.astype(float), "logbase": np.log10(de[base].clip(lower=1))}))
            m = sm.Logit(de.induced.astype(float), X).fit(disp=0); out["adj_logit_OR"], out["adj_logit_p"] = float(np.exp(m.params["bound"])), float(m.pvalues["bound"])
        except Exception as e: say(f"   logit failed for {label}: {e}")
    return out, de

RESP = {}
try:
    int_rows = []
    primary = "pooled_supported_by_both_reps_q05" if len(BOUND.get("pooled_supported_by_both_reps_q05", [])) else "pooled_q01"
    say(f"primary bound-gene definition for integration/Hallmark: {primary}")
    if need(DE_CSV, "DE_CSV (liver TZD)"):
        de = pd.read_csv(DE_CSV); de["sym"] = de.GeneName.str.upper(); de = de.sort_values("baseMean", ascending=False).drop_duplicates("sym")
        RESP["liver TZD (mouse, GSE162276)"] = (de, "log2FoldChange", "stat", "padj", "baseMean")
        say(f"liver DE: Pnpla3 baseMean={de.loc[de.sym=='PNPLA3','baseMean'].iloc[0]:.0f}, padj={de.loc[de.sym=='PNPLA3','padj'].iloc[0]:.3g}" if (de.sym == "PNPLA3").any() else "liver DE: PNPLA3 absent")
    if need(HASC_CSV, "HASC_CSV (hASC rosiglitazone)"):
        h = pd.read_csv(HASC_CSV).drop_duplicates("Gene"); h["sym"] = h.Gene.str.upper()
        RESP["hASC rosiglitazone (human, GSE115421)"] = (h, "log2FC_Rosi_vs_DMSO", "t_stat", "padj", None)
    for k, bs in BOUND.items():
        for lab, (d_, lfc, st, pa, ba) in RESP.items():
            o, _ = integrate({x.upper() for x in bs}, d_, lab, lfc, st, pa, ba); o["bound_definition"] = k; int_rows.append(o)
    if int_rows:
        I = pd.DataFrame(int_rows); I.to_csv(f"{OD}/integration_with_RNAseq.csv", index=False)
        say(I[["bound_definition", "dataset", "n_bound", "pct_bound_induced", "pct_unbound_induced", "fisher_OR", "fisher_p", "adj_logit_p"]].to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    # overlap with hASC PPARG-bound genes (GSE115421): universe = genes in the hASC table with a TSS in the GTF
    if "hASC rosiglitazone (human, GSE115421)" in RESP:
        h = RESP["hASC rosiglitazone (human, GSE115421)"][0]
        hb = set(h.loc[h.bound.astype(bool), "sym"]); uni = set(h.sym) & {g.upper() for g in main.index}
        ov_rows = []
        for k, bs in BOUND.items():
            hb_u = hb & uni; hep = {x.upper() for x in bs} & uni; ov = hb_u & hep
            M, n, Nn = len(uni), len(hb_u), len(hep); p = stats.hypergeom.sf(len(ov) - 1, M, n, Nn) if Nn else np.nan
            ov_rows.append(dict(bound_definition=k, universe=M, hASC_bound=n, HepG2_bound=Nn, overlap=len(ov), expected=n * Nn / M if M else np.nan,
                                hypergeom_p=p, PNPLA3_in_hASC=("PNPLA3" in hb_u), PNPLA3_in_HepG2=("PNPLA3" in hep)))
        O = pd.DataFrame(ov_rows); O.to_csv(f"{OD}/overlap_with_hASC_bound_genes.csv", index=False); fig_data["overlap"] = O
        say(O.to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    fig_data["RESP"] = RESP; fig_data["primary"] = primary; status["C-integration"] = "ok" if RESP else "skipped (no RNA-seq tables)"
except Exception as e: status["C-integration"] = "FAILED"; say("C (integration) failed: " + repr(e)); traceback.print_exc()

# ================================================================ D. Hallmark ORA
try:
    say("\n=== D. Hallmark over-representation ===")
    hp = HALLMARK_PKL
    if not need(hp, "HALLMARK_PKL"): raise RuntimeError("no Hallmark gene sets")
    H = pickle.load(open(hp, "rb"))
    bs = {x.upper() for x in BOUND[fig_data.get("primary", "pooled_q01")]}
    if "liver TZD (mouse, GSE162276)" in RESP:
        d_ = RESP["liver TZD (mouse, GSE162276)"][0]; bg = set(d_.loc[d_.baseMean >= 10, "sym"]); bgname = "genes with baseMean>=10 in liver RNA-seq"
    else: bg = {g.upper() for g in main.index}; bgname = "all GTF genes"
    bsb = bs & bg; res = []
    for term, genes in H.items():
        K = set(genes) & bg
        if len(K) < 5: continue
        k = len(K & bsb); res.append(dict(Term=term, K=len(K), k=k, expected=len(K) * len(bsb) / len(bg),
                                          pval=stats.hypergeom.sf(k - 1, len(bg), len(K), len(bsb)), genes=";".join(sorted(K & bsb))))
    R = pd.DataFrame(res).sort_values("pval")
    o = R.pval.values; m = len(o); q = np.minimum.accumulate((o * m / np.arange(1, m + 1))[::-1])[::-1]; R["padj"] = np.minimum(q, 1)
    R.to_csv(f"{OD}/Hallmark_ORA_HepG2_bound_genes.csv", index=False); fig_data["hall"] = R
    say(f"background: {bgname} (n={len(bg):,}); bound genes in background: {len(bsb):,}")
    say(R[["Term", "K", "k", "expected", "padj"]].head(10).to_string(index=False, float_format=lambda x: f"{x:.3g}")); status["D"] = "ok"
except Exception as e: status["D"] = "skipped/FAILED: " + str(e); say("D: " + repr(e))

# ================================================================ E. motif scan
PFM = {"A": [10, 7, 19, 13, 17, 0, 0, 0, 2, 41, 41, 39, 0, 0, 0, 2, 38, 10, 13, 7], "C": [12, 13, 4, 1, 0, 0, 0, 2, 37, 0, 0, 0, 0, 0, 2, 36, 1, 12, 15, 10],
       "G": [7, 8, 6, 27, 24, 41, 39, 0, 2, 0, 0, 2, 41, 39, 0, 1, 2, 4, 9, 12], "T": [12, 13, 12, 0, 0, 0, 2, 39, 0, 0, 0, 0, 0, 2, 39, 2, 0, 15, 4, 12]}
BG = {"A": .29, "C": .21, "G": .21, "T": .29}
pwm = np.zeros((20, 4))
for p in range(20):
    tot = sum(PFM[b][p] for b in "ACGT")
    for j, b in enumerate("ACGT"): pwm[p, j] = np.log2(((PFM[b][p] + .8 * BG[b]) / (tot + .8)) / BG[b])
pwm_rc = pwm[::-1, ::-1]             # reverse complement (A<->T, C<->G, position reversed)
pwm_max = pwm.max(1).sum()
def best_score(codes):               # codes: (R, L) ints 0-3, 4 = N ; returns best %max over both strands and positions
    R, L = codes.shape; n = L - 20 + 1; sc = np.zeros((R, n)); sc2 = np.zeros((R, n)); bad = np.zeros((R, n), bool)
    pw = np.vstack([pwm.T, np.full((1, 20), -9.0)]); pr = np.vstack([pwm_rc.T, np.full((1, 20), -9.0)])   # (5,20), N penalised
    for p in range(20):
        seg = codes[:, p:p + n]; sc += pw[seg, p]; sc2 += pr[seg, p]
    return 100 * np.maximum(sc, sc2).max(1) / pwm_max
def fasta_chroms(path, wanted):
    op_ = gzip.open if path.endswith(".gz") else open; name, buf = None, []
    with op_(path, "rt") as f:
        for ln in f:
            if ln[0] == ">":
                if name in wanted and buf: yield name, "".join(buf)
                name = ln[1:].split()[0]; buf = []
            elif name in wanted: buf.append(ln.strip())
        if name in wanted and buf: yield name, "".join(buf)
LUT = np.full(256, 4, np.uint8)
for i, b in enumerate("ACGT"): LUT[ord(b)] = i; LUT[ord(b.lower())] = i
try:
    say("\n=== E. PPARG::RXRA motif scan (JASPAR MA0065.1) ===")
    if not need(FASTA, "FASTA"): raise RuntimeError("no FASTA")
    L_ = 300; half = L_ // 2
    named = {}
    for g, w in [("PNPLA3", "intron1_site"), ("PNPLA3", "site_plus4.8kb"), ("PNPLA3", "upstream_PPRE"), ("PLIN1", "promoter"), ("PCK1", "promoter")]:
        if g not in main.index: continue
        ch, lo, hi = site(g, w); o = overlaps(pool, ch, lo, hi)
        if w == "upstream_PPRE": c = (lo + hi) // 2; lab = f"{g}_{w}(window,no peak)"
        elif len(o): r = o.sort_values("nlq").iloc[-1]; c = int(r.pos); lab = f"{g}_{w}(pooled peak summit)"
        else: c = (lo + hi) // 2; lab = f"{g}_{w}(window centre,no peak)"
        named[lab] = (ch, max(c - half, 0), max(c - half, 0) + L_)
    top = pool.sort_values("nlq", ascending=False).head(20000).reset_index(drop=True)
    want = set(top.chr) | {v[0] for v in named.values()}
    pk_codes, pk_gc, rnd_codes, named_codes = [], [], [], {}
    ar = np.arange(L_)
    for ch, seq in fasta_chroms(FASTA, want):
        arr = LUT[np.frombuffer(seq.encode(), np.uint8)]; Lc = len(arr)
        sub = top[top.chr == ch]
        if len(sub):
            st = sub.pos.values.astype(int) - half; st = st[(st >= 0) & (st + L_ <= Lc)]
            if len(st):
                w_ = arr[st[:, None] + ar]; okw = ~(w_ == 4).any(1); w_ = w_[okw]
                pk_codes.append(w_); gc_p = np.isin(w_, (1, 2)).mean(1); pk_gc.append(gc_p)
                cs = rng.integers(0, max(Lc - L_, 1), size=len(st) * 20 + 50); wc = arr[cs[:, None] + ar]; okc = ~(wc == 4).any(1); wc = wc[okc]
                gcc = np.isin(wc, (1, 2)).mean(1); order = np.argsort(gcc); sg_ = gcc[order]
                for off in (-1, 0, 1):
                    pos_ = np.clip(np.searchsorted(sg_, gc_p) + off, 0, len(sg_) - 1); rnd_codes.append(wc[order[pos_]])
        for lab, (c2, a, b) in named.items():
            if c2 == ch: w_ = arr[a:b]; named_codes[lab] = w_ if len(w_) == L_ else None
        del arr
    pk_codes = np.concatenate(pk_codes); rnd_codes = np.concatenate(rnd_codes)
    PK = np.concatenate([best_score(pk_codes[i:i + 5000]) for i in range(0, len(pk_codes), 5000)])
    RN = np.concatenate([best_score(rnd_codes[i:i + 5000]) for i in range(0, len(rnd_codes), 5000)])
    thr = np.percentile(RN, 95); frac_pk = float((PK >= thr).mean()); mw = stats.mannwhitneyu(PK, RN, alternative="greater").pvalue
    say(f"{len(PK):,} pooled peaks (top by -log10 q) vs {len(RN):,} GC-matched random {L_}-bp windows; best-site score %max: peaks median {np.median(PK):.1f}, random median {np.median(RN):.1f}")
    say(f"fraction of peaks above the random 95th percentile ({thr:.1f}%): {frac_pk:.3f} (expected 0.05; {frac_pk/0.05:.1f}-fold); Mann-Whitney p={mw:.2g}")
    nm_rows = []
    for lab, c_ in named_codes.items():
        if c_ is None: continue
        s_ = float(best_score(c_[None, :])[0]); nm_rows.append(dict(site=lab, best_score_pct_max=s_, percentile_vs_random=100 * float((RN < s_).mean()), percentile_vs_peaks=100 * float((PK < s_).mean())))
    NM = pd.DataFrame(nm_rows); NM.to_csv(f"{OD}/motif_scan_named_sites.csv", index=False); say(NM.to_string(index=False, float_format=lambda x: f"{x:.1f}"))
    fig_data["motif"] = (PK, RN, NM, thr, frac_pk); status["E"] = "ok"
except Exception as e: status["E"] = "skipped/FAILED: " + str(e); say("E: " + repr(e)); traceback.print_exc()

# ================================================================ F. figure
def panel(ax, fn, title):
    try: fn(ax); ax.set_title(title, loc="left", fontsize=9, fontweight="bold")
    except Exception as e: ax.text(.5, .5, f"not available\n({str(e)[:60]})", ha="center", va="center", fontsize=7, transform=ax.transAxes); ax.set_title(title, loc="left", fontsize=9)
colors = {s: c for s, c in zip(samples, ["#4c6a9c", "#c0504d", "#2d2d2d"])}
fig = plt.figure(figsize=(17, 19)); gs = GridSpec(4, 4, figure=fig, hspace=.6, wspace=.4, top=.91, bottom=.05)
fig.suptitle("HepG2 PPARG ChIP-seq (ENCSR130VQL, hg38): pooled-replicate and input-based analysis", fontsize=12, y=.935)

def pA(ax):
    qc = fig_data["qc"]; t = qc.drop_duplicates(["sample", "q"]); lab = [f"{s}\n{q.replace('q0', 'q<0.0')}" for s, q in zip(t["sample"], t.q)]
    ax.bar(range(len(t)), t.n_peaks, color=[colors[s] for s in t["sample"]]); ax.set_xticks(range(len(t))); ax.set_xticklabels(lab, fontsize=6, rotation=60); ax.set_ylabel("peaks")
def pB(ax):
    for s in samples:
        d = P[(s, "q01")]; w = np.sort((d.e - d.s).values); ax.plot(w, np.arange(1, len(w) + 1) / max(len(w), 1), color=colors[s], label=s)
    ax.set_xlabel("peak width (bp), q<0.01"); ax.set_ylabel("ECDF"); ax.legend(fontsize=6)
def pC(ax):
    qc = fig_data["qc"]; cats = ["promoter (TSS+/-3 kb)", "exon", "intron", "intergenic"]; t = qc[qc.q == "q01"].pivot(index="sample", columns="feature", values="frac")[cats]
    bottom = np.zeros(len(t))
    for c, col in zip(cats, ["#8e44ad", "#e67e22", "#3498db", "#95a5a6"]): ax.bar(range(len(t)), t[c].values, bottom=bottom, label=c, color=col); bottom += t[c].values
    ax.set_xticks(range(len(t))); ax.set_xticklabels(t.index, fontsize=6, rotation=30); ax.set_ylabel("fraction of peak summits"); ax.legend(fontsize=5, loc="upper center", bbox_to_anchor=(.5, -.18), ncol=2)
def pD(ax):
    a, b, pr, sr = fig_data["corr"]; ax.hexbin(a, b, gridsize=60, bins="log", cmap="viridis"); ax.set_xlabel(f"{reps[0]} log2(ChIP/input), 1-kb tiles", fontsize=7)
    ax.set_ylabel(f"{reps[1]}", fontsize=7); ax.text(.05, .92, f"r={pr:.2f}  rho={sr:.2f}", transform=ax.transAxes, fontsize=8, color="w")
panel(fig.add_subplot(gs[0, 0]), pA, "a  peak numbers"); panel(fig.add_subplot(gs[0, 1]), pB, "b  peak widths")
panel(fig.add_subplot(gs[0, 2]), pC, "c  genomic annotation"); panel(fig.add_subplot(gs[0, 3]), pD, "d  replicate concordance")

def pE(ax_unused):
    ch, tss, strand = main.loc["PNPLA3", "chr"], int(main.loc["PNPLA3", "tss"]), main.loc["PNPLA3", "strand"]
    a, b = tss - 2500, tss + 7500
    inner = GridSpecFromSubplotSpec(5, 1, subplot_spec=gs[1, :], hspace=.08, height_ratios=[1, 1, 1, 1, 1.2]); axs = []
    for i, (lab, bw, c) in enumerate([(s, trk[s], colors[s]) for s in samples] + [("input RPGC", ctrl_bw, "#7f8c8d")]):
        ax = fig.add_subplot(inner[i]); axs.append(ax); nb = 500; v = bw.stats(ch, a, b, nBins=nb, type="mean"); v = np.array([0 if x is None else x for x in v])
        ax.fill_between(np.linspace(a, b, nb), 0, v, color=c, step="mid"); ax.set_ylabel(("log2 " if i < 3 else "") + lab.replace("ChIP_", ""), fontsize=6, rotation=0, ha="right", va="center"); ax.set_xticks([])
        ax.axvline(tss, color="k", lw=.8)
        for w, col in (("upstream_PPRE", "tab:red"), ("intron1_site", "tab:green"), ("site_plus4.8kb", "tab:orange")):
            _, lo, hi = site("PNPLA3", w); ax.axvspan(lo, hi, color=col, alpha=.15)
    ax = fig.add_subplot(inner[4]); axs.append(ax)
    for j, (s, q) in enumerate([(reps[0], "q01"), (reps[0], "q05"), (reps[1], "q01"), (reps[1], "q05"), ("ChIP_pooled", "q01"), ("ChIP_pooled", "q05")]):
        d = P[(s, q)]
        if d is None: continue
        for r in overlaps(d, ch, a, b).itertuples(): ax.add_patch(plt.Rectangle((r.s, j + .1), r.e - r.s, .8, color=colors[s], alpha=.9 if q == "q01" else .45))
    ax.set_ylim(0, 6); ax.set_yticks(np.arange(6) + .5); ax.set_yticklabels([f"{s.replace('ChIP_','')} {q.replace('q0','q<0.0')}" for s, q in [(reps[0], "q01"), (reps[0], "q05"), (reps[1], "q01"), (reps[1], "q05"), ("ChIP_pooled", "q01"), ("ChIP_pooled", "q05")]], fontsize=5)
    ax.set_xlim(a, b); ax.set_xlabel(f"PNPLA3 {ch} (hg38). Black = TSS; red = upstream PPREs (-850/-650); green = intron-1 site (+1.0/+2.3 kb); orange = +4.8 kb peak", fontsize=7)
    axs[0].set_title("e  PNPLA3 locus: log2(ChIP/input), input and called peaks", loc="left", fontsize=9, fontweight="bold")
    for x in axs[:4]: x.set_xlim(a, b)
try: pE(None)
except Exception as e: fig.text(.1, .55, f"panel e failed: {e}", fontsize=8)

def pF(ax):
    W = fig_data["W"]; sites = [("PNPLA3", w) for w in WIN] + [(g, "promoter") for g in ("PLIN1", "PCK1", "FABP4", "MYOD1")]
    x = np.arange(len(sites)); wd = .26
    for i, s in enumerate(samples):
        d = [W[(W.gene == g) & (W.window == w) & (W["sample"] == s)].iloc[0] for g, w in sites]
        ax.bar(x + (i - 1) * wd, [r.log2_enrichment for r in d], wd, color=colors[s], label=s.replace("ChIP_", ""))
        for xi, r in zip(x, d): ax.text(xi + (i - 1) * wd, max(r.log2_enrichment, 0) + .05, "***" if r.binom_p < 1e-6 else "**" if r.binom_p < 1e-3 else "*" if r.binom_p < .05 else "", ha="center", fontsize=7)
    ax.axhline(0, color="k", lw=.5); ax.set_xticks(x); ax.set_xticklabels([f"{g}\n{w.replace('_',' ')}" if g == "PNPLA3" else g for g, w in sites], fontsize=6, rotation=45, ha="right")
    ax.set_ylabel("log2 (ChIP / input read density)"); ax.legend(fontsize=6)
def pG(ax):
    W = fig_data["W"]; sites = [("PNPLA3", "intron1_site"), ("PNPLA3", "site_plus4.8kb"), ("PNPLA3", "promoter"), ("PLIN1", "promoter"), ("PCK1", "promoter")]
    x = np.arange(len(sites)); wd = .26
    for i, s in enumerate(samples):
        v = [W[(W.gene == g) & (W.window == w) & (W["sample"] == s)].percentile_log2ratio.iloc[0] for g, w in sites]; ax.bar(x + (i - 1) * wd, v, wd, color=colors[s])
    iv = [W[(W.gene == g) & (W.window == w) & (W["sample"] == "ChIP_pooled")].input_percentile.iloc[0] for g, w in sites]
    ax.plot(x, iv, "k_", ms=18, mew=2, label="input percentile")
    ax.axhline(99, color="r", ls="--", lw=.8); ax.set_ylim(0, 105); ax.set_xticks(x); ax.set_xticklabels([f"{g}\n{w.replace('_',' ')}" for g, w in sites], fontsize=6, rotation=45, ha="right")
    ax.set_ylabel("genome-wide percentile\n(same-width tiles)"); ax.legend(fontsize=6)
def pH(ax):
    PK, RN, NM, thr, frac = fig_data["motif"]; bins = np.linspace(min(PK.min(), RN.min()), 100, 50)
    ax.hist(RN, bins, density=True, alpha=.5, color="#7f8c8d", label="GC-matched random"); ax.hist(PK, bins, density=True, alpha=.5, color="#8e44ad", label=f"pooled peaks (n={len(PK):,})")
    ax.axvline(thr, color="k", ls=":", lw=.8)
    for r, c in zip(NM.itertuples(), plt.cm.tab10.colors): ax.axvline(r.best_score_pct_max, color=c, lw=1.2, label=r.site.split("(")[0])
    ax.set_xlabel("best PPARG::RXRA match (% of max PWM score), 300-bp windows"); ax.legend(fontsize=5); ax.text(.02, .95, f"{frac/0.05:.1f}x peaks above random 95th pct", transform=ax.transAxes, fontsize=7, va="top")
panel(fig.add_subplot(gs[2, 0:2]), pF, "f  read-count enrichment over input (binomial test)"); panel(fig.add_subplot(gs[2, 2]), pG, "g  genome-wide percentile")
panel(fig.add_subplot(gs[2, 3]), pH, "h  motif scan")

def ecdf_panel(label):
    def f(ax):
        d_, lfc, st, pa, ba = fig_data["RESP"][label]; bs = {x.upper() for x in BOUND[fig_data["primary"]]}
        d_ = d_.dropna(subset=[st]); b_ = d_[d_.sym.isin(bs)][st].values; u_ = d_[~d_.sym.isin(bs)][st].values
        for v, c, lb in [(u_, "#7f8c8d", f"unbound (n={len(u_):,})"), (b_, "#8e44ad", f"HepG2-bound (n={len(b_):,})")]:
            v = np.sort(v); ax.plot(v, np.arange(1, len(v) + 1) / len(v), color=c, label=lb)
        ax.set_xlim(-6, 6); ax.set_xlabel("response statistic (t / Wald)"); ax.set_ylabel("ECDF"); ax.legend(fontsize=6)
        I = pd.read_csv(f"{OD}/integration_with_RNAseq.csv"); r = I[(I.dataset == label) & (I.bound_definition == fig_data["primary"])].iloc[0]
        ax.text(.03, .92, f"induced: bound {r.pct_bound_induced:.1f}% vs unbound {r.pct_unbound_induced:.1f}%\nFisher p={r.fisher_p:.2g}; adj. p={r.adj_logit_p:.2g}", transform=ax.transAxes, fontsize=6, va="top")
    return f
labs = list(fig_data.get("RESP", {}))
for i, lab in enumerate(labs[:2]): panel(fig.add_subplot(gs[3, i]), ecdf_panel(lab), f"{'ij'[i]}  {lab}")
def pK(ax):
    O = fig_data["overlap"]; r = O[O.bound_definition == fig_data["primary"]].iloc[0]
    ax.bar([0, 1], [r.expected, r.overlap], color=["#7f8c8d", "#8e44ad"]); ax.set_xticks([0, 1]); ax.set_xticklabels(["expected\nby chance", "observed"]); ax.set_ylabel("genes bound in both HepG2 and hASC")
    ax.text(.5, .9, f"HepG2-bound {r.HepG2_bound:,}; hASC-bound {r.hASC_bound:,}\nhypergeometric p={r.hypergeom_p:.2g}\nPNPLA3 bound: hASC={r.PNPLA3_in_hASC}, HepG2={r.PNPLA3_in_HepG2}", transform=ax.transAxes, ha="center", va="top", fontsize=6)
def pL(ax):
    R = fig_data["hall"].head(8)[::-1]; ax.barh([t.replace("HALLMARK_", "").replace("_", " ").title() for t in R.Term], -np.log10(R.padj.clip(lower=1e-300)), color="#8e44ad")
    ax.axvline(-np.log10(.05), color="r", ls="--", lw=.8); ax.set_xlabel("-log10 FDR"); ax.tick_params(axis="y", labelsize=6)
panel(fig.add_subplot(gs[3, 2]), pK, "k  overlap with hASC-bound genes"); panel(fig.add_subplot(gs[3, 3]), pL, "l  Hallmark, HepG2-bound genes")
fig.savefig(f"{OD}/hepg2_full_analysis_figure.pdf", bbox_inches="tight"); fig.savefig(f"{OD}/hepg2_full_analysis_figure.png", dpi=200, bbox_inches="tight")

# ================================================================ verdict
try:
    W = fig_data["W"]; say("\n=== VERDICT FOR PNPLA3 (automatic; read with the caveats in the text) ===")
    for w in ("intron1_site", "site_plus4.8kb", "promoter", "upstream_PPRE"):
        d = W[(W.gene == "PNPLA3") & (W.window == w)].set_index("sample")
        rep_ok = [bool(d.loc[r, "peak_q05"] == True) for r in reps]
        enr_ok = [bool(d.loc[r, "fold_enrichment"] >= 2 and d.loc[r, "binom_p"] < 1e-5) for r in reps]
        pooled = bool(d.loc["ChIP_pooled", "peak_q01"]); pc = d.loc["ChIP_pooled", "percentile_log2ratio"]; ip = d.loc["ChIP_pooled", "input_percentile"]
        call = "REPRODUCIBLE binding evidence" if all(rep_ok) and all(enr_ok) and pooled else \
               "PARTIAL: supported in some analyses only" if pooled or any(rep_ok) or any(enr_ok) else "NO significant enrichment over input"
        say(f"{w:16s}: {call}  [replicate peaks q<0.05 {rep_ok}; replicate enrichment>=2x & p<1e-5 {enr_ok}; pooled q<0.01 peak {pooled}; "
            f"pooled log2-ratio percentile {pc:.1f}; input percentile {ip:.1f}{' (INPUT ALSO HIGH: possible artefact)' if ip > 95 else ''}]")
    c = [g for g in POS + NEG if g in main.index]
    say("Controls bound by the primary definition (positive/negative): " + f"{sum(g in BOUND[fig_data['primary']] for g in POS)}/{len(POS)} / {sum(g in BOUND[fig_data['primary']] for g in NEG)}/{len(NEG)}")
except Exception as e: say("verdict failed: " + repr(e))
say("\nsection status: " + ", ".join(f"{k}={v}" for k, v in status.items()))
open(f"{OD}/hepg2_full_analysis_report.txt", "w").write("\n".join(report) + "\n")
print(f"\nDone. Outputs in {OD}/")

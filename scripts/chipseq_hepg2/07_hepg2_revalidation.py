#!/usr/bin/env python3
"""Step 7: revalidation of the HepG2 PPARG ChIP-seq evidence and of the binding-vs-response enrichment (Fig panels f-i).

Answers the points a reviewer will raise:
  1. peak-level statistics straight from the MACS3 output (fold change at summit, -log10 P, -log10 q, rank among all peaks)
     for every named window, per replicate and pooled - instead of tests on aggregated window counts
  2. an independent MACS-style local-background test (max of 400 bp / 1 kb / 5 kb / 10 kb input density and genome average)
  3. coverage / mappability / blacklist / optional extra tracks (DNase, H3K27ac...) at each window
  4. the binding-vs-response enrichment redone with 1:1 human-mouse orthologs on ONE common gene universe, with raw counts,
     odds ratios with 95% CI, two-sided P, an expression-adjusted model, and FORMAL tests that the two datasets differ
     (Woolf z, GEE interaction, gene-label permutation incl. expression-stratified), plus sensitivity analyses over
     alternative bound-gene definitions (promoter +-1 kb, +-3 kb [primary, prespecified], distal 3-50 kb nearest-TSS)
  5. effect sizes with 95% CI and equivalence (TOST) for PNPLA3 and the control genes in both RNA-seq datasets

Run from the project root after steps 1-5 (conda env active):
  OUT=results/hepg2_ppar GTF=ref/gencode.v46.basic.annotation.gtf.gz \
  DE_CSV=inputs/GSE162276_TZDvsControl_DE_results.csv HASC_CSV=inputs/ASCRosi_DE_with_ChIP_annotation.csv \
  python3 scripts/chipseq_hepg2/07_hepg2_revalidation.py

Optional environment variables
  ORTHO_FILE   MGI HOM_MouseHumanSequence.rpt  OR a CSV with columns mouse,human (1:1 orthologs). If absent the script tries to
               download https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt into $OUT/../ref
  MAP_BW       Umap multi-read mappability bigWig (e.g. k50.Umap.MultiTrackMappability.bw, hg38)
  BLACKLIST_BED  ENCODE blacklist v2 (hg38-blacklist.v2.bed[.gz])
  EXTRA_BW     name=path,name=path   e.g. HepG2_DNase=ENCFF....bigWig,HepG2_H3K27ac=ENCFF....bigWig (signal in each window)
  NPERM        permutations (default 10000)
  DE_RECOMP_CSV  DESeq2 re-run of the GSE162276 read counts (made by step 00); supplies the SE for the liver confidence intervals in section 5
                 (default: GSE162276_TZDvsControl_DE_recomputed_DESeq2.csv next to DE_CSV)
  REVAL_TAG    suffix for the output folder, for the sensitivity run: DE_CSV=<recomputed table> REVAL_TAG=_DESeq2rerun
Outputs in $OUT/full_analysis/revalidation/ (CSV tables + revalidation_report.txt)
"""
import os, re, gzip, glob, sys, traceback, urllib.request, warnings
import numpy as np, pandas as pd
from scipy import stats
warnings.filterwarnings("ignore")
try: import pysam
except ImportError: sys.exit("pysam is required:  pip install pysam")
import pyBigWig
try: import statsmodels.api as sm, statsmodels.formula.api as smf
except ImportError: sm = smf = None

OUT = os.environ.get("OUT", "results/hepg2_ppar"); GTF = os.environ["GTF"]
DE_CSV, HASC_CSV = os.environ.get("DE_CSV", ""), os.environ.get("HASC_CSV", "")
CTRL_BAM = os.environ.get("CTRL_BAM", f"{OUT}/bam/control_sub20.bam")
DE_RECOMP_CSV = os.environ.get("DE_RECOMP_CSV", os.path.join(os.path.dirname(DE_CSV), "GSE162276_TZDvsControl_DE_recomputed_DESeq2.csv") if DE_CSV else "")
TAG = os.environ.get("REVAL_TAG", "")      # e.g. REVAL_TAG=_DESeq2rerun for the sensitivity run with DE_CSV=<recomputed table>
OD = f"{OUT}/full_analysis"; RV = f"{OD}/revalidation{TAG}"; os.makedirs(RV, exist_ok=True)
NPERM = int(os.environ.get("NPERM", 10000)); PROM = 3000; GENOME = 2.913e9
rng = np.random.default_rng(7); report = []; status = {}
def say(s=""): print(s); report.append(s)
def need(f, what):
    if not f or not os.path.exists(f): say(f"!! {what} not found: '{f}'"); return False
    return True

# ---------------------------------------------------------------- inputs (same conventions as step 5)
reps = sorted(os.path.basename(p).replace(".RPGC.bw", "") for p in glob.glob(f"{OUT}/tracks/*.RPGC.bw")
              if not os.path.basename(p).startswith(("control", "ChIP_pooled"))
              and os.path.exists(f"{OUT}/peaks/" + os.path.basename(p).replace(".RPGC.bw", "") + "_peaks.narrowPeak"))
if len(reps) < 2: sys.exit(f"need >=2 replicates in {OUT}/tracks (found {reps})")
samples = reps + ["ChIP_pooled"]
for f in [f"{OUT}/bam/{s}.bam" for s in samples] + [CTRL_BAM, f"{OD}/window_enrichment_tests.csv", f"{OUT}/tracks/control.RPGC.bw"]:
    if not os.path.exists(f): sys.exit(f"Missing {f}. Run steps 1-5 first.")
NPC = ["chr", "s", "e", "name", "score", "strand", "fc", "nlp", "nlq", "summit"]
def read_peaks(f):
    if not os.path.exists(f): return None
    try: d = pd.read_csv(f, sep="\t", header=None, usecols=range(10), names=NPC)
    except pd.errors.EmptyDataError: d = pd.DataFrame(columns=NPC)
    d["pos"] = d.s + d.summit; return d
P = {}
for s in samples:
    P[(s, "q01")] = read_peaks(f"{OUT}/peaks/{s}_peaks.narrowPeak"); P[(s, "q05")] = read_peaks(f"{OUT}/peaks/relaxed/{s}_q05_peaks.narrowPeak")
for f in [f"{OUT}/bam/{s}.bam" for s in samples] + [CTRL_BAM]:
    if not (os.path.exists(f + ".bai") or os.path.exists(f[:-4] + ".bai")): pysam.index(f)
bam = {s: pysam.AlignmentFile(f"{OUT}/bam/{s}.bam") for s in samples}; cbam = pysam.AlignmentFile(CTRL_BAM)
libsize = lambda b: sum(x.mapped for x in b.get_index_statistics())
N = {s: libsize(bam[s]) for s in samples}; NC = libsize(cbam)
say("Library sizes: " + ", ".join(f"{s}={N[s]:,}" for s in samples) + f", control={NC:,}")
W = pd.read_csv(f"{OD}/window_enrichment_tests.csv")
WIN = W.drop_duplicates(["gene", "window"])[["gene", "window", "chr", "start", "end"]].reset_index(drop=True)

# ================================================================ 1. peak-level statistics from the MACS3 output
try:
    say("\n=== 1. Peak-level statistics (MACS3 narrowPeak; fold change at summit, -log10 P, -log10 q) ===")
    rows = []
    for s in samples:
        for q in ("q01", "q05"):
            pk = P[(s, q)]
            if pk is None or not len(pk): continue
            nl = pk.nlq.values
            for w in WIN.itertuples():
                hit = pk[(pk.chr == w.chr) & (pk.s < w.end) & (pk.e > w.start)]
                if hit.empty: rows.append(dict(gene=w.gene, window=w.window, sample=s, q_level=q, peak=False)); continue
                for r in hit.itertuples():
                    rk = 1 + int((nl > r.nlq).sum())
                    rows.append(dict(gene=w.gene, window=w.window, sample=s, q_level=q, peak=True, peak_chr=r.chr, peak_start=r.s, peak_end=r.e, width=r.e - r.s,
                                     summit=int(r.pos), fold_change_at_summit=r.fc, neg_log10_p=r.nlp, neg_log10_q=r.nlq, q_value=10.0 ** -r.nlq,
                                     rank_among_peaks=rk, n_peaks=len(pk), top_percent_of_peaks=100 * rk / len(pk)))
    PL = pd.DataFrame(rows); PL.to_csv(f"{RV}/peak_level_statistics.csv", index=False)
    sh = PL[(PL.peak) & (PL.q_level == "q01") & (PL["sample"] == "ChIP_pooled")][["gene", "window", "peak_start", "peak_end", "fold_change_at_summit", "neg_log10_q", "rank_among_peaks", "n_peaks"]]
    say("pooled MACS3 peaks (q<0.01) overlapping each window:\n" + sh.to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    rr = PL[(PL.peak) & (PL.q_level == "q05") & (PL.gene == "PNPLA3") & (PL["sample"] != "ChIP_pooled")][["window", "sample", "fold_change_at_summit", "neg_log10_q"]]
    say("per-replicate MACS3 peaks (q<0.05) at PNPLA3 windows:\n" + rr.to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    status["1"] = "ok"
except Exception as e: status["1"] = "FAILED"; say("1 failed: " + repr(e)); traceback.print_exc()

# ================================================================ 2. local-background (MACS-style) test at the summits
def local_test(s, ch, c):
    k = bam[s].count(ch, max(c - 200, 0), c + 200); sc = N[s] / NC; cand = {}
    for nm, wd in (("input_400bp", 400), ("input_1kb", 1000), ("input_5kb", 5000), ("input_10kb", 10000)):
        a, b = max(c - wd // 2, 0), c + wd // 2; cand[nm] = cbam.count(ch, a, b) * sc * 400 / (b - a)
    cand["genome_average"] = N[s] / GENOME * 400
    lam = max(cand.values()); lim = max(cand, key=cand.get)
    return dict(chip_reads_400bp=k, lambda_local=lam, limiting_background=lim, fold_vs_local=k / lam if lam > 0 else np.nan,
                poisson_p=float(stats.poisson.sf(k - 1, lam)) if lam > 0 else np.nan, **{f"lambda_{n}": v for n, v in cand.items()})
try:
    say("\n=== 2. Local-background test (max of 400 bp / 1 kb / 5 kb / 10 kb input density and genome average; Poisson) ===")
    ref_pk = P[("ChIP_pooled", "q05")] if P[("ChIP_pooled", "q05")] is not None else P[("ChIP_pooled", "q01")]
    rows = []
    for w in WIN.itertuples():
        hit = ref_pk[(ref_pk.chr == w.chr) & (ref_pk.s < w.end) & (ref_pk.e > w.start)]
        if len(hit): c = int(hit.loc[hit.nlq.idxmax(), "pos"]); basis = "summit of strongest pooled peak"
        else: c = (w.start + w.end) // 2; basis = "window centre (no peak)"
        for s in samples: rows.append(dict(gene=w.gene, window=w.window, sample=s, centre=c, centre_basis=basis, **local_test(s, w.chr, c)))
    LB = pd.DataFrame(rows); LB.to_csv(f"{RV}/local_background_test.csv", index=False)
    say(LB[(LB.gene.eq("PNPLA3") | LB.window.eq("promoter")) & (LB["sample"] == "ChIP_pooled")][["gene", "window", "chip_reads_400bp", "lambda_local", "limiting_background", "fold_vs_local", "poisson_p"]].to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    status["2"] = "ok"
except Exception as e: status["2"] = "FAILED"; say("2 failed: " + repr(e)); traceback.print_exc()

# ================================================================ 3. coverage, mappability, blacklist, extra tracks
try:
    say("\n=== 3. Coverage / mappability / blacklist / extra tracks per window ===")
    ctrl_bw = pyBigWig.open(f"{OUT}/tracks/control.RPGC.bw")
    pooled_bw = pyBigWig.open(f"{OUT}/tracks/ChIP_pooled.RPGC.bw") if os.path.exists(f"{OUT}/tracks/ChIP_pooled.RPGC.bw") else None
    mapbw = pyBigWig.open(os.environ["MAP_BW"]) if os.environ.get("MAP_BW") and os.path.exists(os.environ["MAP_BW"]) else None
    bl = None
    if os.environ.get("BLACKLIST_BED") and os.path.exists(os.environ["BLACKLIST_BED"]):
        bl = pd.read_csv(os.environ["BLACKLIST_BED"], sep="\t", header=None, usecols=[0, 1, 2], names=["chr", "s", "e"])
    extra = {}
    for kv in filter(None, os.environ.get("EXTRA_BW", "").split(",")):
        k, v = kv.split("=", 1)
        if os.path.exists(v): extra[k] = pyBigWig.open(v)
    def covfrac(bw, ch, a, b):
        v = np.nan_to_num(np.array(bw.values(ch, a, b), float)); return float((v > 0).mean()) if len(v) else np.nan
    rows = []
    for w in WIN.itertuples():
        r = dict(gene=w.gene, window=w.window, chr=w.chr, start=w.start, end=w.end, width=w.end - w.start,
                 input_bp_covered_fraction=covfrac(ctrl_bw, w.chr, w.start, w.end),
                 pooled_chip_bp_covered_fraction=covfrac(pooled_bw, w.chr, w.start, w.end) if pooled_bw else np.nan,
                 input_reads=cbam.count(w.chr, w.start, w.end), pooled_chip_reads=bam["ChIP_pooled"].count(w.chr, w.start, w.end))
        r["mean_mappability"] = (mapbw.stats(w.chr, w.start, w.end, type="mean")[0] if mapbw else np.nan)
        r["overlaps_blacklist"] = (bool(((bl.chr == w.chr) & (bl.s < w.end) & (bl.e > w.start)).any()) if bl is not None else np.nan)
        for k, bw in extra.items(): r[f"mean_{k}"] = bw.stats(w.chr, w.start, w.end, type="mean")[0]
        rows.append(r)
    CV = pd.DataFrame(rows); CV.to_csv(f"{RV}/window_coverage_mappability.csv", index=False)
    say(CV[CV.gene.eq("PNPLA3")].to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    if mapbw is None: say("   (MAP_BW not given: mappability not computed - 'input_bp_covered_fraction' is a sequencing-coverage proxy only)")
    if bl is None: say("   (BLACKLIST_BED not given: blacklist overlap not checked)")
    status["3"] = "ok"
except Exception as e: status["3"] = "FAILED"; say("3 failed: " + repr(e)); traceback.print_exc()

# ================================================================ 4. orthologs, universe, enrichment, formal interaction tests
def load_orthologs():
    f = os.environ.get("ORTHO_FILE", "")
    if not f or not os.path.exists(f):
        f = os.path.join(os.path.dirname(os.path.abspath(OUT)), "ref", "HOM_MouseHumanSequence.rpt"); os.makedirs(os.path.dirname(f), exist_ok=True)
        if not os.path.exists(f):
            try:
                urllib.request.urlretrieve("https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt", f); say(f"downloaded MGI ortholog report -> {f}")
            except Exception as e: say(f"!! could not download the MGI ortholog report ({e}). Provide ORTHO_FILE (MGI report or CSV 'mouse,human')."); return None
    try:
        d = pd.read_csv(f, sep="\t", dtype=str)
        if "DB Class Key" in d.columns:
            d = d[["DB Class Key", "Common Organism Name", "Symbol"]].dropna(); d["org"] = np.where(d["Common Organism Name"].str.lower().str.startswith("mouse"), "m", np.where(d["Common Organism Name"].str.lower().str.startswith("human"), "h", "x"))
            g = d[d.org != "x"].groupby(["DB Class Key", "org"]).Symbol.apply(list).unstack("org").dropna()
            g = g[(g.m.str.len() == 1) & (g.h.str.len() == 1)]
            return {m[0].upper(): h[0].upper() for m, h in zip(g.m, g.h)}
        d = pd.read_csv(f, dtype=str); cm = next(c for c in d.columns if "mouse" in c.lower()); ch = next(c for c in d.columns if "human" in c.lower())
        d = d[[cm, ch]].dropna().drop_duplicates(); d = d[~d[cm].duplicated(keep=False) & ~d[ch].duplicated(keep=False)]
        return {a.upper(): b.upper() for a, b in zip(d[cm], d[ch])}
    except Exception as e: say(f"!! could not parse ortholog file {f}: {e}"); return None

def or_ci(a, b, c, d):
    t = [x + .5 for x in (a, b, c, d)] if min(a, b, c, d) == 0 else [a, b, c, d]
    lo = np.log(t[0] * t[3] / (t[1] * t[2])); se = np.sqrt(sum(1 / x for x in t)); return np.exp(lo), np.exp(lo - 1.96 * se), np.exp(lo + 1.96 * se), lo, se
def counts(bound, ind):
    bound, ind = np.asarray(bound, bool), np.asarray(ind, bool)
    return int((bound & ind).sum()), int((bound & ~ind).sum()), int((~bound & ind).sum()), int((~bound & ~ind).sum())
def describe_tab(label, bound, ind, base=None):
    a, b, c, d = counts(bound, ind); o, l, h, lo, se = or_ci(a, b, c, d)
    r = dict(dataset=label, n_bound=a + b, bound_induced=a, bound_not_induced=b, unbound_induced=c, unbound_not_induced=d,
             pct_bound_induced=100 * a / max(a + b, 1), pct_unbound_induced=100 * c / max(c + d, 1), OR=o, OR_ci_lo=l, OR_ci_hi=h,
             fisher_p_two_sided=stats.fisher_exact([[a, b], [c, d]])[1], fisher_p_one_sided_greater=stats.fisher_exact([[a, b], [c, d]], alternative="greater")[1],
             adj_OR_expression=np.nan, adj_OR_ci_lo=np.nan, adj_OR_ci_hi=np.nan, adj_p_expression=np.nan)
    if sm is not None and base is not None and np.asarray(bound).sum() > 5:
        try:
            X = sm.add_constant(pd.DataFrame({"bound": np.asarray(bound, float), "logbase": np.log10(np.clip(np.asarray(base, float), 1, None))}))
            m = sm.Logit(np.asarray(ind, float), X).fit(disp=0); ci = np.exp(m.conf_int().loc["bound"].values)
            r.update(adj_OR_expression=float(np.exp(m.params["bound"])), adj_OR_ci_lo=float(ci[0]), adj_OR_ci_hi=float(ci[1]), adj_p_expression=float(m.pvalues["bound"]))
        except Exception as e: say(f"   adjusted logit failed ({label}): {e}")
    return r, lo, se

def interaction(U, bound, label):
    """U: DataFrame with columns gene, ind_liver, ind_hasc, base (liver baseMean). bound: bool array aligned with U."""
    bound = np.asarray(bound, bool); out = dict(definition=label)
    tl = counts(bound, U.ind_liver.values); th = counts(bound, U.ind_hasc.values)
    _, _, _, l1, s1 = or_ci(*tl); _, _, _, l2, s2 = or_ci(*th); d_obs = l2 - l1
    out.update(logOR_liver=l1, logOR_hASC=l2, ratio_of_ORs_hASC_over_liver=float(np.exp(d_obs)), woolf_z=d_obs / np.sqrt(s1 ** 2 + s2 ** 2), woolf_p=float(2 * stats.norm.sf(abs(d_obs / np.sqrt(s1 ** 2 + s2 ** 2)))))
    if smf is not None:
        try:
            L = pd.concat([pd.DataFrame(dict(gene=U.gene.values, ds=0, bound=bound.astype(int), induced=U.ind_liver.astype(int).values)),
                           pd.DataFrame(dict(gene=U.gene.values, ds=1, bound=bound.astype(int), induced=U.ind_hasc.astype(int).values))])
            m = smf.gee("induced ~ bound * ds", "gene", data=L, family=sm.families.Binomial(), cov_struct=sm.cov_struct.Exchangeable()).fit()
            out.update(gee_interaction_logOR=float(m.params["bound:ds"]), gee_interaction_p=float(m.pvalues["bound:ds"]))
        except Exception as e: out.update(gee_interaction_logOR=np.nan, gee_interaction_p=np.nan); say(f"   GEE failed: {e}")
    nb = int(bound.sum()); il, ih = U.ind_liver.values.astype(bool), U.ind_hasc.values.astype(bool); n = len(U)
    def perm(strata=None):
        ds = np.empty(NPERM)
        if strata is None: idx_pool = [(np.arange(n), nb)]
        else: idx_pool = [(np.where(strata == k)[0], int(bound[strata == k].sum())) for k in np.unique(strata)]
        for i in range(NPERM):
            sel = np.concatenate([rng.choice(ix, k, replace=False) for ix, k in idx_pool if k > 0])
            tb = np.zeros(n, bool); tb[sel] = True
            ds[i] = or_ci(*counts(tb, ih))[3] - or_ci(*counts(tb, il))[3]
        return (1 + np.sum(np.abs(ds) >= abs(d_obs))) / (NPERM + 1)
    out["perm_p_difference"] = perm()
    try:
        strata = pd.qcut(U.base.rank(method="first"), 10, labels=False).values; out["perm_p_difference_expression_stratified"] = perm(strata)
    except Exception as e: out["perm_p_difference_expression_stratified"] = np.nan
    return out

try:
    say("\n=== 4. Binding vs response: ortholog-based enrichment, adjusted models and formal between-dataset tests ===")
    orth = load_orthologs()
    if not (need(DE_CSV, "DE_CSV") and need(HASC_CSV, "HASC_CSV")): raise RuntimeError("RNA-seq tables missing")
    liver = pd.read_csv(DE_CSV).sort_values("baseMean", ascending=False).drop_duplicates("GeneName"); liver["mu"] = liver.GeneName.str.upper()
    hasc = pd.read_csv(HASC_CSV).drop_duplicates("Gene"); hasc["hu"] = hasc.Gene.str.upper()
    # ---- TSS table (most upstream TSS per gene) and the bound-gene definitions
    op = gzip.open if GTF.endswith(".gz") else open; tss_rows = []
    with op(GTF, "rt") as f:
        for ln in f:
            if ln[0] == "#": continue
            c = ln.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] != "transcript": continue
            m = re.search(r'gene_name "([^"]+)"', c[8])
            if m: tss_rows.append((m.group(1).upper(), c[0], int(c[3]) - 1 if c[6] == "+" else int(c[4]) - 1, c[6]))
    tss_df = pd.DataFrame(tss_rows, columns=["gene", "chr", "tss", "strand"]).drop_duplicates()
    pool = P[("ChIP_pooled", "q01")]; keep = np.ones(len(pool), bool)
    for r_ in reps:
        x = P[(r_, "q05")]
        if x is None: keep &= False; continue
        k = np.zeros(len(pool), bool)
        for ch in pool.chr.unique():
            xs = x[x.chr == ch].sort_values("s"); sel = (pool.chr == ch).values
            if not len(xs): continue
            S_, E_ = xs.s.values, np.maximum.accumulate(xs.e.values); i = np.searchsorted(S_, pool.e.values[sel], "right") - 1
            k[sel] = (i >= 0) & (E_[np.clip(i, 0, None)] > pool.s.values[sel])
        keep &= k
    hc = pool[keep]; say(f"high-confidence peak set (pooled q<0.01 supported by q<0.05 in every replicate): {len(hc):,} peaks")
    def near(hcp, lo, hi):
        genes = set()
        for ch, g in tss_df.groupby("chr"):
            ps = hcp[hcp.chr == ch].pos.values
            if not len(ps): continue
            o = np.argsort(g.tss.values); t = g.tss.values[o]; names = g.gene.values[o]
            j = np.clip(np.searchsorted(t, ps), 1, len(t) - 1) if len(t) > 1 else np.zeros(len(ps), int)
            cand = np.stack([j - 1, j], 1) if len(t) > 1 else j[:, None]
            dist = np.abs(t[cand] - ps[:, None]); best = cand[np.arange(len(ps)), dist.argmin(1)]; dmin = dist.min(1)
            genes |= set(names[best[(dmin > lo) & (dmin <= hi)]])
        return genes
    def promoter_set(hcp, w):
        genes = set()
        for ch, g in tss_df.groupby("chr"):
            ps = np.sort(hcp[hcp.chr == ch].pos.values)
            if not len(ps): continue
            i0 = np.searchsorted(ps, g.tss.values - w); i1 = np.searchsorted(ps, g.tss.values + w, "right"); genes |= set(g.gene.values[i1 > i0])
        return genes
    # peaks (not summits) within the promoter window, as in step 5, for the primary definition
    prim_f = f"{OD}/bound_genes_pooled_supported_by_both_reps_q05.csv"
    D = {}
    D["A_primary_prespecified: high-confidence peak within TSS+-3kb"] = {g.upper() for g in pd.read_csv(prim_f).gene} if os.path.exists(prim_f) else promoter_set(hc, PROM)
    D["B_pooled_q01 peak within TSS+-3kb"] = {g.upper() for g in pd.read_csv(f"{OD}/bound_genes_pooled_q01.csv").gene} if os.path.exists(f"{OD}/bound_genes_pooled_q01.csv") else promoter_set(pool, PROM)
    cf = f"{OD}/bound_genes_consensus_q01_both_replicates.csv"
    if os.path.exists(cf): D["C_both_replicates_q01 peak within TSS+-3kb"] = {g.upper() for g in pd.read_csv(cf).gene}
    D["D_high-confidence peak within TSS+-1kb"] = promoter_set(hc, 1000)
    dist_only = near(hc, PROM, 50000) - D["A_primary_prespecified: high-confidence peak within TSS+-3kb"]
    D["E_distal_only: high-confidence peak 3-50kb from nearest TSS (nearest-gene assignment, approximate)"] = dist_only
    D["F_promoter_or_distal (A union E)"] = D["A_primary_prespecified: high-confidence peak within TSS+-3kb"] | dist_only
    for k, v in D.items(): say(f"  {k}: {len(v):,} genes")
    # ---- datasets
    rows_sym, rows_orth, rows_int, rows_sens = [], [], [], []
    liver["ind"] = (liver.padj < .05) & (liver.log2FoldChange > 0); hasc["ind"] = (hasc.padj < .05) & (hasc.log2FC_Rosi_vs_DMSO > 0)
    # (i) symbol-matched analysis, as in the first version of panel i
    for k, v in D.items():
        for lab, df, key, base in (("liver TZD (mouse, symbols upper-cased)", liver.dropna(subset=["padj"]), "mu", "baseMean"), ("hASC rosiglitazone (human)", hasc.dropna(subset=["padj"]), "hu", None)):
            r, _, _ = describe_tab(lab, df[key].isin(v).values, df.ind.values, df[base].values if base else None); r.update(definition=k, analysis="symbol-matched"); rows_sym.append(r)
    pd.DataFrame(rows_sym).to_csv(f"{RV}/enrichment_symbol_matched.csv", index=False)
    if orth is None: raise RuntimeError("no ortholog table: ortholog-based analysis skipped (symbol-matched table written)")
    liver["human"] = liver.mu.map(orth); L1 = liver.dropna(subset=["human", "padj"]).drop_duplicates("human")
    H1 = hasc.dropna(subset=["padj"]); universe = sorted(set(L1.human) & set(H1.hu) & set(tss_df.gene))
    U = pd.DataFrame({"gene": universe}).merge(L1[["human", "ind", "baseMean", "log2FoldChange", "stat"]].rename(columns={"human": "gene", "ind": "ind_liver", "baseMean": "base", "stat": "stat_liver"}), on="gene") \
                                         .merge(H1[["hu", "ind", "t_stat"]].rename(columns={"hu": "gene", "ind": "ind_hasc", "t_stat": "stat_hasc"}), on="gene")
    say(f"common universe: {len(U):,} genes with a 1:1 human-mouse ortholog, tested in both RNA-seq datasets and annotated in the HepG2 genome; "
        f"induced: liver {int(U.ind_liver.sum())}, hASC {int(U.ind_hasc.sum())}")
    for k, v in D.items():
        b = U.gene.isin(v).values
        rl, _, _ = describe_tab("liver TZD (mouse; 1:1 orthologs)", b, U.ind_liver.values, U.base.values); rh, _, _ = describe_tab("hASC rosiglitazone (human)", b, U.ind_hasc.values, None)
        for r in (rl, rh): r.update(definition=k, analysis="1:1 orthologs, common universe"); rows_orth.append(r)
        it = interaction(U, b, k); it["n_bound_in_universe"] = int(b.sum()); rows_int.append(it)
        say(f"  {k}\n     liver OR {rl['OR']:.2f} [{rl['OR_ci_lo']:.2f}-{rl['OR_ci_hi']:.2f}] P(2-sided)={rl['fisher_p_two_sided']:.3g}; adj-for-expression OR {rl['adj_OR_expression']:.2f} [{rl['adj_OR_ci_lo']:.2f}-{rl['adj_OR_ci_hi']:.2f}] P={rl['adj_p_expression']:.3g}"
            f"\n     hASC  OR {rh['OR']:.2f} [{rh['OR_ci_lo']:.2f}-{rh['OR_ci_hi']:.2f}] P(2-sided)={rh['fisher_p_two_sided']:.3g}"
            f"\n     difference (hASC/liver OR ratio {it['ratio_of_ORs_hASC_over_liver']:.2f}): Woolf P={it['woolf_p']:.3g}; GEE interaction P={it.get('gee_interaction_p', np.nan):.3g}; permutation P={it['perm_p_difference']:.3g}; expression-stratified permutation P={it['perm_p_difference_expression_stratified']:.3g}")
    U.assign(bound_A=U.gene.isin(D["A_primary_prespecified: high-confidence peak within TSS+-3kb"])).to_csv(f"{RV}/common_universe_table.csv", index=False)
    EO = pd.DataFrame(rows_orth); EO.to_csv(f"{RV}/enrichment_orthologs_common_universe.csv", index=False)
    pd.DataFrame(rows_int).to_csv(f"{RV}/enrichment_between_dataset_tests.csv", index=False)
    EO[EO.definition.str.startswith("A_")].to_csv(f"{RV}/panel_i_statistics.csv", index=False)       # read by step 6
    pd.DataFrame(rows_int)[lambda d: d.definition.str.startswith("A_")].to_csv(f"{RV}/panel_i_interaction.csv", index=False)
    status["4"] = "ok"
except Exception as e: status["4"] = "FAILED/PARTIAL"; say("4: " + repr(e)); traceback.print_exc()

# ================================================================ 5. effect sizes with CI and equivalence for named genes
try:
    say("\n=== 5. Effect sizes (log2 fold change, 95% CI) and equivalence (TOST, 90% CI inside the margin) ===")
    names = sorted(set(W.gene) | {"CD36", "PNPLA2", "PPARG", "PPARGC1A", "SLC2A4", "LPL", "ANGPTL4"}); rows = []
    if need(DE_CSV, "DE_CSV"):
        L = pd.read_csv(DE_CSV); L["u"] = L.GeneName.str.upper()
        RE = pd.read_csv(DE_RECOMP_CSV) if (DE_RECOMP_CSV and os.path.exists(DE_RECOMP_CSV)) else None
        if RE is not None: RE["u"] = RE.GeneName.str.upper()
        else: say("  NOTE: no DESeq2 re-run table found -> liver CIs cannot be computed (the authors' xlsx has no standard error); only authors' log2FC/padj are listed")
        for g in names:
            r = L[L.u == g]
            if r.empty: continue
            r = r.sort_values("baseMean", ascending=False).iloc[0]
            row = {"dataset": "liver TZD vs control (mouse)", "gene": g, "baseMean": r.baseMean, "authors_log2FC": r.log2FoldChange, "authors_padj": r.padj}
            q = RE[RE.u == g].sort_values("baseMean", ascending=False) if RE is not None else pd.DataFrame()
            if len(q) and np.isfinite(q.iloc[0].lfcSE):
                q = q.iloc[0]; b, se = q.log2FoldChange, q.lfcSE
                row.update({"log2FC": b, "ci95_lo": b - 1.96 * se, "ci95_hi": b + 1.96 * se, "fold": 2 ** b, "fold_lo": 2 ** (b - 1.96 * se), "fold_hi": 2 ** (b + 1.96 * se), "padj": q.padj,
                            "equivalent_within_1.5x": bool(b - 1.645 * se > -np.log2(1.5) and b + 1.645 * se < np.log2(1.5)), "equivalent_within_2x": bool(b - 1.645 * se > -1 and b + 1.645 * se < 1),
                            "se_source": "pydeseq2 Wald SE, DESeq2 re-run of the read counts in the same xlsx (re-analysis, not the authors' file)"})
            else:
                row.update({"log2FC": r.log2FoldChange, "padj": r.padj, "se_source": "none (authors' xlsx has no SE)"})
            rows.append(row)
    if need(HASC_CSV, "HASC_CSV"):
        H = pd.read_csv(HASC_CSV); H["u"] = H.Gene.str.upper(); DF = 4          # 5 donors, paired -> t distribution with 4 df
        t975, t95 = stats.t.ppf(0.975, DF), stats.t.ppf(0.95, DF)
        for g in names:
            r = H[H.u == g]
            if r.empty: continue
            r = r.iloc[0]; b = r.log2FC_Rosi_vs_DMSO
            if "log2FC_SE" in H.columns: se, src = r.log2FC_SE, "SD of the 5 donor differences / sqrt(5); t(4 df) intervals"
            else: se, src = (abs(b / r.t_stat) if (r.t_stat and abs(r.t_stat) >= 0.5) else np.nan), "log2FC / t (older table without log2FC_SE)"
            rows.append({"dataset": "hASC rosiglitazone vs DMSO (human)", "gene": g, "baseMean": np.nan, "log2FC": b, "ci95_lo": b - t975 * se, "ci95_hi": b + t975 * se, "fold": 2 ** b, "fold_lo": 2 ** (b - t975 * se), "fold_hi": 2 ** (b + t975 * se),
                         "padj": r.padj, "equivalent_within_1.5x": bool(b - t95 * se > -np.log2(1.5) and b + t95 * se < np.log2(1.5)), "equivalent_within_2x": bool(b - t95 * se > -1 and b + t95 * se < 1), "se_source": src})
    ES = pd.DataFrame(rows); ES.to_csv(f"{RV}/gene_effect_sizes.csv", index=False)
    say(ES[ES.gene.isin(["PNPLA3", "CD36", "FABP4", "PNPLA2", "PCK1", "PLIN1"])][["dataset", "gene", "authors_log2FC", "authors_padj", "log2FC", "ci95_lo", "ci95_hi", "fold", "fold_lo", "fold_hi", "padj", "equivalent_within_1.5x", "equivalent_within_2x", "se_source"]].to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    status["5"] = "ok"
except Exception as e: status["5"] = "FAILED"; say("5 failed: " + repr(e)); traceback.print_exc()

say("\nsection status: " + ", ".join(f"{k}={v}" for k, v in status.items()))
open(f"{RV}/revalidation_report.txt", "w").write("\n".join(report)); print(f"\nDone. Outputs in {RV}/")

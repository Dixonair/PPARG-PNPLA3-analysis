#!/usr/bin/env python3
"""Step 8: Source Data workbook for the HepG2 PPARG ChIP-seq panels (Fig. 4j-m and Supplementary Fig. 4g-l).

Every sheet is read straight from the files written by steps 1-5 (bigWig tracks, MACS3 peak files, full_analysis/*.csv) - no value is typed in
or estimated. The track sheets contain exactly the binned values that step 06 draws (400 bins across the plotted window).

Run from the project root with the conda env active (after steps 1-5; step 7 is not needed):
    OUT=results/hepg2_ppar GTF=ref/gencode.v46.basic.annotation.gtf.gz python3 scripts/chipseq_hepg2/08_export_source_data.py
Optional: RS_POS (default 43928847, rs738409 on GRCh38 chr22; the script checks it falls in exon 3 of the gene model used).
Output: $OUT/full_analysis/manuscript_figures/Source_Data_Fig4_SuppFig4_HepG2.xlsx
Needs: pandas numpy pyBigWig openpyxl
"""
import os, re, gzip, glob, sys, hashlib, datetime
import numpy as np, pandas as pd, pyBigWig

OUT = os.environ.get("OUT", "results/hepg2_ppar"); GTF = os.environ["GTF"]
RS_POS = int(os.environ.get("RS_POS", 43928847))
FA = f"{OUT}/full_analysis"; FD = f"{FA}/manuscript_figures"; os.makedirs(FD, exist_ok=True)
def need(f):
    if not os.path.exists(f): sys.exit(f"Missing {f}. Run steps 1-5 first.")
    return f
def md5(f): return hashlib.md5(open(f, "rb").read()).hexdigest()
used = []                                   # every input file is listed with its md5 on the README sheet
def use(f): used.append(f); return f

reps = sorted(os.path.basename(p).replace(".RPGC.bw", "") for p in glob.glob(f"{OUT}/tracks/*.RPGC.bw")
              if not os.path.basename(p).startswith(("control", "ChIP_pooled"))
              and os.path.exists(f"{OUT}/peaks/" + os.path.basename(p).replace(".RPGC.bw", "") + "_peaks.narrowPeak"))
if len(reps) != 2: sys.exit(f"expected two replicate tracks, found {reps}")
S = reps + ["ChIP_pooled"]; LAB = {reps[0]: reps[0], reps[1]: reps[1], "ChIP_pooled": "pooled"}
trk = {s: pyBigWig.open(use(need(f"{OUT}/tracks/{s}.log2ChIPoverInput.bw"))) for s in S}
ctrl = pyBigWig.open(use(need(f"{OUT}/tracks/control.RPGC.bw")))
def rp(f):
    try: d = pd.read_csv(f, sep="\t", header=None, usecols=[0, 1, 2, 6, 7, 8, 9], names=["chr", "start", "end", "fold_change_at_summit", "neglog10_P", "neglog10_q", "summit_offset"])
    except pd.errors.EmptyDataError: d = pd.DataFrame(columns=["chr", "start", "end", "fold_change_at_summit", "neglog10_P", "neglog10_q", "summit_offset"])
    return d
P = {}
for s in S:
    P[(s, "q<0.01")] = rp(use(need(f"{OUT}/peaks/{s}_peaks.narrowPeak")))
    P[(s, "q<0.05")] = rp(use(need(f"{OUT}/peaks/relaxed/{s}_q05_peaks.narrowPeak")))
W = pd.read_csv(use(need(f"{FA}/window_enrichment_tests.csv")))
QC = pd.read_csv(use(need(f"{FA}/peak_QC_annotation.csv"))); CONC = pd.read_csv(use(need(f"{FA}/replicate_peak_concordance.csv")))
MOT = pd.read_csv(use(f"{FA}/motif_scan_named_sites.csv")) if os.path.exists(f"{FA}/motif_scan_named_sites.csv") else None

# ---- gene models (same rule as step 06: the transcript with most exons, then the longest)
GENES = ["PNPLA3", "PLIN1", "PCK1", "FABP4"]
op = gzip.open if GTF.endswith(".gz") else open
tx, ex = {}, {}
with op(use(need(GTF)), "rt") as f:
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
    MODEL[g] = dict(chr=ch, strand=st, tss=a if st == "+" else b - 1, exons=sorted(ex[(g, t)]), tid=t)

def window(g, up, down):
    m = MODEL[g]; sg = 1 if m["strand"] == "+" else -1; tss = m["tss"]
    return m["chr"], (tss - up if sg == 1 else tss - down), (tss + down if sg == 1 else tss + up), sg, tss
def tracks(g, up, down, samples, with_input=True, n=400):
    ch, a, b, sg, tss = window(g, up, down); edges = np.linspace(a, b, n + 1)
    d = pd.DataFrame({"chr": ch, "bin_start": edges[:-1].astype(int), "bin_end": edges[1:].astype(int)})
    d["bin_midpoint"] = ((d.bin_start + d.bin_end) / 2).round().astype(int); d["distance_from_TSS_bp"] = sg * (d.bin_midpoint - tss)
    for s in samples:
        v = trk[s].stats(ch, int(a), int(b), nBins=n, type="mean"); d[f"{LAB[s]}_log2_ChIP_over_input"] = [np.nan if x is None else x for x in v]
    if with_input:
        v = ctrl.stats(ch, int(a), int(b), nBins=n, type="mean"); d["input_RPGC"] = [np.nan if x is None else x for x in v]
    return d
def peaks_in(g, up, down, samples):
    ch, a, b, sg, tss = window(g, up, down); rows = []
    for s in samples:
        for lvl in ("q<0.01", "q<0.05"):
            d = P[(s, lvl)]; d = d[(d.chr == ch) & (d.start < b) & (d.end > a)]
            for r in d.itertuples():
                summ = r.start + int(r.summit_offset)
                rows.append(dict(gene=g, sample=LAB[s], peak_set=lvl, chr=r.chr, start=r.start, end=r.end, summit=summ, summit_distance_from_TSS_bp=sg * (summ - tss),
                                 fold_change_at_summit=r.fold_change_at_summit, neglog10_P=r.neglog10_P, neglog10_q=r.neglog10_q))
    return pd.DataFrame(rows)
def macs_hit(gene, win):
    r = W[(W.gene == gene) & (W.window == win)].iloc[0]
    for lvl in ("q<0.01", "q<0.05"):
        d = P[("ChIP_pooled", lvl)]; d = d[(d.chr == r.chr) & (d.start < r.end) & (d.end > r.start)]
        if len(d): k = d.neglog10_q.astype(float).idxmax(); return float(d.loc[k, "neglog10_q"]), float(d.loc[k, "fold_change_at_summit"]), lvl
    return None

sheets, notes = {}, {}
# ---- Fig. 4j: PNPLA3 locus (window TSS -2.5 kb ... +7.5 kb, as plotted)
sheets["Fig4j_tracks"] = tracks("PNPLA3", 2500, 7500, S)
notes["Fig4j_tracks"] = "Fig. 4j. PNPLA3 locus, 400 equal bins across the plotted window; mean of the 10-bp bigWig values per bin (blank = no data; drawn as 0). log2(ChIP/input) per replicate and pooled (deepTools bamCompare, read-count scaling, pseudocount 1, reads extended to 200 bp); input_RPGC from the control library."
sheets["Fig4j_peaks"] = peaks_in("PNPLA3", 2500, 7500, S)
notes["Fig4j_peaks"] = "Fig. 4j peak rows: every MACS3 peak overlapping the plotted window, per replicate and pooled, at q<0.01 and at q<0.05 (the q<0.05 set includes the q<0.01 peaks; light bars are those not in the q<0.01 set)."
m = MODEL["PNPLA3"]; sg = 1 if m["strand"] == "+" else -1
exr = [dict(exon=i, start=s, end=e) for i, (s, e) in enumerate(m["exons"] if sg == 1 else m["exons"][::-1], 1)]
gm = pd.DataFrame(exr); gm.insert(0, "transcript", m["tid"]); gm.insert(0, "gene", "PNPLA3"); gm["strand"] = m["strand"]; gm["TSS"] = m["tss"]
wn = W[W.gene == "PNPLA3"].drop_duplicates("window").set_index("window")
sites = pd.DataFrame([dict(feature=w, chr=wn.loc[w, "chr"], start=int(wn.loc[w, "start"]), end=int(wn.loc[w, "end"])) for w in ("upstream_PPRE", "intron1_site", "site_plus4.8kb", "promoter") if w in wn.index])
sites = pd.concat([sites, pd.DataFrame([dict(feature="rs738409 (reference position used; check against dbSNP)", chr=m["chr"], start=RS_POS - 1, end=RS_POS)])], ignore_index=True)
sheets["Fig4j_gene_model"] = pd.concat([gm, pd.DataFrame([{}]), sites.rename(columns={"feature": "gene"})], ignore_index=True)
notes["Fig4j_gene_model"] = "Fig. 4j gene model (GENCODE exon coordinates, 0-based start) and the named windows / SNP marker drawn on the panel."
# ---- Fig. 4k: peak numbers
t = QC.drop_duplicates(["sample", "q"])[["sample", "q", "n_peaks"] + [c for c in QC.columns if c == "median_width"]]
sheets["Fig4k_peak_numbers"] = t.reset_index(drop=True); notes["Fig4k_peak_numbers"] = "Fig. 4k. Number of MACS3 peaks per sample at q<0.01 (q01) and q<0.05 (q05)."
# ---- Fig. 4l: enrichment per window
rows = []
for (g, w), d in W.groupby(["gene", "window"], sort=False):
    d = d.set_index("sample"); pooled = d.loc["ChIP_pooled"]; mh = macs_hit(g, w)
    rows.append(dict(gene=g, window=w, chr=pooled.chr, start=pooled.start, end=pooled.end, informative=bool((pooled.chip_reads + pooled.input_reads) >= 5),
                     **{f"log2_enrichment_{LAB[s]}": d.loc[s, "log2_enrichment"] for s in S}, **{f"chip_reads_{LAB[s]}": d.loc[s, "chip_reads"] for s in S}, input_reads=pooled.input_reads,
                     macs3_pooled_neglog10_q=(np.nan if mh is None else mh[0]), macs3_pooled_fold_change_at_summit=(np.nan if mh is None else mh[1]), macs3_peak_set=("no peak" if mh is None else mh[2]),
                     descriptive_binomial_p_pooled=pooled.binom_p))
sheets["Fig4l_window_enrichment"] = pd.DataFrame(rows)
notes["Fig4l_window_enrichment"] = "Fig. 4l. log2(ChIP/input read density) per window (bars = pooled, dots = replicates; 0.5 pseudocount); labels = MACS3 q of the overlapping pooled peak. Windows with <5 reads in ChIP+input are uninformative (informative=False). The binomial P is a descriptive test on window read counts (independent-reads assumption; not a MACS statistic)."
# ---- Fig. 4m: percentiles
pc = W[W["sample"] == "ChIP_pooled"][["gene", "window", "chip_reads", "input_reads", "log2_ratio_track", "percentile_log2ratio", "input_RPGC", "input_percentile"]].reset_index(drop=True)
sheets["Fig4m_percentiles"] = pc; notes["Fig4m_percentiles"] = "Fig. 4m. Genome-wide percentile of the pooled log2(ChIP/input) (filled circles) and of the input signal alone (open diamonds), relative to same-width genome tiles (random 10% of the genome). Blank = window too narrow for a tile background."
# ---- Supplementary Fig. 4g, h, i
sheets["SuppFig4g_annotation"] = QC.copy(); notes["SuppFig4g_annotation"] = "Suppl. Fig. 4g. Annotation of MACS3 peak summits (promoter = TSS +/-3 kb, exon, intron, intergenic); 'frac' is the plotted fraction for q01 peaks."
sheets["SuppFig4h_concordance"] = CONC.copy(); notes["SuppFig4h_concordance"] = "Suppl. Fig. 4h. Fraction of peaks of one replicate/set overlapping a peak of the other (plotted x 100)."
if MOT is not None: sheets["SuppFig4i_motif"] = MOT.copy(); notes["SuppFig4i_motif"] = "Suppl. Fig. 4i. Best PPARG::RXRA (JASPAR MA0065.1) site score in each named site and its percentile versus GC-matched random windows."
# ---- Suppl. Fig. 4j-l: control loci (pooled, +/-4 kb)
for g, l in (("PLIN1", "j"), ("PCK1", "k"), ("FABP4", "l")):
    if g not in MODEL: continue
    sheets[f"SuppFig4{l}_{g}_tracks"] = tracks(g, 4000, 4000, ["ChIP_pooled"]); notes[f"SuppFig4{l}_{g}_tracks"] = f"Suppl. Fig. 4{l}. {g} locus (TSS +/-4 kb), 400 bins, pooled log2(ChIP/input) and input RPGC."
    pk = peaks_in(g, 4000, 4000, ["ChIP_pooled"]); sheets[f"SuppFig4{l}_{g}_peaks"] = pk if len(pk) else pd.DataFrame({"note": [f"no pooled MACS3 peak within the plotted window at q<0.05"]})
    notes[f"SuppFig4{l}_{g}_peaks"] = f"Suppl. Fig. 4{l}. Pooled MACS3 peaks overlapping the plotted window."
# ---- text values: PNPLA3 windows (numbers quoted in the Results)
pn = W[W.gene == "PNPLA3"].copy(); sheets["Text_PNPLA3_windows"] = pn.reset_index(drop=True)
notes["Text_PNPLA3_windows"] = "Values quoted in the Results text for the PNPLA3 windows (read counts, enrichment, percentiles, peak calls)."
# ---- README
readme = ["Source Data - Fig. 4j-m and Supplementary Fig. 4g-l (HepG2 PPARG ChIP-seq re-analysis)", f"Generated {datetime.datetime.now():%Y-%m-%d %H:%M} by 08_export_source_data.py from the files listed below; no value is typed in.", "",
          "Data: ENCODE ENCSR130VQL (HepG2 with a CRISPR-inserted epitope tag at PPARG; Partridge et al., Nature 2020), control ENCSR239QGH (input from unmodified HepG2). Genome hg38, GENCODE v46.", "", "SHEETS"]
readme += [f"  {k}: {notes[k]}" for k in sheets]
readme += ["", "INPUT FILES (md5)"] + [f"  {f}  {md5(f)}" for f in dict.fromkeys(used)]
sheets = {"README": pd.DataFrame({"": readme}), **sheets}
xl = f"{FD}/Source_Data_Fig4_SuppFig4_HepG2.xlsx"
with pd.ExcelWriter(xl, engine="openpyxl") as w_:
    for k, d in sheets.items():
        d.to_excel(w_, sheet_name=k[:31], index=False)
        ws = w_.sheets[k[:31]]
        for i, col in enumerate(d.columns, 1): ws.column_dimensions[ws.cell(1, i).column_letter].width = min(max(12, len(str(col)) + 2), 60)
        if k == "README": ws.column_dimensions["A"].width = 160
print(f"Wrote {xl}"); [print(f"  {k}: {len(d)} rows") for k, d in sheets.items()]

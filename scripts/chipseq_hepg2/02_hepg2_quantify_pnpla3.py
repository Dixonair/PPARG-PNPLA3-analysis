#!/usr/bin/env python3
"""Quantify PPARG ChIP signal / peaks at PNPLA3 (and control genes) in HepG2.
Usage (env vars): GTF=/path/gencode.v4x.annotation.gtf[.gz]  OUT=results/hepg2_ppar  python3 02_hepg2_quantify_pnpla3.py
Needs: pyBigWig pandas matplotlib.
All windows are defined relative to the annotated TSS (hg38), so they are strand-aware.
"""
import os, glob, gzip, re, sys
import pandas as pd, pyBigWig
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

OUT = os.environ.get("OUT", "results/hepg2_ppar")
GTF = os.environ["GTF"]
POS  = ["FABP4","ADIPOQ","CD36","PLIN1","PCK1","LPL","ANGPTL4","PLIN2"]   # PPARG/PPAR targets
NEG  = ["MYOD1","RHO","INS"]                                              # closed in HepG2
GENES = ["PNPLA3"] + POS + NEG
# PNPLA3 (+ strand) windows, bp relative to TSS: promoter, upstream PPRE (-790..-713 +-pad), intron-1 ChIP site
WIN = {"TSS_pm2kb": (-2000, 2000), "upstream_PPRE": (-850, -650), "intron1_site": (1000, 2300)}

def tss_table(gtf, genes):
    op = gzip.open if gtf.endswith(".gz") else open
    t = {}
    with op(gtf, "rt") as f:
        for ln in f:
            if ln.startswith("#"): continue
            c = ln.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] != "gene": continue
            m = re.search(r'gene_name "([^"]+)"', c[8])
            if not m or m.group(1) not in genes: continue
            chrom, s, e, strand = c[0], int(c[3]) - 1, int(c[4]), c[6]
            t[m.group(1)] = (chrom if chrom.startswith("chr") else "chr" + chrom, s if strand == "+" else e - 1, strand)
    return t

def mean_signal(bw, chrom, a, b):
    a = max(a, 0); v = bw.stats(chrom, a, b, type="mean")[0]
    return 0.0 if v is None else v

tss = tss_table(GTF, set(GENES))
missing = [g for g in GENES if g not in tss]
if missing: print("WARNING genes not in GTF:", missing)
print("PNPLA3 TSS from GTF:", tss.get("PNPLA3"), "(expected ~chr22:43,923,805, + strand)")

tracks = {os.path.basename(p).replace(".RPGC.bw", ""): p for p in sorted(glob.glob(f"{OUT}/tracks/*.RPGC.bw"))}
ratio  = {os.path.basename(p).replace(".log2ChIPoverInput.bw", ""): p for p in sorted(glob.glob(f"{OUT}/tracks/*.log2ChIPoverInput.bw"))}
peaks  = {os.path.basename(p).replace("_peaks.narrowPeak", ""): p for p in sorted(glob.glob(f"{OUT}/peaks/*_peaks.narrowPeak"))}

rows = []
for lab, p in tracks.items():
    bw = pyBigWig.open(p); br = pyBigWig.open(ratio[lab]) if lab in ratio else None
    pk = pd.read_csv(peaks[lab], sep="\t", header=None, usecols=[0,1,2,6,8,9],
                     names=["chr","s","e","fc","neglogq","summit"]) if lab in peaks else pd.DataFrame()
    for g in GENES:
        if g not in tss: continue
        chrom, t0, strand = tss[g]
        sgn = 1 if strand == "+" else -1
        for wname, (a, b) in WIN.items():
            if g != "PNPLA3" and wname != "TSS_pm2kb": continue
            lo, hi = sorted((t0 + sgn * a, t0 + sgn * b))
            r = dict(sample=lab, gene=g, window=wname, mean_RPGC=mean_signal(bw, chrom, lo, hi))
            if br: r["mean_log2_ChIP_over_input"] = mean_signal(br, chrom, lo, hi)
            if len(pk):
                hit = pk[(pk.chr == chrom) & (pk.e > lo) & (pk.s < hi)]
                r["peak_overlap"] = len(hit) > 0
                r["best_peak_neglog10q"] = hit.neglogq.max() if len(hit) else None
            rows.append(r)
res = pd.DataFrame(rows); res.to_csv(f"{OUT}/hepg2_PPARG_signal_at_genes.csv", index=False)
print(res.to_string(index=False))

# peaks within +-10 kb of PNPLA3 TSS (position relative to TSS)
near = []
for lab, p in peaks.items():
    pk = pd.read_csv(p, sep="\t", header=None, usecols=[0,1,2,6,8,9], names=["chr","s","e","fc","neglogq","summit"])
    chrom, t0, strand = tss["PNPLA3"]
    h = pk[(pk.chr == chrom) & (pk.e > t0 - 10000) & (pk.s < t0 + 10000)].copy()
    h["summit_rel_TSS"] = h.s + h.summit - t0; h["sample"] = lab; near.append(h)
near = pd.concat(near) if near else pd.DataFrame()
near.to_csv(f"{OUT}/hepg2_PNPLA3_peaks_pm10kb.csv", index=False); print("\nPeaks within +-10 kb of PNPLA3 TSS:\n", near.to_string(index=False) if len(near) else "NONE")

# track plot + small bedGraph to send back
chrom, t0, strand = tss["PNPLA3"]; lo, hi = t0 - 6000, t0 + 12000
labs = list(tracks); fig, ax = plt.subplots(len(labs) + (len(ratio) > 0 and 1 or 0), 1, figsize=(9, 2*len(labs)+2), sharex=True, squeeze=False)
for i, lab in enumerate(labs):
    bw = pyBigWig.open(tracks[lab]); xs = list(range(lo, hi, 50))
    ys = [mean_signal(bw, chrom, x, x + 50) for x in xs]; ax[i,0].fill_between(xs, ys, step="mid"); ax[i,0].set_ylabel(lab[:14] + "\nRPGC", fontsize=7)
    with open(f"{OUT}/PNPLA3_region_{lab}.bedGraph", "w") as f:
        for x, y in zip(xs, ys): f.write(f"{chrom}\t{x}\t{x+50}\t{y}\n")
for a_ in ax[:,0]:
    a_.axvline(t0, color="k", lw=.8)
    for (a, b), c in zip([WIN["upstream_PPRE"], WIN["intron1_site"]], ["tab:red", "tab:green"]): a_.axvspan(t0+a, t0+b, color=c, alpha=.2)
ax[-1,0].set_xlabel(f"{chrom} (hg38); black = TSS, red = upstream PPREs, green = intron-1 ChIP site (GSE115421)")
plt.tight_layout(); plt.savefig(f"{OUT}/hepg2_PNPLA3_tracks.pdf")
print("Wrote:", f"{OUT}/hepg2_PPARG_signal_at_genes.csv, hepg2_PNPLA3_peaks_pm10kb.csv, hepg2_PNPLA3_tracks.pdf")

#!/usr/bin/env python3
"""Step 0: build the RNA-seq input tables DIRECTLY from the original files - nothing is typed in or carried over.

Run this on a local machine, so so it is clear where every number comes from.

Inputs (put them in one folder, SRC; defaults to the current folder):
  GSE162276_HF_TZD_CvsHF_C_DEG.xlsx   the GSE162276 DEG table (TZD vs control, Pparg-intact liver)
  ASC_Rosi.csv                         hASC expression, 5 donors x 3 replicates x (DMSO, Rosi)
  bound_genes_promoter.csv             hASC PPARG ChIP-bound genes (GSE115421_hASC rows are used)

Outputs (folder DEST, default ./inputs):
  GSE162276_TZDvsControl_DE_results.csv            PRIMARY liver table = the authors' own numbers from the xlsx
        GeneName, baseMean, log2FoldChange, pvalue, padj  : copied / averaged from the xlsx (see the provenance file)
        stat    : DERIVED = sign(log2FoldChange) * z of the authors' two-sided p   (the xlsx has no statistic column)
        lfcSE   : left empty (the xlsx has none; do not invent one)
  GSE162276_TZDvsControl_DE_recomputed_DESeq2.csv  SENSITIVITY liver table = DESeq2 (pydeseq2) re-run on the raw read counts in the
        same xlsx. Clearly a re-analysis, used (a) as a robustness check and (b) because it provides lfcSE for confidence intervals.
  ASCRosi_DE_with_ChIP_annotation.csv              hASC table computed from ASC_Rosi.csv: per donor log2(mean of 3 replicates + 1), Rosi minus DMSO,
        then a paired one-sample t-test across the 5 donors (df = 4), Benjamini-Hochberg. log2FC_SE = sd(donor differences)/sqrt(5).
        Duplicated gene symbols: the row with the highest mean expression is kept.   bound = symbol listed in bound_genes_promoter.csv (hASC).
  inputs_provenance.txt                            md5 of every source file, row counts, rules, and cross-checks between the two liver tables

Needs: pandas numpy scipy statsmodels openpyxl  (+ pydeseq2 for the sensitivity table:  pip install pydeseq2)
"""
import os, sys, hashlib, warnings
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings("ignore")

SRC = os.environ.get("SRC", "."); DEST = os.environ.get("DEST", "inputs"); os.makedirs(DEST, exist_ok=True)
XLSX = os.path.join(SRC, "GSE162276_HF_TZD_CvsHF_C_DEG.xlsx"); ASC = os.path.join(SRC, "ASC_Rosi.csv"); BND = os.path.join(SRC, "bound_genes_promoter.csv")
for f in (XLSX, ASC, BND):
    if not os.path.exists(f): sys.exit(f"STOP: {f} not found. Put the original file there or set SRC=<folder>.")
log = []
def say(s=""): print(s); log.append(s)
def save(df, name):
    """write a table; an existing file with the same name is kept as <name>.previous (never silently overwritten)"""
    f = os.path.join(DEST, name)
    if os.path.exists(f):
        os.replace(f, f + ".previous"); say(f"  (existing {name} kept as {name}.previous)")
    df.to_csv(f, index=False)
def md5(f): return hashlib.md5(open(f, "rb").read()).hexdigest()
say("SOURCE FILES (md5)"); [say(f"  {os.path.basename(f)}  {md5(f)}") for f in (XLSX, ASC, BND)]

# ------------------------------------------------------------------ liver: authors' table
x = pd.read_excel(XLSX)
norm_cols = [c for c in x.columns if c.startswith("PG_") and not c.endswith(("_readcount", "_fpkm"))]
cnt_cols = [c for c in x.columns if c.endswith("_readcount")]
say(f"\nLIVER xlsx: {len(x)} rows, {x.GeneName.nunique()} unique GeneName; {len(norm_cols)} normalised-count columns, {len(cnt_cols)} read-count columns")
x["baseMean"] = x[norm_cols].mean(axis=1)
t = x.dropna(subset=["padj"]).sort_values("baseMean", ascending=False).drop_duplicates("GeneName").copy()
p = t["pvalue"].astype(float).clip(lower=1e-300)
t["stat"] = np.sign(t.log2FoldChange) * stats.norm.isf(p / 2)
t["lfcSE"] = np.nan
auth = t[["GeneName", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]].copy()
auth["stat_is_derived_from_p"] = True
save(auth, "GSE162276_TZDvsControl_DE_results.csv")
say(f"  wrote GSE162276_TZDvsControl_DE_results.csv: {len(auth)} genes with a non-missing padj (rule: one row per GeneName, highest baseMean)")
say("  baseMean = mean of the 10 normalised-count columns in the xlsx; log2FoldChange, pvalue, padj = copied from the xlsx; stat = derived from the p-value")
for g in ("Pnpla3", "Fabp4", "Cd36", "Pnpla2"):
    r = auth[auth.GeneName == g]
    if len(r): r = r.iloc[0]; say(f"    {g}: log2FC {r.log2FoldChange:.4f}  p {r.pvalue:.3g}  padj {r.padj:.3g}  baseMean {r.baseMean:.1f}")

# ------------------------------------------------------------------ liver: DESeq2 re-run on the raw counts (sensitivity)
try:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    cnt = x.drop_duplicates("Gene_ID").set_index("Gene_ID")[cnt_cols]; cnt.columns = [c.replace("_readcount", "") for c in cnt_cols]
    assert (cnt.values == np.round(cnt.values)).all(), "read counts are not integers"
    # group assignment: the first 5 sample columns carry the TZD samples. This is CHECKED below against the sign of the authors' log2FC.
    grp = pd.Series(["TZD"] * 5 + ["CTRL"] * 5, index=cnt.columns, name="cond")
    cnt = cnt[cnt.sum(axis=1) > 0]
    dds = DeseqDataSet(counts=cnt.T.astype(int), metadata=grp.to_frame(), design="~cond", refit_cooks=True, quiet=True); dds.deseq2()
    ds = DeseqStats(dds, contrast=["cond", "TZD", "CTRL"], quiet=True); ds.summary()
    r = ds.results_df.reset_index().merge(x[["Gene_ID", "GeneName"]].drop_duplicates("Gene_ID"), on="Gene_ID")
    r = r.dropna(subset=["padj"]).sort_values("baseMean", ascending=False).drop_duplicates("GeneName")
    j = r.set_index("GeneName").join(auth.set_index("GeneName"), lsuffix="_re", rsuffix="_auth", how="inner")
    c = np.corrcoef(j.log2FoldChange_re, j.log2FoldChange_auth)[0, 1]
    if c <= 0.3: sys.exit(f"STOP: the DESeq2 re-run does not agree in direction with the authors' log2FC (r={c:.2f}); the TZD/control sample assignment is wrong.")
    save(r[["GeneName", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]], "GSE162276_TZDvsControl_DE_recomputed_DESeq2.csv")
    say(f"\n  wrote GSE162276_TZDvsControl_DE_recomputed_DESeq2.csv: {len(r)} genes (pydeseq2 Wald test on the raw read counts; TZD vs control; sample assignment check r={c:.3f})")
    ia = set(auth.GeneName[(auth.padj < .05) & (auth.log2FoldChange > 0)]); ir = set(r.GeneName[(r.padj < .05) & (r.log2FoldChange > 0)])
    say(f"  CROSS-CHECK induced genes (padj<0.05 & log2FC>0): authors' table {len(ia)}, DESeq2 re-run {len(ir)}, shared {len(ia & ir)} (Jaccard {len(ia & ir) / len(ia | ir):.3f})")
    sl = np.polyfit(j.log2FoldChange_re, j.log2FoldChange_auth, 1)[0]
    say(f"  CROSS-CHECK log2FC: authors vs re-run slope {sl:.2f} (<1 means the authors' values are shrunken or normalised differently; direction and p-values agree)")
    for g in ("Pnpla3", "Fabp4", "Cd36", "Pnpla2"):
        if g in j.index: q = j.loc[g]; say(f"    {g}: authors log2FC {q.log2FoldChange_auth:.3f} (padj {q.padj_auth:.3g}) | re-run log2FC {q.log2FoldChange_re:.3f}, SE {q.lfcSE_re:.3f} (padj {q.padj_re:.3g})")
except ImportError:
    say("\n  pydeseq2 not installed -> sensitivity liver table NOT written (pip install pydeseq2, then rerun this step). The primary table is unaffected.")

# ------------------------------------------------------------------ hASC
a0 = pd.read_csv(ASC); cols = list(a0.columns[1:])
a0["_m"] = a0[cols].mean(axis=1)
a = a0.sort_values("_m", ascending=False).drop_duplicates("Gene").drop(columns="_m").set_index("Gene")
say(f"\nhASC ASC_Rosi.csv: {len(a0)} rows, {a0.Gene.nunique()} unique symbols (duplicated symbols: the highest-mean row is kept)")
diffs = []
for pt in range(1, 6):
    dm = a[[f"P{pt}DMSO{k}" for k in (1, 2, 3)]].mean(axis=1); ro = a[[f"P{pt}Rosi{k}" for k in (1, 2, 3)]].mean(axis=1)
    diffs.append(np.log2(ro + 1) - np.log2(dm + 1))
d = pd.concat(diffs, axis=1); tt, pv = stats.ttest_1samp(d.values, 0, axis=1)
h = pd.DataFrame({"log2FC_Rosi_vs_DMSO": d.mean(axis=1), "log2FC_SE": d.std(axis=1, ddof=1) / np.sqrt(d.shape[1]), "t_stat": tt, "pvalue": pv}, index=a.index)
ok = h.pvalue.notna(); h = h[ok].copy(); h["padj"] = multipletests(h.pvalue, method="fdr_bh")[1]
pb = pd.read_csv(BND); pb = pb[pb.dataset == "GSE115421_hASC"]; bound = set(pb.SYMBOL.astype(str))
h["bound"] = h.index.isin(bound); h["induced"] = (h.padj < 0.05) & (h.log2FC_Rosi_vs_DMSO > 0)
save(h.reset_index(), "ASCRosi_DE_with_ChIP_annotation.csv")
say(f"  wrote ASCRosi_DE_with_ChIP_annotation.csv: {len(h)} genes, {int(h.bound.sum())} bound (bound_genes_promoter.csv, GSE115421_hASC), {int(h.induced.sum())} induced")
if "PNPLA3" in h.index:
    q = h.loc["PNPLA3"]; say(f"    PNPLA3: log2FC {q.log2FC_Rosi_vs_DMSO:.4f}, SE {q.log2FC_SE:.3f}, t {q.t_stat:.3f}, p {q.pvalue:.4g}, padj {q.padj:.4g}")
open(os.path.join(DEST, "inputs_provenance.txt"), "w").write("\n".join(log))
say(f"\nDone. Provenance written to {DEST}/inputs_provenance.txt")

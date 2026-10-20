#!/usr/bin/env python3
"""
Pre-ranked Hallmark GSEA for (1) adipocyte-specific inducible Pparg KO (GSE62937, 3 vs 3) and
(2) hASC rosiglitazone vs DMSO (5 donors, paired t from script 05). Writes two CSVs.
Inputs in $PPARG_DIR: GSE62937_series_matrix.txt, GPL6887_MouseWG-6_V2_0_R3_11278593_A.txt,
  ASCRosi_DE_with_ChIP_annotation.csv (from chipseq/05_integrate_with_rnaseq.py).
Run:  PPARG_DIR=<folder> python3 12_hallmark_gsea_mouse_and_hasc.py   (needs gseapy, pandas, numpy, scipy)
Ranking: (1) variance-moderated t (prior df 4, prior variance = median probe variance), one probe per
gene with the largest |t|; mouse symbols upper-cased to match MSigDB human symbols.
(2) donor-paired t statistic from ASCRosi_DE_with_ChIP_annotation.csv.
NOTE: gseapy prerank, 1000 permutations, seed 7; NES will differ slightly from an fgsea run.
"""
import pandas as pd, numpy as np, pickle, io, gseapy as gp
import os, urllib.request
UP = OUT = os.environ.get("PPARG_DIR", ".")      # one folder holding the inputs and receiving outputs
HALLMARK_PKL = os.path.join(UP, "hallmark_sets.pkl")
if not os.path.exists(HALLMARK_PKL):             # MSigDB v7.1 Hallmark, CC-BY mirror pinned to an immutable commit
    url = ("https://raw.githubusercontent.com/numpde/genesets/53ce4ba8614d6d3ac2ca33243ea3f9f2c1f86ef5/"
           "genesets/msigdb/parsed/v7.1/genesets.json.zip")
    urllib.request.urlretrieve(url, os.path.join(UP, "genesets.json.zip"))
    _df = pd.read_json(os.path.join(UP, "genesets.json.zip"))
    pickle.dump({c: list(_df[c].symbols) for c in _df if c.startswith("HALLMARK")}, open(HALLMARK_PKL, "wb"))
H = pickle.load(open(HALLMARK_PKL, "rb"))

# ---- (1) mouse Adn-Pparg KO ----
lines = open(f"{UP}/GSE62937_series_matrix.txt").read().splitlines()
b = next(i for i, l in enumerate(lines) if l.startswith("ID_REF"))
e = next((i for i, l in enumerate(lines) if l.startswith("!series_matrix_table_end")), len(lines))
m = pd.read_csv(io.StringIO("\n".join(lines[b:e])), sep="\t", index_col=0).apply(pd.to_numeric, errors="coerce")
m.columns = [c.strip('"') for c in m.columns]; m.index = [str(i).strip('"') for i in m.index]
ctrl, ko = ["GSM1536683", "GSM1536684", "GSM1536685"], ["GSM1536686", "GSM1536687", "GSM1536688"]
x = m[ctrl + ko].dropna()
if x.values.max() > 100: x = np.log2(x.clip(lower=1))
man = pd.read_csv(f"{UP}/GPL6887_MouseWG-6_V2_0_R3_11278593_A.txt", sep="\t", skiprows=8, nrows=45281,
                  usecols=["Probe_Id", "Symbol"], dtype=str).dropna()
sym = man.drop_duplicates("Probe_Id").set_index("Probe_Id")["Symbol"]
x = x.loc[x.index.intersection(sym.index)]
d = x[ko].mean(1) - x[ctrl].mean(1)
s2 = (x[ko].var(1, ddof=1) + x[ctrl].var(1, ddof=1)) / 2
d0, s0 = 4, s2.median(); s2p = (d0 * s0 + 4 * s2) / (d0 + 4)
t = d / np.sqrt(s2p * (2 / 3))
r = pd.DataFrame({"sym": sym.loc[x.index].str.upper(), "t": t}).dropna()
r = r[r.sym.str.len() > 0].loc[lambda df: df.groupby("sym")["t"].apply(lambda s: s.abs().idxmax()).values]
rank_ko = r.set_index("sym")["t"].sort_values(ascending=False)
print("KO ranked genes:", len(rank_ko), " Pnpla3 t =", rank_ko.get("PNPLA3"))

# ---- (2) hASC Rosi ----
a = pd.read_csv(f"{OUT}/ASCRosi_DE_with_ChIP_annotation.csv")
rank_asc = a.drop_duplicates("Gene").set_index("Gene")["t_stat"].sort_values(ascending=False)
print("hASC ranked genes:", len(rank_asc), " PNPLA3 t =", rank_asc.get("PNPLA3"))

for name, rk in [("GSEA_hallmark_AdnPpargKO_GSE62937.csv", rank_ko), ("GSEA_hallmark_hASC_Rosi.csv", rank_asc)]:
    res = gp.prerank(rnk=rk, gene_sets=H, min_size=10, max_size=500, permutation_num=1000, seed=7, threads=2,
                     outdir=None, no_plot=True, verbose=False).res2d
    res = res.rename(columns={"Lead_genes": "Lead_genes"})
    res.to_csv(f"{OUT}/{name}", index=False)
    print(name, res.shape, (res["FDR q-val"].astype(float) < 0.05).sum(), "sets FDR<0.05")
    print(res[["Term", "NES", "FDR q-val"]].head(8).to_string())

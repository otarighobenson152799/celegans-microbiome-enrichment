import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import pearsonr
from scipy.stats import false_discovery_control

PROC = Path("/mnt/d/ce_microbiome/data_compost/processed")
clr = pd.read_csv(PROC / "clr_matrix.tsv", sep="\t", index_col=0)
tax = pd.read_csv(PROC / "taxa_classified.tsv", sep="\t")
seq2genus = dict(zip(tax["sequence"], tax["Genus"]))

genus_map = pd.Series(clr.index.map(lambda s: seq2genus.get(s, "Unclassified")), index=clr.index)
clr_genus = clr.groupby(genus_map).mean()

targets = ["Brucella", "Ochrobactrum", "Escherichia-Shigella", "Lysinibacillus",
           "Rummeliibacillus", "Membranicola", "Sporomusa", "Streptococcus",
           "Leadbetterella", "Parapedobacter"]

present = [t for t in targets if t in clr_genus.index]
X = clr_genus.loc[present]
n = len(present)

corr = np.zeros((n, n))
pvals = np.zeros((n, n))
for i in range(n):
    for j in range(i, n):
        if i == j:
            corr[i,j] = 1.0; pvals[i,j] = 0.0
        else:
            r, p = pearsonr(X.iloc[i], X.iloc[j])
            corr[i,j] = r; corr[j,i] = r
            pvals[i,j] = p; pvals[j,i] = p

off = pvals[np.triu_indices(n, k=1)]
q = false_discovery_control(off)
qmat = np.zeros((n, n))
qmat[np.triu_indices(n, k=1)] = q
qmat = qmat + qmat.T

print("\n=== CO-OCCURRENCE WITH BRUCELLA / OCHROBACTRIUM (q<0.05, |r|>0.2) ===")
edge_count = 0
for idx_target in [i for i, g in enumerate(present) if g in ("Brucella", "Ochrobactrum")]:
    g_name = present[idx_target]
    print(f"\n  {g_name}:")
    for j in range(n):
        if j == idx_target: continue
        if qmat[idx_target, j] < 0.05 and abs(corr[idx_target, j]) > 0.2:
            print(f"    {present[j]:22s}  r = {corr[idx_target, j]:+.3f}   q = {qmat[idx_target, j]:.3g}")
            edge_count += 1

print(f"\nTotal significant edges involving Brucella/Ochrobactrum: {edge_count}")

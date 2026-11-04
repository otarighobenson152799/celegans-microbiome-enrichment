import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import pearsonr, spearmanr
from scipy.cluster.hierarchy import dendrogram
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

# Paths
COMPOST = Path("/mnt/d/ce_microbiome/data_compost/processed")
KEGG = Path("/mnt/d/ce_microbiome/results_compost/revision_kegg_expanded")
ASR = Path("/mnt/d/ce_microbiome/results_compost/revision_asr")
OUT = Path("/mnt/d/ce_microbiome/figures_final/Figure_8_Master.png")

# Setup 2x3 grid
fig = plt.figure(figsize=(22, 14))
gs = fig.add_gridspec(2, 3, height_ratios=[1.2, 1], hspace=0.3, wspace=0.25)

# ==========================================
# ROW 1: MECHANISM
# ==========================================
# --- A: Network ---
axA = fig.add_subplot(gs[0, 0])
clr = pd.read_csv(COMPOST / "clr_matrix.tsv", sep="\t", index_col=0)
tax = pd.read_csv(COMPOST / "taxa_classified.tsv", sep="\t")
seq2g = dict(zip(tax["sequence"], tax["Genus"]))
g_map = pd.Series(clr.index.map(lambda s: seq2g.get(s, "Unclassified")), index=clr.index)
cg = clr.groupby(g_map).mean()
tgts = ["Brucella", "Ochrobactrum", "Escherichia-Shigella", "Lysinibacillus", "Rummeliibacillus", "Membranicola", "Sporomusa", "Streptococcus", "Leadbetterella", "Parapedobacter"]
pres = [t for t in tgts if t in cg.index]
X = cg.loc[pres]
n = len(pres)
corr = np.zeros((n, n)); pval = np.zeros((n, n))
for i in range(n):
    for j in range(i, n):
        if i == j: corr[i,j] = 1.0
        else:
            r, p = pearsonr(X.iloc[i], X.iloc[j])
            corr[i,j] = corr[j,i] = r; pval[i,j] = pval[j,i] = p
ang = np.linspace(0, 2*np.pi, n, endpoint=False)
pos = {i: (np.cos(ang[i]), np.sin(ang[i])) for i in range(n)}
for i in range(n):
    for j in range(i+1, n):
        if pval[i,j] < 0.01 and abs(corr[i,j]) > 0.2:
            axA.plot([pos[i][0], pos[j][0]], [pos[i][1], pos[j][1]], color="#d62728" if corr[i,j]>0 else "#1f77b4", lw=1+3*abs(corr[i,j]), alpha=0.6)
for i, g in enumerate(pres):
    is_f = g in ("Brucella", "Ochrobactrum")
    axA.scatter(*pos[i], s=900 if is_f else 600, color="#FFD700" if is_f else "#4A4A8A", edgecolor="black", zorder=5)
    axA.text(pos[i][0]*1.15, pos[i][1]*1.15, g, ha="center", va="center", fontsize=9, fontweight="bold" if is_f else "normal", bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="gray", alpha=0.9))
axA.set_xlim(-1.4, 1.4); axA.set_ylim(-1.4, 1.4); axA.axis("off")
axA.set_title("A. Co-occurrence Network", fontweight="bold", fontsize=14, loc='left')

# --- B: ASR ---
axB = fig.add_subplot(gs[0, 1])
Z = np.load(ASR / "asr_linkage.npy")
s = pd.read_csv(KEGG / "m00115_survey_categorized.tsv", sep="\t")
s = s[s.clean_category != "exclude_outgroup"].drop_duplicates(subset=["org_code"]).reset_index(drop=True)
cds = s.org_code.tolist()
cat = dict(zip(s.org_code, s.clean_category))
m = dict(zip(s.org_code, s.has_M00115.astype(str).str.upper()=="YES"))
lbls = [f"{c} {'[M00115+]' if m[c] else '[M00115-]'}" for c in cds]
dendrogram(Z, labels=lbls, orientation="right", ax=axB, leaf_font_size=6, color_threshold=0)
for lbl in axB.get_yticklabels():
    code = lbl.get_text().split()[0]
    if cat[code] == "environmental_named": lbl.set_color("#2ca02c")
    elif cat[code].startswith("outgroup") or cat[code] == "uncharacterized_brucella_sp": lbl.set_color("#9467bd")
    else: lbl.set_color("#d62728")
axB.set_xlabel("Jaccard distance", fontsize=9)
axB.set_title("B. Ancestral State Reconstruction", fontweight="bold", fontsize=14, loc='left')

# --- C: Guild Distribution ---
axC = fig.add_subplot(gs[0, 2])
others = ["Escherichia", "Lysinibacillus", "Leadbetterella", "Sporomusa", "Rummeliibacillus", "Streptococcus", "Shigella"]
pct = [100, 100, 100, 33, 0, 0, 0]
axC.bar(["Brucella/Ochro\n(Env)"], [100], color="#2ca02c", edgecolor="black", width=0.6)
axC.bar(others, pct, color="#4A4A8A", edgecolor="black", width=0.6, alpha=0.7)
axC.set_ylabel("% Genomes carrying M00115", fontweight="bold")
axC.set_ylim(0, 115)
axC.axhline(100, color="gray", linestyle="--", alpha=0.5)
axC.tick_params(axis='x', rotation=45)
axC.set_title("C. Guild Functional Distribution", fontweight="bold", fontsize=14, loc='left')

# ==========================================
# ROW 2: DIVERSITY (Only A and B)
# ==========================================
# --- D: Alpha Diversity ---
axD = fig.add_subplot(gs[1, 0])
# Reconstruct metadata for Shannon
clr_d = pd.read_csv(COMPOST / "clr_matrix.tsv", sep="\t", index_col=0)
rel = np.exp(clr_d); rel = rel.div(rel.sum(axis=0), axis=1)
shannon = -np.sum(rel * np.log(rel + 1e-9), axis=0)
samples = list(clr_d.columns)
meta = pd.DataFrame({
    "env": [s.split("_")[0] for s in samples],
    "type": ["Worm" if "_W" in s else "Soil" for s in samples],
    "shannon": shannon.values,
})
sns.boxplot(data=meta, x='env', y='shannon', hue='type', order=['AP','BP','OR','POT'], palette={'Worm': '#d62728', 'Soil': '#1f77b4'}, ax=axD)
axD.set_xlabel("Environment", fontweight='bold')
axD.set_ylabel("Shannon Index", fontweight='bold')
axD.set_title("D. Alpha Diversity (Shannon)", fontweight="bold", fontsize=14, loc='left')
axD.legend(title="", loc='upper right')

# --- E: Beta Diversity PCoA ---
axE = fig.add_subplot(gs[1, 1])
X_pca = clr_d.T.values
Xc = X_pca - X_pca.mean(axis=0, keepdims=True)
U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
pcs = U[:, :2] * S[:2]
var_exp = (S[:2] ** 2) / (S ** 2).sum() * 100
types = ["Worm" if "_W" in s else "Soil" for s in samples]
colors = {"Worm": "#d62728", "Soil": "#1f77b4"}
markers = {"Worm": "o", "Soil": "^"}
for t in ["Worm", "Soil"]:
    mask = np.array([tt == t for tt in types])
    axE.scatter(pcs[mask, 0], pcs[mask, 1], c=colors[t], marker=markers[t], s=130, edgecolor="black", alpha=0.85, label=t)
axE.set_xlabel(f"PC1 ({var_exp[0]:.1f}%)", fontweight='bold')
axE.set_ylabel(f"PC2 ({var_exp[1]:.1f}%)", fontweight='bold')
axE.set_title("E. Beta Diversity (PCoA of CLR)", fontweight="bold", fontsize=14, loc='left')
axE.legend(title="")

# Hide the last panel (bottom right)
ax_empty = fig.add_subplot(gs[1, 2])
ax_empty.axis("off")

plt.tight_layout()
fig.savefig(OUT, dpi=300, bbox_inches="tight", facecolor="white")
print(f"✅ Master merged figure saved to: {OUT}")

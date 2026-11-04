import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# --- Paths ---
BASE = Path("/mnt/d/ce_microbiome")
OUT_DIR = BASE / "Manuscript final supplementary figure"
OUT_DIR.mkdir(exist_ok=True, parents=True)
OUT_FILE = OUT_DIR / "Supplementary_Figure_1.png"

# --- Dataset Mapping ---
datasets = [
    ("data_ceMbio", "CeMbio"),
    ("data_bigbiome", "BIGbiome"),
    ("data_aging", "Aging"),
    ("data_dirksen", "Dirksen"),
    ("data_infection", "Environ"),
    ("data_compost", "Compost")
]

richness_data = []
phylum_data = []

print("Extracting data for Supplementary Figure 1...")

for folder, label in datasets:
    proc_dir = BASE / folder / "processed"
    if not proc_dir.exists():
        print(f"  ⚠️ Skipping {label}: folder {proc_dir} not found")
        continue
        
    clr_file = proc_dir / "clr_matrix.tsv"
    if not clr_file.exists():
        print(f"  ⚠️ Skipping {label}: clr_matrix.tsv not found")
        continue
        
    try:
        clr = pd.read_csv(clr_file, sep="\t", index_col=0)
        
        # 1. Per-sample richness (ASVs > 0)
        sample_richness = (np.exp(clr) > 0).sum(axis=0)
        for val in sample_richness.values:
            richness_data.append({"Dataset": label, "Richness": val})
            
        # 2. Phylum composition
        tax_file = proc_dir / "taxa_classified.tsv"
        if not tax_file.exists():
            tax_file = proc_dir / "taxa.tsv"
            
        if tax_file.exists():
            tax = pd.read_csv(tax_file, sep="\t")
            
            # Dynamically find the Phylum column
            phy_col = None
            for c in ["Phylum", "phylum", "Taxonomy"]:
                if c in tax.columns:
                    phy_col = c
                    break
            
            if phy_col:
                if phy_col == "Taxonomy":
                    # If raw string, extract phylum (e.g., "d__Bacteria;p__Proteobacteria...")
                    def extract_phy(t):
                        for part in str(t).split(";"):
                            if part.strip().startswith("p__"):
                                return part.strip()[3:]
                        return "Unclassified"
                    phylum_series = tax["Taxonomy"].apply(extract_phy)
                else:
                    phylum_series = tax[phy_col]
                
                # Match sequence index
                seq_col = "sequence" if "sequence" in tax.columns else tax.index
                seq2phy = dict(zip(seq_col, phylum_series))
                phyla = [seq2phy.get(idx, "Unclassified") for idx in clr.index]
                
                # Aggregate relative abundance
                rel = np.exp(clr); rel = rel.div(rel.sum(axis=0), axis=1)
                rel.index = phyla
                phylum_rel = rel.groupby(level=0).mean().mean(axis=1)
                
                # Top 10 + Other
                top_phyla = phylum_rel.nlargest(10)
                other_sum = phylum_rel[~phylum_rel.index.isin(top_phyla.index)].sum()
                
                for phy, ab in top_phyla.items():
                    phylum_data.append({"Dataset": label, "Phylum": phy, "Abundance": ab * 100})
                if other_sum > 0:
                    phylum_data.append({"Dataset": label, "Phylum": "Other", "Abundance": other_sum * 100})
                    
    except Exception as e:
        print(f"  ❌ Error processing {label}: {e}")

# --- Plotting ---
df_rich = pd.DataFrame(richness_data)
df_phy = pd.DataFrame(phylum_data)

fig, axes = plt.subplots(1, 2, figsize=(18, 7), gridspec_kw={'width_ratios': [1, 1.2]})

# Panel A: Richness Boxplot
axA = axes[0]
order = [d[1] for d in datasets if d[1] in df_rich["Dataset"].unique()]
if not df_rich.empty:
    sns.boxplot(data=df_rich, x="Dataset", y="Richness", order=order, ax=axA, palette="Set2", width=0.6, showfliers=False)
    sns.stripplot(data=df_rich, x="Dataset", y="Richness", order=order, ax=axA, color="black", alpha=0.3, size=3, jitter=True)
axA.set_title("A. Per-sample ASV richness across datasets", fontweight="bold", loc="left", fontsize=16)
axA.set_ylabel("Observed ASVs per sample", fontsize=14)
axA.set_xlabel("")
axA.tick_params(axis='x', rotation=30, labelsize=12)
axA.tick_params(axis='y', labelsize=12)
axA.spines['top'].set_visible(False)
axA.spines['right'].set_visible(False)

# Panel B: Phylum Composition
axB = axes[1]
if not df_phy.empty:
    pivot_phy = df_phy.pivot(index="Dataset", columns="Phylum", values="Abundance").fillna(0)
    pivot_phy = pivot_phy.reindex([d for d in order if d in pivot_phy.index])
    pivot_phy.plot(kind="bar", stacked=True, ax=axB, colormap="tab20", edgecolor="black", linewidth=0.5, width=0.7)
    axB.set_title("B. Mean phylum-level composition per dataset", fontweight="bold", loc="left", fontsize=16)
    axB.set_ylabel("Mean Relative Abundance (%)", fontsize=14)
    axB.set_xlabel("")
    axB.tick_params(axis='x', rotation=30, labelsize=12)
    axB.tick_params(axis='y', labelsize=12)
    axB.legend(title="Phylum", bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=11)
    axB.spines['top'].set_visible(False)
    axB.spines['right'].set_visible(False)
else:
    axB.text(0.5, 0.5, "Taxonomy data not available", ha='center', va='center', fontsize=14)

plt.tight_layout()
plt.savefig(OUT_FILE, dpi=300, bbox_inches="tight", facecolor="white")
print(f"\n✅ SUCCESS! Supplementary Figure 1 saved to:\n{OUT_FILE}")

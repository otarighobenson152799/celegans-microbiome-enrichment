import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# Data
data = {
    "Dataset": ["CeMbio", "BIGbiome", "Aging", "Dirksen", "Environ", "Compost"],
    "Accession": ["PRJEB37101/35", "PRJNA540192", "PRJNA982115", "ERP014569", "PRJNA1002096", "PRJNA856419"],
    "Host / design": [
        "N2, 12-strain defined", "43 wild genotypes, Day1/Day3", 
        "N2, 8 timepoints (L4→D8)", "3 genotypes, worm+lawn", 
        "High-diversity natural-like", "4 envs × worm/soil × 3"
    ],
    "Samples": [146, 274, 31, 109, 114, 24],
    "ASVs": [424, 44, 4485, 143, 11503, 5201],
    "Role in paper": [
        "Controlled baseline", "Genotype variation", "Longitudinal",
        "Experimental ground truth", "High-complexity context", "Natural environments (novel)"
    ]
}
df = pd.DataFrame(data)

fig = plt.figure(figsize=(20, 14))
gs = GridSpec(3, 2, figure=fig, height_ratios=[1, 1.2, 1.2], hspace=0.4, wspace=0.2)

# ==========================================
# Panel A: Table
# ==========================================
axA = fig.add_subplot(gs[0, :])
axA.axis('off')
axA.set_title("A. Six independent C. elegans microbiome datasets analyzed", loc="left", fontweight="bold", fontsize=20, pad=20)

# Create the table
table_data = df.copy()
table_data["Samples"] = table_data["Samples"].astype(str)
table_data["ASVs"] = table_data["ASVs"].apply(lambda x: f"{x:,}")

table = axA.table(cellText=table_data.values, colLabels=table_data.columns, loc='center', cellLoc='center')
table.auto_set_font_size(False)
table.set_fontsize(14)
table.scale(1, 2) # Stretch cells vertically

# Style the table
for (row, col), cell in table.get_celld().items():
    cell.set_edgecolor('black')
    if row == 0:
        cell.set_facecolor('#333333')
        cell.set_text_props(color='white', fontweight='bold')
    elif row == 6: # Highlight Compost
        cell.set_facecolor('#e6f2ff')
    else:
        cell.set_facecolor('white')

# ==========================================
# Panel B: Sample Size
# ==========================================
axB = fig.add_subplot(gs[1, 0])
colors = ['#d3d3d3'] * 5 + ['#e7298a'] # Highlight Compost
bars = axB.bar(df["Dataset"], df["Samples"], color=colors, edgecolor='black', linewidth=1.5)
for bar in bars:
    yval = bar.get_height()
    axB.text(bar.get_x() + bar.get_width()/2, yval + 5, int(yval), ha='center', va='bottom', fontweight='bold', fontsize=14)

axB.set_title("B. Sample size", loc="left", fontweight="bold", fontsize=18, pad=15)
axB.set_ylabel("Samples", fontsize=16)
axB.set_ylim(0, 320)
axB.tick_params(axis='x', rotation=30, labelsize=14)
axB.tick_params(axis='y', labelsize=14)
axB.spines['top'].set_visible(False)
axB.spines['right'].set_visible(False)

# ==========================================
# Panel C: Community Richness
# ==========================================
axC = fig.add_subplot(gs[1, 1])
bars = axC.bar(df["Dataset"], df["ASVs"], color=colors, edgecolor='black', linewidth=1.5)
axC.set_yscale("log")
for bar in bars:
    yval = bar.get_height()
    axC.text(bar.get_x() + bar.get_width()/2, yval * 1.15, f"{int(yval):,}", ha='center', va='bottom', fontweight='bold', fontsize=14)

axC.set_title("C. Community richness", loc="left", fontweight="bold", fontsize=18, pad=15)
axC.set_ylabel("ASVs (log scale)", fontsize=16)
axC.set_ylim(20, 25000)
axC.tick_params(axis='x', rotation=30, labelsize=14)
axC.tick_params(axis='y', labelsize=14)
axC.spines['top'].set_visible(False)
axC.spines['right'].set_visible(False)

# ==========================================
# Panel D: Design Space
# ==========================================
axD = fig.add_subplot(gs[2, :])

# Define coordinates (X: Controlled -> Natural, Y: Single host -> Multi-host)
points = {
    "CeMbio": (1, 1),
    "Aging": (4, 1),
    "Dirksen": (2.5, 3),
    "Environ": (8, 1),
    "Compost": (8.5, 4.5),
    "BIGbiome": (3, 8)
}

for name, (x, y) in points.items():
    is_compost = (name == "Compost")
    color = '#e7298a' if is_compost else '#555555'
    size = 600 if is_compost else 400
    axD.scatter(x, y, s=size, color=color, edgecolor='black', linewidth=1.5, zorder=5)
    
    # Position labels dynamically to avoid overlap
    if name == "Compost":
        axD.text(x, y + 0.6, name, ha='center', fontweight='bold', fontsize=14)
    elif name == "BIGbiome":
        axD.text(x, y + 0.6, name, ha='center', fontweight='bold', fontsize=14)
    elif name == "Environ":
        axD.text(x, y + 0.6, name, ha='center', fontweight='bold', fontsize=14)
    else:
        axD.text(x, y - 0.8, name, ha='center', fontweight='bold', fontsize=14)

axD.set_xlim(-1, 11)
axD.set_ylim(-1, 10)

# Add dividing dashed lines
axD.axhline(2.5, color='gray', linestyle='--', alpha=0.6)
axD.axvline(5.5, color='gray', linestyle='--', alpha=0.6)

# Titles and labels
axD.set_title("D. Design space: how each dataset tests host-microbiome specificity", loc="left", fontweight="bold", fontsize=18, pad=15)
axD.set_xlabel("Controlled / defined  →  Natural / high-diversity", fontsize=16, fontweight="bold", labelpad=15)
axD.set_ylabel("Single host  →  Multi-host", fontsize=16, fontweight="bold", labelpad=15)

# Remove ticks
axD.set_xticks([])
axD.set_yticks([])

# Add annotation text cleanly at the bottom
axD.text(0.5, 0.05, "* Compost = only dataset spanning natural environments", 
         transform=axD.transAxes, ha='center', fontsize=16, fontweight="bold", color="#e7298a")

# Thick border for panel D
for spine in axD.spines.values():
    spine.set_linewidth(2)

plt.tight_layout()
out_path = "/mnt/d/ce_microbiome/figures_final/Figure_1_Perfect.png"
fig.savefig(out_path, dpi=300, bbox_inches='tight', facecolor='white')
print(f"✅ Perfect Figure 1 saved to: {out_path}")

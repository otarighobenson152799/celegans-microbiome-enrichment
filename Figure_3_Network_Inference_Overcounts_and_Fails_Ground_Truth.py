import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import matplotlib.patches as mpatches

BASE = Path("/mnt/d/ce_microbiome")
OUT = BASE / "figures_final" / "Figure_3_Complete.png"

print("Generating the complete 8-panel Figure 3...")

# Setup 4x2 grid
fig = plt.figure(figsize=(22, 26))
gs = fig.add_gridspec(4, 2, hspace=0.4, wspace=0.25)

# ==========================================
# Panel A: Overcounting gradient
# ==========================================
axA = fig.add_subplot(gs[0, 0])
datasets = ['BIGbiome', 'CeMbio', 'Aging', 'Environ']
labels = ['BIGbiome\n(43 geno)', 'CeMbio\n(defined)', 'Aging\n(N2)', 'Environ\n(high-div)']
values = [1.3, 3.0, 5.0, 12.0]
colors = ['#6baed6', '#31a354', '#e6ab02', '#d95f02']

bars = axA.bar(labels, values, color=colors, edgecolor='black', width=0.6)
for bar in bars:
    yval = bar.get_height()
    axA.text(bar.get_x() + bar.get_width()/2, yval + 0.5, f"{yval:.1f}x", ha='center', fontweight='bold', fontsize=14)

axA.set_ylabel("Overcounting (sig ASVs / genera)", fontsize=14, fontweight="bold")
axA.set_title("A. Overcounting gradient: 1.3x to 12x", loc="left", fontweight="bold", fontsize=18)
axA.set_ylim(0, 14)
axA.spines['top'].set_visible(False); axA.spines['right'].set_visible(False)

# ==========================================
# Panel B: Network fails Dirksen ground truth
# ==========================================
axB = fig.add_subplot(gs[0, 1])
taxa = ["Rhodococcus", "Chryseobacterium", "Stenotrophomonas", "Achromobacter", "Brucella", 
        "Corynebacterium", "Rothia", "Granulicatella", "Actinomyces", "Streptococcus", "Kingella", "Capnocytophaga"]
q_vals = [1.0, 1.0, 1.0, 1.0, 1.0, 0.04, 0.03, 0.03, 0.02, 0.02, 0.01, 0.01]
colors_B = ["#2ca02c"]*5 + ["#d62728"]*7

axB.barh(taxa, q_vals, color=colors_B, edgecolor='black')
axB.axvline(0.05, color='red', linestyle='--', linewidth=2)
axB.invert_yaxis()
axB.set_xlabel("BH q-value", fontsize=14, fontweight="bold")
axB.set_title("B. Network fails Dirksen ground truth", loc="left", fontweight="bold", fontsize=18)
axB.spines['top'].set_visible(False); axB.spines['right'].set_visible(False)

# ==========================================
# Panel C: Dirksen sample classification
# ==========================================
axC = fig.add_subplot(gs[1, 0])
classes = ["Worm\n(>10%)", "Intermediate\n(1-10%)", "Lawn\n(<1%)"]
counts = [102, 5, 1]
colors_C = ['#2ca02c', '#f0ad4e', '#6c757d']

bars = axC.bar(classes, counts, color=colors_C, edgecolor='black', width=0.6)
for bar in bars:
    yval = bar.get_height()
    axC.text(bar.get_x() + bar.get_width()/2, yval + 2, int(yval), ha='center', fontweight='bold', fontsize=14)

axC.set_ylabel("# samples", fontsize=14, fontweight="bold")
axC.set_title("C. Dirksen sample classification (108 worm)", loc="left", fontweight="bold", fontsize=18)
axC.set_ylim(0, 115)
axC.spines['top'].set_visible(False); axC.spines['right'].set_visible(False)

# ==========================================
# Panel D: Overcounting scales with richness
# ==========================================
axD = fig.add_subplot(gs[1, 1])
richness = [44, 424, 4485, 11503]
overcounting = [1.3, 3.0, 5.0, 12.0]
names = ["BIGbiome", "CeMbio", "Aging", "Environ"]

axD.scatter(richness, overcounting, s=300, color='#d62728', edgecolor='black', zorder=5)
for i, name in enumerate(names):
    axD.text(richness[i] * 1.2, overcounting[i], name, fontsize=12, fontweight='bold', va='center')

axD.set_xscale("log")
axD.set_xlabel("Community richness (ASVs, log)", fontsize=14, fontweight="bold")
axD.set_ylabel("Overcounting", fontsize=14, fontweight="bold")
axD.set_title("D. Overcounting scales with richness", loc="left", fontweight="bold", fontsize=18)
axD.spines['top'].set_visible(False); axD.spines['right'].set_visible(False)

# ==========================================
# Panel E: Mean abundance across datasets
# ==========================================
axE = fig.add_subplot(gs[2, 0])
datasets_E = ["Aging", "CeMbio", "Compost", "Dirksen", "Environ", "BIGbiome"]
abundance = [0.0, 0.17, 4.34, 6.09, 9.37, 24.35]
colors_E = ['#1f77b4', '#d62728', '#e377c2', '#ff7f0e', '#8c564b', '#1f77b4']

bars = axE.bar(datasets_E, abundance, color=colors_E, edgecolor='black', width=0.6)
for bar in bars:
    yval = bar.get_height()
    text = f"{yval:.2f}%" if yval > 0 else "ND"
    axE.text(bar.get_x() + bar.get_width()/2, yval + 0.5, text, ha='center', fontweight='bold', fontsize=12)

axE.set_ylabel("Mean abundance (%)", fontsize=14, fontweight="bold")
axE.set_title("E. Mean abundance of Brucella/Ochrobactrum", loc="left", fontweight="bold", fontsize=18)
axE.set_ylim(0, 28)
axE.spines['top'].set_visible(False); axE.spines['right'].set_visible(False)

# ==========================================
# Panel F: Log2FC worm vs soil in Compost
# ==========================================
axF = fig.add_subplot(gs[2, 1])
envs = ["AP", "BP", "OR", "POT"]
brucella_lfc = [8.0, 7.2, 7.3, 15.4]
ochro_lfc = [3.6, 4.5, 3.5, 2.6]
x = np.arange(len(envs))
width = 0.35

axF.bar(x - width/2, brucella_lfc, width, label='Brucella', color='#8b0000', edgecolor='black')
axF.bar(x + width/2, ochro_lfc, width, label='Ochrobactrum', color='#a0522d', edgecolor='black')
axF.axhline(0, color='gray', linestyle='--')

axF.set_xticks(x)
axF.set_xticklabels(envs, fontweight='bold')
axF.set_ylabel("log2FC (worm/soil)", fontsize=14, fontweight="bold")
axF.set_title("F. Log2FC enrichment in Compost", loc="left", fontweight="bold", fontsize=18)
axF.legend(frameon=False)
axF.spines['top'].set_visible(False); axF.spines['right'].set_visible(False)

# ==========================================
# Panel G: Brucella abundance in worm vs lawn
# ==========================================
axG = fig.add_subplot(gs[3, 0])
types = ["Worm\n(n=102)", "Lawn\n(n=1)"]
abund = [6.36, 0.08]
colors_G = ['#2ca02c', '#8c564b']

bars = axG.bar(types, abund, color=colors_G, edgecolor='black', width=0.5)
for bar in bars:
    yval = bar.get_height()
    axG.text(bar.get_x() + bar.get_width()/2, yval + 0.2, f"{yval:.2f}%", ha='center', fontweight='bold', fontsize=14)

axG.set_ylabel("Brucella abundance (%)", fontsize=14, fontweight="bold")
axG.set_title("G. Brucella abundance in Dirksen", loc="left", fontweight="bold", fontsize=18)
axG.set_ylim(0, 7.5)
axG.spines['top'].set_visible(False); axG.spines['right'].set_visible(False)

# ==========================================
# Panel H: Detection vs enrichment matrix (FIXED)
# ==========================================
axH = fig.add_subplot(gs[3, 1])
axH.axis('off')
axH.set_title("H. Detection vs enrichment evidence matrix", loc="left", fontweight="bold", fontsize=18, pad=15)

# Using \n to wrap text so it fits inside the cells
matrix_data = [
    ["Dataset", "Brucella", "Ochrobactrum"],
    ["Compost", "Enriched\n(worm vs soil)", "Enriched\n(worm vs soil)"],
    ["Dirksen", "Enriched\n(worm vs lawn)", "Enriched\n(worm vs lawn)"],
    ["BIGbiome", "Detected +\ngenotype-assoc.", "Detected +\ngenotype-assoc."],
    ["Environ", "Detected only", "Detected only"],
    ["CeMbio", "Detected only", "Detected only"],
    ["Aging", "Absent (0.00%)", "Absent (0.00%)"]
]

table = axH.table(cellText=matrix_data, loc='center', cellLoc='center')
table.auto_set_font_size(False)
table.set_fontsize(12) 
table.scale(1, 3) # Increased height to accommodate the \n line breaks

for (row, col), cell in table.get_celld().items():
    cell.set_edgecolor('black')
    
    # Set custom column widths (Dataset=0.25, Brucella=0.375, Ochrobactrum=0.375)
    if col == 0:
        cell.set_width(0.25)
    else:
        cell.set_width(0.375)

    if row == 0:
        cell.set_facecolor('#333333')
        cell.set_text_props(color='white', fontweight='bold')
    elif col == 0:
        cell.set_facecolor('#f2f2f2')
        cell.set_text_props(fontweight='bold')
    else:
        text = cell.get_text().get_text()
        if "Enriched" in text:
            cell.set_facecolor('#2ca02c')
            cell.set_text_props(color='white', fontweight='bold')
        elif "Detected +" in text:
            cell.set_facecolor('#1f77b4')
            cell.set_text_props(color='white', fontweight='bold')
        elif "Detected only" in text:
            cell.set_facecolor('#f0ad4e')
        elif "Absent" in text:
            cell.set_facecolor('#d62728')
            cell.set_text_props(color='white', fontweight='bold')

plt.tight_layout()
plt.savefig(OUT, dpi=300, bbox_inches="tight", facecolor="white")
print(f"✅ SUCCESS! Complete 8-panel Figure 3 saved to:\n{OUT}")

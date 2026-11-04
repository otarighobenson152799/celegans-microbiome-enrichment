import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

BASE = Path("/mnt/d/ce_microbiome")
SUPP_DIR = BASE / "Manuscript final supplementary figure"

print("Generating Perfect Supplementary Figure 5...")

fig, ax = plt.subplots(figsize=(14, 8))

modules = ["M00023\n(Tryptophan)", "M00115\n(NAD+)", "M00027\n(GABA)", "M00022\n(Shikimate)"]
env_pct = [100, 100, 11, 100]
path_pct = [95, 0, 63, 100]

x = np.arange(len(modules))
width = 0.35

ax.bar(x - width/2, env_pct, width, label="Environmental (n=9)", color="#2ca02c", edgecolor="black")
ax.bar(x + width/2, path_pct, width, label="Pathogenic (n=43)", color="#d62728", edgecolor="black")

ax.set_xticks(x)
ax.set_xticklabels(modules, fontweight="bold", fontsize=14)
ax.set_ylabel("% of genomes carrying pathway", fontweight="bold", fontsize=16)

# THE FIX: Set the title with a massive pad so it sits at the very top
ax.set_title("Targeted survey of neuroactive precursor pathways", loc="left", fontweight="bold", fontsize=20, pad=50)

# THE FIX: Use fig.legend() to place it relative to the figure, not the axes
# This places it between the title and the plot area
fig.legend(frameon=False, fontsize=14, loc="upper center", bbox_to_anchor=(0.5, 0.92), ncol=2)

# Give the chart breathing room
ax.set_ylim(0, 140)

# Add bar values
for i in range(len(modules)):
    if env_pct[i] > 0: ax.text(i - width/2, env_pct[i] + 3, f"{env_pct[i]}%", ha='center', fontweight='bold', fontsize=12)
    if path_pct[i] > 0: ax.text(i + width/2, path_pct[i] + 3, f"{path_pct[i]}%", ha='center', fontweight='bold', fontsize=12)

# Clean annotation for GABA
ax.annotate("Pathogen-enriched\n(opposite direction)", 
            xy=(2 + width/2, 85),      
            xytext=(2 + width/2, 120), 
            ha='center', fontsize=12, style='italic', color='#d62728', fontweight='bold',
            arrowprops=dict(arrowstyle="->", color="#d62728", lw=2))

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Manually adjust layout so the title, legend, and plot don't overlap
fig.subplots_adjust(top=0.82, bottom=0.1, left=0.1, right=0.95)

# Save with bbox_inches="standard" instead of tight to respect our manual layout
fig.savefig(SUPP_DIR / "Supplementary_Figure_5.png", dpi=300, facecolor="white")
plt.close()
print("✅ Supplementary Figure 5 is now perfectly clean and non-overlapping.")

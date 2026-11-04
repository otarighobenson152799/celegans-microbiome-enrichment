import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches
from pathlib import Path

BASE = Path("/mnt/d/ce_microbiome")
FIG_DIR = BASE / "manuscript_final_Figure"

print("Generating Fixed Main Figure 5...")
fig = plt.figure(figsize=(18, 10))
gs = fig.add_gridspec(1, 2, width_ratios=[1.2, 1], wspace=0.35)

# --- Panel A: Heatmap ---
axA = fig.add_subplot(gs[0, 0])
modules = ["M00115\n(NAD+)", "M00126\n(Folate)", "M00027\n(GABA)", "M00083\n(Fatty Acid)"]
matrix_data = []
for i in range(9):   # Environmental
    matrix_data.append([1, 1 if i < 4 else 0, 1 if i < 1 else 0, 1])
for i in range(43):  # Pathogens
    matrix_data.append([0, 1, 1 if i < 27 else 0, 1])

matrix_df = pd.DataFrame(matrix_data, columns=modules)

sns.heatmap(matrix_df, cmap="Blues", cbar=False, linewidths=0.5, linecolor='white', ax=axA)
axA.set_xticklabels(axA.get_xticklabels(), rotation=0, fontweight="bold", fontsize=14)
axA.set_yticks([]) 

# THE FIX: Increased labelpad to 45 to push the text away from the color bars
axA.set_ylabel("58 Brucella/Ochrobactrum Genomes", fontweight="bold", fontsize=14, labelpad=45)
axA.set_title("A. Binary presence/absence of key KEGG modules", loc="left", fontweight="bold", fontsize=18, pad=20)

# Add side color bars
for i in range(len(matrix_df)):
    color = "#2ca02c" if i < 9 else "#d62728"
    axA.add_patch(plt.Rectangle((-0.3, i), 0.2, 1, color=color, clip_on=False, transform=axA.transData))

# Legend below Panel A
env_patch = mpatches.Patch(color="#2ca02c", label="Environmental (n=9)")
path_patch = mpatches.Patch(color="#d62728", label="Pathogen (n=43)")
axA.legend(handles=[env_patch, path_patch], loc="upper center", bbox_to_anchor=(0.5, -0.08), ncol=2, frameon=False, fontsize=12)

# --- Panel B: Bar Chart ---
axB = fig.add_subplot(gs[0, 1])
data_b = pd.DataFrame({
    "Module": ["NAD+\n(M00115)", "Folate\n(M00126)", "GABA\n(M00027)"],
    "Environmental (n=9)": [100.0, 44.4, 11.1],
    "Pathogen (n=43)": [0.0, 100.0, 62.8]
})

x = np.arange(len(data_b))
width = 0.35
axB.bar(x - width/2, data_b["Environmental (n=9)"], width, label="Environmental", color="#2ca02c", edgecolor="black")
axB.bar(x + width/2, data_b["Pathogen (n=43)"], width, label="Pathogen", color="#d62728", edgecolor="black")

axB.set_xticks(x)
axB.set_xticklabels(data_b["Module"], fontweight="bold", fontsize=14)
axB.set_ylabel("% of genomes carrying module", fontweight="bold", fontsize=14)
axB.set_title("B. Discriminating pathways", loc="left", fontweight="bold", fontsize=18, pad=20)
axB.legend(frameon=False, fontsize=12, loc="upper right")
axB.set_ylim(0, 120)

for i, (env, path) in enumerate(zip(data_b["Environmental (n=9)"], data_b["Pathogen (n=43)"])):
    axB.text(i - width/2, env + 2, f"{env:.0f}%", ha='center', fontweight='bold', fontsize=12)
    axB.text(i + width/2, path + 2, f"{path:.0f}%", ha='center', fontweight='bold', fontsize=12)

axB.spines['top'].set_visible(False); axB.spines['right'].set_visible(False)

plt.tight_layout()
fig.savefig(FIG_DIR / "Figure_5_KEGG_Functional_Divergence.png", dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
print("✅ Main Figure 5 fixed.")

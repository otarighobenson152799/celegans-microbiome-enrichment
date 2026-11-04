import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from PIL import Image, ImageDraw, ImageFont
from pathlib import Path

BASE = Path("/mnt/d/ce_microbiome")
OUT_DIR = BASE / "Manuscript final supplementary figure"
OUT_DIR.mkdir(exist_ok=True, parents=True)
OUT_FILE = OUT_DIR / "Supplementary_Figure_6.png"

TEMP_A = Path("/tmp/supp6_A_permanova.png")
TEMP_B = Path("/tmp/supp6_B_permdisp.png")
SRC_C = BASE / "figures_final/Supp_Fig_11/Supp_Fig_11.png"

print("Generating Supplementary Figure 6...")

# ==========================================
# Generate Panel A: PERMANOVA
# ==========================================
print("  Generating Panel A (PERMANOVA)...")
perm_data = pd.DataFrame({
    "Term": ["Environment", "Worm/Soil", "Interaction"],
    "R2": [0.212, 0.136, 0.194],
    "p_value": ["p = 0.001", "p = 0.001", "p = 0.001"]
})

fig, ax = plt.subplots(figsize=(10, 6))
sns.barplot(data=perm_data, x="Term", y="R2", ax=ax, palette="Set1", edgecolor="black")
ax.set_ylabel("R² (Variance Explained)", fontweight="bold", fontsize=14)
ax.set_xlabel("Model Term", fontweight="bold", fontsize=14)
ax.set_title("A. PERMANOVA by terms", loc="left", fontweight="bold", fontsize=18, pad=15)

for i, row in perm_data.iterrows():
    ax.text(i, row["R2"] + 0.01, row["p_value"], ha="center", fontsize=12, fontweight="bold")
    
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_ylim(0, 0.25)
plt.tight_layout()
plt.savefig(TEMP_A, dpi=300, bbox_inches="tight", facecolor="white")
plt.close()

# ==========================================
# Generate Panel B: PERMDISP
# ==========================================
print("  Generating Panel B (PERMDISP)...")
clr_file = BASE / "data_compost/processed/clr_matrix.tsv"
clr = pd.read_csv(clr_file, sep="\t", index_col=0)
X = clr.T.values
types = ["Worm" if "_W" in s else "Soil" for s in clr.columns]

dist_centroids = []
for t in ["Worm", "Soil"]:
    mask = np.array([tt == t for tt in types])
    centroid = X[mask].mean(axis=0)
    dists = np.sqrt(((X[mask] - centroid) ** 2).sum(axis=1))
    dist_centroids.extend(dists)
    
df_disp = pd.DataFrame({"Type": types, "Distance": dist_centroids})

fig, ax = plt.subplots(figsize=(10, 6))
sns.boxplot(data=df_disp, x="Type", y="Distance", ax=ax, palette={"Worm": "#d62728", "Soil": "#1f77b4"}, width=0.5)
sns.stripplot(data=df_disp, x="Type", y="Distance", ax=ax, color="black", alpha=0.5, jitter=True)
ax.set_ylabel("Distance to Centroid", fontweight="bold", fontsize=14)
ax.set_xlabel("")
ax.set_title("B. PERMDISP (Distances to Centroids)", loc="left", fontweight="bold", fontsize=18, pad=15)
ax.text(0.5, 0.95, "PERMDISP F = 1.44, p = 0.263 (n.s.)", transform=ax.transAxes, ha="center", fontweight="bold", 
        bbox=dict(facecolor='white', edgecolor='gray', alpha=0.9, boxstyle='round,pad=0.3'))
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig(TEMP_B, dpi=300, bbox_inches="tight", facecolor="white")
plt.close()

# ==========================================
# Stitch A, B, and C vertically
# ==========================================
print("  Stitching final 3-panel figure...")

imgs = []
for p in [TEMP_A, TEMP_B, SRC_C]:
    if p.exists():
        imgs.append(Image.open(p).convert("RGB"))
    else:
        print(f"⚠️ Missing {p}, using blank.")
        imgs.append(Image.new("RGB", (1000, 800), "white"))

TARGET_WIDTH = 2000
PADDING = 80

resized = []
for img in imgs:
    w, h = img.size
    new_h = int(h * (TARGET_WIDTH / w))
    resized.append(img.resize((TARGET_WIDTH, new_h), Image.Resampling.LANCZOS))

total_h = sum(i.height for i in resized) + PADDING * (len(resized) + 1)
canvas = Image.new("RGB", (TARGET_WIDTH + 2 * PADDING, total_h), "white")
draw = ImageDraw.Draw(canvas)

try:
    FONT = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 100)
except IOError:
    FONT = ImageFont.load_default()

labels = ["A", "B", "C"]
y_offset = PADDING
for i, (img, label) in enumerate(zip(resized, labels)):
    canvas.paste(img, (PADDING, y_offset))
    # Skip drawing 'C' label if the source image already has it internally
    if label != "C":
        draw.text((PADDING, y_offset - 110), label, fill="black", font=FONT)
    y_offset += img.height + PADDING

canvas.save(OUT_FILE, dpi=(300, 300))
print(f"\n✅ SUCCESS! Supplementary Figure 6 saved to:\n{OUT_FILE}")

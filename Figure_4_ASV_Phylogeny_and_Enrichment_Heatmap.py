import os
from PIL import Image, ImageDraw, ImageFont
from pathlib import Path

SRC = Path("/mnt/d/ce_microbiome/figures_final/Supp_Fig_9/Supp_Fig_9.png")
DEST = Path("/mnt/d/ce_microbiome/Manuscript final supplementary figure")
DEST.mkdir(exist_ok=True, parents=True)
OUT = DEST / "Supplementary_Figure_4.png"

print("Generating correctly labeled Supplementary Figure 4...")

if not SRC.exists():
    print(f"❌ Missing {SRC}. Cannot generate figure.")
    exit(1)

img = Image.open(SRC).convert("RGB")
w, h = img.size

# Create a slightly taller canvas to put the A and B labels at the top
PADDING_TOP = 150
new_h = h + PADDING_TOP
canvas = Image.new("RGB", (w, new_h), "white")
canvas.paste(img, (0, PADDING_TOP))

draw = ImageDraw.Draw(canvas)

try:
    FONT = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 120)
except IOError:
    FONT = ImageFont.load_default()

# The tree is on the left (roughly 0-30% of width)
# The heatmap is on the right (roughly 30-100% of width)
# Place "A" over the tree, "B" over the heatmap
draw.text((50, 10), "A", fill="black", font=FONT)
draw.text((int(w * 0.32), 10), "B", fill="black", font=FONT)

canvas.save(OUT, dpi=(300, 300))
print(f"\n✅ SUCCESS! Properly labeled Supplementary Figure 4 saved to:\n{OUT}")

import os
from PIL import Image, ImageDraw, ImageFont
from pathlib import Path

SRC = Path("/mnt/d/ce_microbiome/figures_final")
DEST = Path("/mnt/d/ce_microbiome/Manuscript final supplementary figure")
DEST.mkdir(exist_ok=True, parents=True)
OUT_FILE = DEST / "Supplementary_Figure_3.png"

try:
    FONT = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 100)
except IOError:
    FONT = ImageFont.load_default()

print("Generating clean 2-panel Supplementary Figure 3...")

# Removed the redundant Dirksen panel. Keeping only Environ and Compost.
img_paths = [
    SRC / "Supp_Fig_2/Panel_A.png",  # A: Environ Overcounting
    SRC / "Supp_Fig_3/Panel_A.png"   # B: Compost ASV Validation
]
labels = ["A", "B"]

TARGET_WIDTH = 2200
PADDING = 100

# Load and resize
imgs = []
for p, label in zip(img_paths, labels):
    if p.exists():
        img = Image.open(p).convert("RGB")
        w, h = img.size
        new_h = int(h * (TARGET_WIDTH / w))
        imgs.append(img.resize((TARGET_WIDTH, new_h), Image.Resampling.LANCZOS))
    else:
        print(f"⚠️ Missing {p.name}. Check the filename!")
        imgs.append(Image.new("RGB", (TARGET_WIDTH, 800), "white"))

# Stack vertically
total_h = sum(i.height for i in imgs) + PADDING * (len(imgs) + 1)
canvas = Image.new("RGB", (TARGET_WIDTH + 2 * PADDING, total_h), "white")
draw = ImageDraw.Draw(canvas)

y_offset = PADDING
for i, (img, label) in enumerate(zip(imgs, labels)):
    canvas.paste(img, (PADDING, y_offset))
    # Draw the A, B label
    draw.text((PADDING, y_offset - 110), label, fill="black", font=FONT)
    y_offset += img.height + PADDING

canvas.save(OUT_FILE, dpi=(300, 300))
print(f"\n✅ SUCCESS! Clean 2-panel Supplementary Figure 3 saved to:\n{OUT_FILE}")

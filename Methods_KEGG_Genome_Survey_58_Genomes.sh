#!/bin/bash
set -uo pipefail

OUT="/mnt/d/ce_microbiome/results_compost/revision_kegg_expanded"
mkdir -p "$OUT"
cd "$OUT"

echo "=========================================="
echo "Stage 2A v2: Expand KEGG M00115 survey"
echo "(using /find/genomes/* endpoints)"
echo "=========================================="

echo "[1/5] Query Brucella genomes..."
curl -s --max-time 30 "https://rest.kegg.jp/find/genomes/brucella" > brucella_raw.txt
echo "    Brucella genomes: $(wc -l < brucella_raw.txt)"

echo ""
echo "[2/5] Query Ochrobactrum genomes..."
curl -s --max-time 30 "https://rest.kegg.jp/find/genomes/ochrobactrum" > ochrobactrum_raw.txt
echo "    Ochrobactrum genomes: $(wc -l < ochrobactrum_raw.txt)"

echo ""
echo "[3/5] Also query 'anthropi' (current KEGG reclassification)..."
curl -s --max-time 30 "https://rest.kegg.jp/find/genomes/anthropi" > anthropi_raw.txt
echo "    'anthropi' genomes: $(wc -l < anthropi_raw.txt)"

echo ""
echo "[4/5] Also query 'Rhizobiaceae' family (broader context)..."
curl -s --max-time 30 "https://rest.kegg.jp/find/genomes/rhizobiaceae" > rhizobiaceae_raw.txt
echo "    Rhizobiaceae genomes: $(wc -l < rhizobiaceae_raw.txt)"

echo ""
echo "[5/5] Building combined organism list..."

# Parse all KEGG result files into (code, name) pairs
# Format from /find: "T00072  bme; Brucella melitensis bv. 1 16M"
python3 << 'PYEOF'
import re
from pathlib import Path

OUT = Path("/mnt/d/ce_microbiome/results_compost/revision_kegg_expanded")
sources = ["brucella_raw.txt", "ochrobactrum_raw.txt",
           "anthropi_raw.txt", "rhizobiaceae_raw.txt"]

organisms = {}  # code -> (name, source_category)
for src in sources:
    path = OUT / src
    if not path.exists() or path.stat().st_size == 0:
        continue
    for line in path.read_text().splitlines():
        m = re.match(r"^T\d+\s+(\w+);\s+(.+)$", line.strip())
        if not m:
            continue
        code = m.group(1)
        name = m.group(2)
        if code in organisms:
            continue
        # Categorize
        if "Brucella" in name and any(x in name for x in ["abortus","melitensis","suis","ovis","canis","microti","pinnipedialis","ceti","neotomae","papionis","vulpis","inopinata"]):
            cat = "Brucella_classical_pathogen"
        elif "anthropi" in name.lower():
            cat = "Brucella_anthropi_environmental"
        elif "Ochrobactrum" in name:
            cat = "Ochrobactrum_environmental"
        elif "Brucella" in name:
            cat = "Brucella_other"
        else:
            cat = "Rhizobiaceae_other"
        organisms[code] = (name, cat, src)

print(f"Unique organisms to query: {len(organisms)}")
for code, (name, cat, src) in sorted(organisms.items()):
    print(f"  {code:6s}  [{cat:35s}]  {name}")

# Save the list for the bash loop
with open(OUT / "organism_list.tsv", "w") as f:
    for code, (name, cat, src) in sorted(organisms.items()):
        f.write(f"{code}\t{name}\t{cat}\n")
print(f"\nSaved organism_list.tsv")
PYEOF

echo ""
echo "=========================================="
echo "Querying M00115 for each organism..."
echo "=========================================="

# Write header
printf "org_code\torganism\tcategory\thas_M00115\tn_genes_M00115\n" > m00115_survey.tsv

TOTAL=$(wc -l < organism_list.tsv)
i=0
while IFS=$'\t' read -r code name category; do
    i=$((i+1))
    clean_name=$(echo "$name" | sed 's/\t/ /g; s/"//g')
    
    # Query KEGG module links for this organism
    n=$(curl -s --max-time 20 "https://rest.kegg.jp/link/module/${code}" 2>/dev/null | grep -c "M00115" || true)
    n=${n:-0}
    
    if [ "$n" -gt 0 ]; then has="YES"; else has="no"; fi
    
    printf "%s\t%s\t%s\t%s\t%s\n" "$code" "$clean_name" "$category" "$has" "$n" >> m00115_survey.tsv
    echo "  [$i/$TOTAL] $code ($category): M00115=$has (${n} genes)"
    sleep 0.4
done < organism_list.tsv

echo ""
echo "=========================================="
echo "Summary"
echo "=========================================="

python3 << 'PYEOF'
import pandas as pd
from pathlib import Path

OUT = Path("/mnt/d/ce_microbiome/results_compost/revision_kegg_expanded")
df = pd.read_csv(OUT / "m00115_survey.tsv", sep="\t")

print("\n===== ORGANISM-LEVEL SURVEY =====")
print(df.to_string(index=False))

print("\n===== BY CATEGORY =====")
summary = df.groupby("category").agg(
    n_organisms=("org_code", "count"),
    n_with_M00115=("has_M00115", lambda x: (x == "YES").sum()),
).reset_index()
summary["pct_with_M00115"] = (100 * summary["n_with_M00115"] / summary["n_organisms"]).round(1)
print(summary.to_string(index=False))

print("\n===== INTERPRETATION =====")
env_cats = ["Brucella_anthropi_environmental", "Ochrobactrum_environmental"]
path_cats = ["Brucella_classical_pathogen"]

env_df = df[df["category"].isin(env_cats)]
path_df = df[df["category"].isin(path_cats)]

if len(env_df) > 0:
    env_pct = 100 * (env_df["has_M00115"] == "YES").mean()
    print(f"Environmental (B.anthropi + Ochrobactrum, n={len(env_df)}): {env_pct:.0f}% carry M00115")
else:
    print("Environmental: no genomes found in KEGG")

if len(path_df) > 0:
    path_pct = 100 * (path_df["has_M00115"] == "YES").mean()
    print(f"Classical pathogens (n={len(path_df)}):                {path_pct:.0f}% carry M00115")
else:
    print("Pathogens: no genomes found")

# Also list which organisms carry M00115
m00115_pos = df[df["has_M00115"] == "YES"]
print(f"\nOrganisms carrying M00115 ({len(m00115_pos)}):")
for _, row in m00115_pos.iterrows():
    print(f"  {row['org_code']:6s}  {row['organism']}")
PYEOF

echo ""
echo "Output: $OUT/m00115_survey.tsv"
echo "=========================================="

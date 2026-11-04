import pandas as pd
from pathlib import Path
from docx import Document

# --- Supp Table 14 ---
SRC = Path("/mnt/d/ce_microbiome/results_compost/revision_kegg_expanded/neuroactive_module_survey.tsv")
DST = Path("/mnt/d/ce_microbiome/Tables_final/Supp_Table_14_Neuroactive_Precursors.xlsx")

df = pd.read_csv(SRC, sep="\t")
df = df.rename(columns={
    "module": "KEGG module ID",
    "description": "Description",
    "env": "Environmental (n=9)",
    "path": "Pathogenic (n=43)",
    "p": "Fisher p",
    "q_BH": "BH q-value",
})
with pd.ExcelWriter(DST, engine="openpyxl") as writer:
    df.to_excel(writer, sheet_name="Neuroactive survey", index=False)
    ws = writer.sheets["Neuroactive survey"]
    for col in ws.columns:
        max_len = max(len(str(c.value)) if c.value else 0 for c in col)
        ws.column_dimensions[col[0].column_letter].width = min(max_len + 3, 50)
print(f"✅ Saved {DST}")
print(f"   {len(df)} rows in Supp Table 14")

# --- Verify FINAL_v7 ---
doc = Document("/mnt/d/ce_microbiome/Manuscript_Final/Celegans_Microbiome_Manuscript_FINAL_v7.docx")
t = "\n".join(p.text for p in doc.paragraphs)

print("\n=== v7 sanity checks ===")
checks = {
    "New title":                     "Conserved Host Enrichment and Functional Divergence",
    "Neuroactive Results paragraph": "neurotransmitter precursor pathways",
    "Neuroactive Discussion":        "do not support a specific neuroactive",
    "Supp Table 14 referenced":      "Supplementary Table 14",
    "Tryptophan mentioned":          "Tryptophan biosynthesis",
    "GLasso trimmed":                "remains future work",
    "KEGG trimmed":                  "would be required to confirm in situ activity",
}
for label, pat in checks.items():
    print(f"  {'✅' if pat in t else '⚠️ MISSING'}: {label}")

print(f"\nWord count: {len(t.split()):,}")

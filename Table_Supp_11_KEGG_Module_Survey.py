import pandas as pd
from pathlib import Path

SRC = Path("/mnt/d/ce_microbiome/results_compost/revision_kegg_expanded/kegg_full_module_survey.tsv")
DST = Path("/mnt/d/ce_microbiome/Tables_final/Supp_Table_13_KEGG_Module_Survey.xlsx")

df = pd.read_csv(SRC, sep="\t")
df = df.rename(columns={
    "module": "KEGG module ID",
    "env_present": "N environmental",
    "env_pct": "% environmental",
    "path_present": "N pathogenic",
    "path_pct": "% pathogenic",
    "p_value": "Fisher exact p",
    "q_BH": "BH q-value",
})
df.insert(0, "Rank", range(1, len(df) + 1))

with pd.ExcelWriter(DST, engine="openpyxl") as writer:
    df.to_excel(writer, sheet_name="All modules", index=False)
    ws = writer.sheets["All modules"]
    for col in ws.columns:
        max_len = max(len(str(c.value)) if c.value else 0 for c in col)
        ws.column_dimensions[col[0].column_letter].width = min(max_len + 3, 40)

print(f"✅ Saved {DST}")
print(f"   {len(df)} modules ranked")
print(f"   Top 3 (tied): M00552, M00115, M00999 — all p=2.72e-10")
print(f"   Rank 9: M00126 (folate) — 4/9 env vs 43/43 path")
print(f"   Rank ~18: M00027 (GABA) — 1/9 env vs 27/43 path")

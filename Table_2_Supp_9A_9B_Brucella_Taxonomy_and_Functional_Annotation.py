import pandas as pd
from pathlib import Path

OUT = Path("/mnt/d/ce_microbiome/Tables_final")
OUT.mkdir(exist_ok=True)

rows = [
    # Module, Description, oan, baa, bme, Interpretation
    ["M00083", "Fatty acid biosynthesis", "YES (14)", "YES (13)", "YES (13)", "Shared — no discriminator"],
    ["M00115", "NAD+ biosynthesis from aspartate", "YES (6)", "no", "no", "Environmental-specific"],
    ["M00126", "Folate biosynthesis", "no", "YES (8)", "YES (8)", "Pathogen-specific"],
    ["M00027", "GABA shunt", "no", "YES (6)", "YES (6)", "Pathogen-specific"],
    ["M00122", "Cobalamin (B12) biosynthesis", "no", "no", "no", "Absent in all (KEGG annotation gap)"],
    ["M00912", "NAD+ biosynthesis from tryptophan", "no", "no", "no", "Absent in all"],
    ["K01592", "Tyrosine decarboxylase (tyramine)", "YES (1)", "YES (1)", "YES (1)", "Shared — no discriminator"],
    ["oan00640", "Propanoate metabolism (SCFA)", "YES", "YES (baa00640)", "YES (bme00640)", "Shared — no discriminator"],
]

df = pd.DataFrame(rows, columns=[
    "KEGG_ID", "Pathway / Module", 
    "Ochrobactrum_anthropi (oan)", 
    "Brucella_abortus (baa)", 
    "Brucella_melitensis (bme)", 
    "Interpretation"
])

xlsx = OUT / "Supp_Table_5_Brucella_Ochrobactrum_Functional_Annotation.xlsx"
with pd.ExcelWriter(xlsx, engine="openpyxl") as writer:
    df.to_excel(writer, sheet_name="Supp Table 5", index=False)
    ws = writer.sheets["Supp Table 5"]
    for col in ws.columns:
        max_len = max(len(str(cell.value)) if cell.value else 0 for cell in col)
        ws.column_dimensions[col[0].column_letter].width = min(max_len + 3, 45)

print(f"✅ Saved {xlsx}")
print(df.to_string(index=False))

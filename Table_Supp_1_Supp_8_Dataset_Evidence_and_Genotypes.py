import pandas as pd
from pathlib import Path

OUT = Path("/mnt/d/ce_microbiome/Tables_final")
OUT.mkdir(exist_ok=True)

# ---- Supp Table 10: Dataset × Evidence matrix ----
supp10 = pd.DataFrame([
    {"Dataset": "CeMbio", "Accession": "PRJEB37101/35", "Detection of Brucella/Ochrobactrum": "Yes (0.17%)",
     "Direct worm-vs-control enrichment evidence": "Not established (no paired control)",
     "Genotype association": "Not applicable",
     "Role in study": "Presence baseline; defined community"},
    {"Dataset": "BIGbiome", "Accession": "PRJNA540192", "Detection of Brucella/Ochrobactrum": "Yes (24.35%)",
     "Direct worm-vs-control enrichment evidence": "Not established (no paired control)",
     "Genotype association": "Yes (43 genotypes)",
     "Role in study": "Host genotype analysis"},
    {"Dataset": "Aging", "Accession": "PRJNA982115", "Detection of Brucella/Ochrobactrum": "Absent (0.00%)",
     "Direct worm-vs-control enrichment evidence": "Not applicable (no control group)",
     "Genotype association": "Not applicable",
     "Role in study": "OP50-diet negative control"},
    {"Dataset": "Dirksen", "Accession": "ERP014569", "Detection of Brucella/Ochrobactrum": "Yes (6.09%)",
     "Direct worm-vs-control enrichment evidence": "Yes (worm vs. lawn; lfc = 6.28)",
     "Genotype association": "Not tested here",
     "Role in study": "Ground-truth validation"},
    {"Dataset": "Environ", "Accession": "PRJNA1002096", "Detection of Brucella/Ochrobactrum": "Yes (9.37%)",
     "Direct worm-vs-control enrichment evidence": "Not established (no paired control)",
     "Genotype association": "Not applicable",
     "Role in study": "High-diversity reference"},
    {"Dataset": "Compost", "Accession": "PRJNA856419", "Detection of Brucella/Ochrobactrum": "Yes (4.34%)",
     "Direct worm-vs-control enrichment evidence": "Yes (worm vs. soil; 4 environments)",
     "Genotype association": "Not applicable",
     "Role in study": "Primary discovery dataset"},
])
supp10.to_excel(OUT / "Supp_Table_10_Dataset_Evidence_Matrix.xlsx", index=False)
print(f"✅ Supp Table 10: {len(supp10)} datasets")

# ---- Supp Table 11: Expanded KEGG survey ----
kegg_file = Path("/mnt/d/ce_microbiome/results_compost/revision_kegg_expanded/m00115_survey_categorized.tsv")
if kegg_file.exists():
    df = pd.read_csv(kegg_file, sep="\t")
    df = df.rename(columns={
        "org_code": "KEGG organism code",
        "organism": "Organism",
        "clean_category": "Category",
        "has_M00115": "Has M00115",
        "n_genes_M00115": "N genes in M00115",
    })
    # Exclude outgroup Luteibacter from main table
    df = df[df["Category"] != "exclude_outgroup"]
    df = df[["Category", "Organism", "KEGG organism code", "Has M00115", "N genes in M00115"]]
    df = df.sort_values(["Category", "Organism"])

    summary = df.groupby("Category").agg(
        n_genomes=("KEGG organism code", "count"),
        n_with_M00115=("Has M00115", lambda x: (x == "YES").sum()),
    ).reset_index()
    summary["pct_with_M00115"] = (100 * summary["n_with_M00115"] / summary["n_genomes"]).round(1)
    summary = summary.rename(columns={
        "Category": "Category", "n_genomes": "N genomes",
        "n_with_M00115": "N with M00115", "pct_with_M00115": "% with M00115",
    })

    with pd.ExcelWriter(OUT / "Supp_Table_11_KEGG_M00115_Expanded.xlsx", engine="openpyxl") as writer:
        summary.to_excel(writer, sheet_name="Summary", index=False)
        df.to_excel(writer, sheet_name="Full survey", index=False)
        for sn in writer.sheets:
            ws = writer.sheets[sn]
            for col in ws.columns:
                max_len = max((len(str(c.value)) if c.value else 0) for c in col)
                ws.column_dimensions[col[0].column_letter].width = min(max_len + 3, 60)
    print(f"✅ Supp Table 11: {len(df)} genomes across {len(summary)} categories")
    print(f"\nCategory summary:")
    print(summary.to_string(index=False))
else:
    print(f"⚠️ KEGG file not found: {kegg_file}")

print(f"\nFiles saved to: {OUT}")

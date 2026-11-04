import subprocess, time, numpy as np, pandas as pd
from pathlib import Path
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage

OUT = Path("/mnt/d/ce_microbiome/results_compost/revision_asr")
OUT.mkdir(parents=True, exist_ok=True)
cache = OUT / "ko_presence.tsv"

s = pd.read_csv("/mnt/d/ce_microbiome/results_compost/revision_kegg_expanded/m00115_survey_categorized.tsv", sep="\t")
s = s[s.clean_category != "exclude_outgroup"].drop_duplicates(subset=["org_code"]).reset_index(drop=True)
cds = s.org_code.tolist()
m = dict(zip(s.org_code, s.has_M00115.astype(str).str.upper()=="YES"))

if cache.exists():
    print("Using cached KO table...")
    ko = pd.read_csv(cache, sep="\t", index_col=0)
    kos = list(ko.columns)
else:
    print(f"Fetching KOs for {len(cds)} genomes...")
    d = {}
    for i, c in enumerate(cds):
        try:
            r = subprocess.check_output(["curl", "-s", "--max-time", "15", f"https://rest.kegg.jp/link/ko/{c}"]).decode()
            d[c] = set(l.split("\t")[1].replace("ko:","").strip() for l in r.split("\n") if "\t" in l and l.split("\t")[1].startswith("ko:"))
            if (i+1)%10==0: print(f"  {i+1}/{len(cds)}")
            time.sleep(0.3)
        except: d[c] = set()
    kos = sorted(set().union(*d.values()))
    arr = np.zeros((len(cds), len(kos)), dtype=int)
    for i, c in enumerate(cds):
        for k in d.get(c, []):
            if k in kos:
                arr[i, kos.index(k)] = 1
    ko = pd.DataFrame(arr, index=cds, columns=kos)
    ko.to_csv(cache, sep="\t")

# Use numpy array directly (avoid pandas .loc)
arr = ko.values.astype(bool)
n = len(cds)
D = np.zeros((n, n))
for i in range(n):
    for j in range(i+1, n):
        u = (arr[i] | arr[j]).sum()
        D[i, j] = D[j, i] = 1 - ((arr[i] & arr[j]).sum() / u if u else 1)

Z = linkage(squareform(D, checks=False), method="average")
np.save(OUT / "asr_linkage.npy", Z)

cur = {i: {1 if m[c] else 0} for i, c in enumerate(cds)}
nid = n
for r in Z:
    a, b = int(r[0]), int(r[1])
    inter = cur[a] & cur[b]
    cur[nid] = inter if inter else (cur[a] | cur[b])
    nid += 1

print("\n=== ROOT STATE (1=present, 0=absent) ===")
print(cur[nid - 1])

import sys
import pandas as pd

# Get input and output file paths from command line arguments
inmark = sys.argv[1]
outmark = sys.argv[2]
lmarkers = sys.argv[3]

# Save markers gene list as markers.xlsx
markers1 = pd.read_csv(lmarkers, sep="\t", encoding='latin1')
markers1.to_excel(outmark, sheet_name="markers", index=False)

# Create dictionary from markers
d = {}
for i in open(lmarkers, encoding='latin1'):
    line1 = i.strip().split("\t")
    if len(line1) == 2:
        w = line1[0]
        code = line1[1]
        d.setdefault(w, []).append(code)

e = {i: str(j) for i, j in d.items()}
print(e)
markers = pd.Series(e).to_frame('subtype')
markers.reset_index(inplace=True)
markers = markers.rename(columns={'index': 'gene'})
print(markers)

# Save output from FindMarkers as marker_genes_list.xlsx
infile = pd.read_csv(inmark, sep="\t", header=0, encoding='unicode_escape')

# Get unique cluster identifiers
unique_clusters = infile['cluster'].unique()
print(f"Unique clusters: {unique_clusters}")

# Process each cluster and save to the Excel file
for cluster in unique_clusters:
    a = infile[infile["cluster"] == cluster]
    b = a.merge(markers, on="gene", how="left")
    print(b)
    with pd.ExcelWriter(outmark, engine="openpyxl", mode="a") as f:
        b.to_excel(f, sheet_name=f"cluster_{cluster}", index=False)

# Optional section for additional processing and saving to a different file
"""
xl = pd.ExcelFile('marker_genes_all.xlsx')
list1 = xl.sheet_names[1:]

markers1 = pd.read_excel(open('marker_genes_all.xlsx', 'rb'), sheet_name='markers')

markers1.to_excel("marker_genes_annotated.xlsx", sheet_name="markers", index=False)

for i in list1:
    a = pd.read_excel(open('marker_genes_all.xlsx', 'rb'), sheet_name=i)
    print(a)
    b = a.merge(markers, on="gene", how="left")
    b.sort_values('avg_log2FC', ascending=False, inplace=True)

    with pd.ExcelWriter("marker_genes_annotated.xlsx", engine="openpyxl", mode="a") as f:
        b.to_excel(f, sheet_name=i, index=False)
"""

# Command to run the script
# python /path/to/annotate.py /path/to/metaprograms_genes_cafs.txt /path/to/metaprograms_genes_cafs_annot_CAFs_Our_andAdel.xlsx /path/to/CAFs_Our_andAdel.txt

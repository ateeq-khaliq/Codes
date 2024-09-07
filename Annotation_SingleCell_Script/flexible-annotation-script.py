import sys
import pandas as pd
import os
from collections import Counter

def read_input_file(file_path):
    """Read input file based on its extension."""
    _, ext = os.path.splitext(file_path)
    if ext.lower() == '.csv':
        return pd.read_csv(file_path, encoding='unicode_escape')
    elif ext.lower() == '.txt':
        return pd.read_csv(file_path, sep="\t", encoding='unicode_escape')
    else:
        raise ValueError(f"Unsupported file format: {ext}. Please use .csv or .txt files.")

def convert_gene_names(gene_series, species):
    """Convert gene names based on species."""
    if species.lower() == 'mouse':
        return gene_series.str.capitalize()
    elif species.lower() == 'human':
        return gene_series.str.upper()
    else:
        raise ValueError("Invalid species. Please specify 'mouse' or 'human'.")

def print_file_info(df, file_name):
    """Print information about the DataFrame."""
    print(f"\nFile: {file_name}")
    print("Columns:", df.columns.tolist())
    print("Shape:", df.shape)
    print("First few rows:")
    print(df.head())

# Get input and output file paths from command line arguments
inmark = sys.argv[1]
outmark = sys.argv[2]
lmarkers = sys.argv[3]
species = sys.argv[4]  # New argument for species

# Read and print info about input files
markers1 = read_input_file(lmarkers)
print_file_info(markers1, lmarkers)

infile = read_input_file(inmark)
print_file_info(infile, inmark)

# Use the correct column names for markers file
gene_column_markers = 'Genes'
subtype_column_markers = 'Subtype'

if gene_column_markers not in markers1.columns or subtype_column_markers not in markers1.columns:
    print(f"Error: Expected columns '{gene_column_markers}' and '{subtype_column_markers}' not found in the markers file.")
    print("Available columns:", markers1.columns.tolist())
    sys.exit(1)

# Convert gene names and save markers
markers1[gene_column_markers] = convert_gene_names(markers1[gene_column_markers], species)
markers1.to_excel(outmark, sheet_name="markers", index=False)

# Create dictionary from markers
d = {}
for _, row in markers1.iterrows():
    gene = row[gene_column_markers]
    subtype = row[subtype_column_markers]
    d.setdefault(gene, []).append(subtype)

e = {i: str(j) for i, j in d.items()}
print("Marker dictionary (first 5 items):", dict(list(e.items())[:5]))

markers = pd.Series(e).to_frame(subtype_column_markers)
markers.reset_index(inplace=True)
markers = markers.rename(columns={'index': gene_column_markers})
print("Processed markers (first 5 rows):")
print(markers.head())

# Use the correct column name for gene in the input file
gene_column_infile = 'gene'

if gene_column_infile not in infile.columns:
    print(f"Error: Expected column '{gene_column_infile}' not found in the input file.")
    print("Available columns:", infile.columns.tolist())
    sys.exit(1)

# Convert gene names in the input file
infile[gene_column_infile] = convert_gene_names(infile[gene_column_infile], species)

# Get unique cluster identifiers
if 'cluster' not in infile.columns:
    print("Error: 'cluster' column not found in the input file.")
    print("Available columns:", infile.columns.tolist())
    sys.exit(1)

unique_clusters = infile['cluster'].unique()
print(f"Unique clusters: {unique_clusters}")

# Initialize counters and lists for summary
total_genes = 0
total_annotated_genes = 0
cluster_gene_counts = {}
cluster_annotated_counts = {}
top_genes_per_cluster = {}

# Process each cluster and save to the Excel file
for cluster in unique_clusters:
    a = infile[infile["cluster"] == cluster]
    b = a.merge(markers, left_on=gene_column_infile, right_on=gene_column_markers, how="left")
    
    # Count genes and annotated genes
    cluster_gene_counts[cluster] = len(b)
    cluster_annotated_counts[cluster] = b[subtype_column_markers].notna().sum()
    total_genes += cluster_gene_counts[cluster]
    total_annotated_genes += cluster_annotated_counts[cluster]
    
    # Get top 5 genes by avg_log2FC
    if 'avg_log2FC' in b.columns:
        top_genes = b.nlargest(5, 'avg_log2FC')[gene_column_infile].tolist()
    else:
        top_genes = ['avg_log2FC column not found']
    top_genes_per_cluster[cluster] = top_genes
    
    print(f"Cluster {cluster} processed")
    with pd.ExcelWriter(outmark, engine="openpyxl", mode="a") as f:
        b.to_excel(f, sheet_name=f"cluster_{cluster}", index=False)

# Generate summary information
summary = f"""
Annotation Summary:
-------------------
Species: {species.capitalize()}
Total number of genes processed: {total_genes}
Total number of annotated genes: {total_annotated_genes}
Overall annotation rate: {(total_annotated_genes/total_genes)*100:.2f}%

Cluster-wise Information:
-------------------------
"""

for cluster in unique_clusters:
    summary += f"""
Cluster {cluster}:
    Total genes: {cluster_gene_counts[cluster]}
    Annotated genes: {cluster_annotated_counts[cluster]}
    Annotation rate: {(cluster_annotated_counts[cluster]/cluster_gene_counts[cluster])*100:.2f}%
    Top 5 genes by avg_log2FC: {', '.join(top_genes_per_cluster[cluster])}
"""

# Add summary to Excel file
with pd.ExcelWriter(outmark, engine="openpyxl", mode="a") as f:
    summary_df = pd.DataFrame([summary], columns=['Summary'])
    summary_df.to_excel(f, sheet_name="Summary", index=False)

print("Annotation complete. Results and summary saved to:", outmark)
print("\nSummary:")
print(summary)

# Command to run the script
# python /path/to/annotate.py /path/to/input_file.csv /path/to/output_file.xlsx /path/to/marker_genes.txt species

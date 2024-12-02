
## Preparing Input Files for cNMF from Seurat

# Load required packages
suppressPackageStartupMessages({
    library(data.table)
    library(Matrix)
    library(Seurat)
})

# Function to prepare input files for cNMF
prepare_myeloid_input <- function(seurat_obj, 
                                output_dir = "./myeloid_cnmf", 
                                min_cells = 3, 
                                min_features = 200) {
    
    # Create main output directory
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }
    
    # Create filtered directory
    filtered_dir <- file.path(output_dir, "filtered")
    if (!dir.exists(filtered_dir)) {
        dir.create(filtered_dir, recursive = TRUE)
    }
    
    # Filter the Seurat object
    print("Filtering cells and features...")
    seurat_obj <- subset(seurat_obj, 
                        subset = nFeature_RNA > min_features)
    
    # Get counts matrix
    print("Extracting counts matrix...")
    counts <- GetAssayData(seurat_obj, slot = "counts", assay = "RNA")
    barcodes <- colnames(counts)
    gene_names <- rownames(counts)
    
    # Save counts matrix in MTX format
    print("Saving counts matrix...")
    writeMM(counts, file.path(filtered_dir, "matrix.mtx"))
    
    # Save cell barcodes
    print("Saving cell barcodes...")
    write.table(as.data.frame(barcodes), 
                file.path(filtered_dir, "barcodes.tsv"),
                col.names = FALSE, 
                row.names = FALSE, 
                sep = "\t")
    
    # Save feature names
    print("Saving gene names...")
    features <- data.frame(
        "gene_id" = gene_names,
        "gene_name" = gene_names,
        type = "Gene Expression"
    )
    write.table(features, 
                file.path(filtered_dir, "genes.tsv"),
                col.names = FALSE, 
                row.names = FALSE, 
                sep = "\t")
    
    # Save Seurat object
    print("Saving filtered Seurat object...")
    saveRDS(seurat_obj, file.path(output_dir, "filtered_seurat.rds"))
    
    # Save metadata about the processing
    metadata <- list(
        date = Sys.time(),
        original_cells = ncol(seurat_obj),
        original_genes = nrow(seurat_obj),
        min_cells = min_cells,
        min_features = min_features,
        filtered_cells = ncol(counts),
        filtered_genes = nrow(counts)
    )
    
    saveRDS(metadata, file.path(output_dir, "processing_metadata.rds"))
    
    # Create a README file
    readme_text <- sprintf(
"cNMF Input Files
================
Created on: %s

This directory contains input files for cNMF analysis:

./filtered/
- matrix.mtx: Sparse count matrix
- barcodes.tsv: Cell barcodes
- genes.tsv: Gene names

Additional files:
- filtered_seurat.rds: Filtered Seurat object
- processing_metadata.rds: Processing parameters and statistics

Original data:
- Cells: %d
- Genes: %d

Filtering parameters:
- Minimum cells per feature: %d
- Minimum features per cell: %d

Filtered data:
- Cells: %d
- Genes: %d
",
        Sys.time(),
        metadata$original_cells,
        metadata$original_genes,
        metadata$min_cells,
        metadata$min_features,
        metadata$filtered_cells,
        metadata$filtered_genes
    )
    
    writeLines(readme_text, file.path(output_dir, "README.txt"))
    
    print(paste("Input files prepared and saved in:", output_dir))
    print("Directory structure created:")
    print(paste("  ", output_dir))
    print(paste("  ├── filtered/"))
    print(paste("  │   ├── matrix.mtx"))
    print(paste("  │   ├── barcodes.tsv"))
    print(paste("  │   └── genes.tsv"))
    print(paste("  ├── filtered_seurat.rds"))
    print(paste("  ├── processing_metadata.rds"))
    print(paste("  └── README.txt"))
    
    return(list(
        output_dir = output_dir,
        filtered_dir = filtered_dir,
        metadata = metadata
    ))
}

# Example usage:
# Assuming 'myeloid' is your Seurat object:
# results <- prepare_myeloid_input(myeloid)

###


## run actualy cNMF

# Load required packages
suppressPackageStartupMessages({
    library(data.table)
    library(Matrix)
    library(Seurat)
    library(dplyr)
})

# Set directories where your files are
data_dir = './myeloid_cnmf/'  # Main output directory
filtered_dir = './myeloid_cnmf/filtered/'  # Where your files are located

# Set run name
runname = "myeloid_cNMF"

# Step 1: Prepare step
print("Running prepare step...")
cmd = paste("cnmf prepare --output-dir", data_dir,
            "--name", runname,
            "-c", paste0(filtered_dir, 'matrix.mtx'),
            "--max-nmf-iter 2000", 
            "-k 5 6 7 8 9 10 --n-iter 20", sep=" ")
print(cmd)
system(cmd)

# Step 2: Factorization step 
print("Running factorization step...")
cmd = paste("cnmf factorize --output-dir", data_dir,
            "--name", runname,
            "--worker-index 0 --total-workers 1", sep=" ")
print(cmd)
system(cmd)

# Step 3: Combine results
print("Combining results...")
cmd = paste("cnmf combine --output-dir", data_dir,
            "--name", runname, sep=" ")
print(cmd)
system(cmd)

# Step 4: Generate k selection plot
print("Generating k selection plot...")
cmd = paste("cnmf k_selection_plot --output-dir", data_dir,
            "--name", runname, sep=" ")
print(cmd)
system(cmd)

# Now you should examine the k selection plot before proceeding
print("Please examine the k selection plot and choose optimal k")
print(paste("Plot saved at:", file.path(data_dir, runname, paste0(runname, ".k_selection.png"))))

# After examining plot, run consensus clustering
k_value = 10  # This should be based on your k selection plot
print(paste("Running consensus clustering with k =", k_value))
cmd = paste("cnmf consensus --output-dir", data_dir,
            "--name", runname,
            '--components', k_value,
            '--local-density-threshold', 0.1,
            '--show-clustering', sep=" ")
print(cmd)
system(cmd)

# Load and process results
print("Loading results...")
usage_file <- file.path(data_dir, runname, 
                       paste(runname, "usages", paste0("k_", k_value), "dt_0_1", 'consensus', 'txt', sep="."))
spectra_score_file <- file.path(data_dir, runname, 
                               paste(runname, "gene_spectra_score", paste0("k_", k_value), "dt_0_1", 'txt', sep="."))

usage <- read.table(usage_file, sep='\t', row.names=1, header=TRUE)
spectra_score <- read.table(spectra_score_file, sep='\t', row.names=1, header=TRUE)

# Normalize usage matrix
print("Normalizing usage matrix...")
usage_norm <- as.data.frame(t(apply(usage, 1, function(x) x / sum(x))))

# Get top genes for each program
print("Getting top genes for each program...")
get_top_genes <- function(row, ntop=20) {
    top_indices <- order(row, decreasing = TRUE)[1:ntop]
    return(colnames(spectra_score)[top_indices])
}

top_genes <- apply(spectra_score, 1, get_top_genes)
top_genes_df <- as.data.frame(top_genes)

# Save results
print("Saving results...")
write.csv(usage_norm, file=file.path(data_dir, "normalized_usage.csv"))
write.csv(top_genes_df, file=file.path(data_dir, "top_genes_per_program.csv"))

print("Analysis complete!")
print(paste("Results saved in:", data_dir))

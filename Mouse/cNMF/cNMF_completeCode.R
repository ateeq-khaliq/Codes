
## Preparing Input Files for cNMF from Seurat

# Load required packages
suppressPackageStartupMessages({
    library(data.table)
    library(Matrix)
    library(Seurat)
})

# Function to prepare input files for cNMF
prepare_myeloid_input <- function(myeloid, 
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
    myeloid <- subset(myeloid, 
                        subset = nFeature_RNA > min_features)
    
    # Get counts matrix
    print("Extracting counts matrix...")
    counts <- GetAssayData(myeloid, slot = "counts", assay = "RNA")
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
    saveRDS(myeloid, file.path(output_dir, "filtered_seurat.rds"))
    
    # Save metadata about the processing
    metadata <- list(
        date = Sys.time(),
        original_cells = ncol(myeloid),
        original_genes = nrow(myeloid),
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

###


# Add program scores to Seurat object
usage_norm <- as.data.frame(t(apply(usage, 1, function(x) x / sum(x))))
myeloid@meta.data <- cbind(myeloid@meta.data, usage_norm)

# Visualize program scores on UMAP
FeaturePlot(myeloid, 
            features = colnames(usage_norm),
            ncol = 3)

# Violin plots across cell types
VlnPlot(myeloid, 
        features = colnames(usage_norm),
        group.by = "cell_type",
        ncol = 3)

##

# Add program scores to Seurat object
program_cols <- paste0("Program_", 1:ncol(usage_norm))
colnames(usage_norm) <- program_cols
myeloid@meta.data <- cbind(myeloid@meta.data, usage_norm)

# Calculate average program usage per cluster
avg_usage_by_cluster <- aggregate(usage_norm, 
                                by = list(cluster = myeloid$seurat_clusters), 
                                FUN = mean)

# Create heatmap of program-cluster associations
library(pheatmap)
library(viridis)

# Prepare matrix for heatmap
cluster_program_matrix <- t(avg_usage_by_cluster[,-1])
colnames(cluster_program_matrix) <- paste0("Cluster_", avg_usage_by_cluster$cluster)

# Generate heatmap
pheatmap(cluster_program_matrix,
         scale = "row",
         clustering_method = "ward.D2",
         color = viridis(100),
         main = "Program Usage Across Clusters",
         fontsize = 10)

# Statistical test for program enrichment in clusters
library(stats)
test_enrichment <- function(program, clusters) {
  kruskal.test(program ~ clusters)$p.value
}

enrichment_pvals <- apply(usage_norm, 2, function(x) {
  test_enrichment(x, myeloid$seurat_clusters)
})

print(data.frame(Program = program_cols, 
                P_value = enrichment_pvals))


###
#program_correlations.pdf: Correlation heatmap
#pathway_enrichment.pdf: Pathway dotplots
#pathway_analysis_report.txt: Detailed pathway results

# Save as 'myeloid_pathway_analysis.R'

# Install and load required packages
#if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")

#BiocManager::install(c("clusterProfiler", "org.Mm.eg.db", "enrichplot"))

library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(dplyr)

library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(org.Mm.eg.db)
library(dplyr)

# Function to perform pathway analysis
analyze_program_pathways <- function(gene_list, program_name) {
   gene_ids <- mapIds(org.Mm.eg.db, 
                     keys = gene_list,
                     keytype = "SYMBOL",
                     column = "ENTREZID")
   
   go_bp <- enrichGO(gene = gene_ids[!is.na(gene_ids)],
                     OrgDb = org.Mm.eg.db,
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05)
   
   return(go_bp)
}

# Function to create correlation heatmap
create_pathway_correlation_plot <- function(correlation_matrix, 
                                         pathway_annotations, 
                                         output_file,
                                         width = 10,
                                         height = 8) {
   
   colors <- colorRampPalette(c("#313695", "#FFFFFF", "#A50026"))(100)
   
   annot_colors <- list(
       Program_Type = setNames(viridis(length(unique(pathway_annotations))), 
                             unique(pathway_annotations))
   )
   
   annotation_df <- data.frame(
       Program_Type = pathway_annotations,
       row.names = colnames(correlation_matrix)
   )
   
   pdf(output_file, width = width, height = height, useDingbats = FALSE)
   
   pheatmap(correlation_matrix,
            display_numbers = TRUE,
            number_format = "%.2f",
            number_color = "black",
            fontsize_number = 8,
            cluster_rows = TRUE,
            cluster_cols = TRUE,
            clustering_distance_rows = "euclidean",
            clustering_method = "complete",
            labels_row = pathway_annotations,
            labels_col = pathway_annotations,
            annotation_row = annotation_df,
            annotation_col = annotation_df,
            annotation_colors = annot_colors,
            fontsize = 10,
            fontsize_row = 8,
            fontsize_col = 8,
            color = colors,
            breaks = seq(-1, 1, length.out = 101),
            border_color = "white",
            main = "Program Correlations with Pathway Annotations",
            legend = TRUE,
            legend_breaks = c(-1, -0.5, 0, 0.5, 1),
            legend_labels = c("-1.0", "-0.5", "0", "0.5", "1.0")
   )
   
   dev.off()
}

# Modify the pathway analysis function to include Hallmark and KEGG
analyze_program_pathways <- function(gene_list, program_name) {
    # Convert gene symbols to ENTREZ IDs
    gene_ids <- mapIds(org.Mm.eg.db, 
                      keys = gene_list,
                      keytype = "SYMBOL",
                      column = "ENTREZID")
    gene_ids <- gene_ids[!is.na(gene_ids)]
    
    # GO BP Analysis
    go_bp <- enrichGO(gene = gene_ids,
                      OrgDb = org.Mm.eg.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05)
    
    # KEGG Analysis
    kegg <- enrichKEGG(gene = gene_ids,
                       organism = 'mmu',
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05)
    
    # Hallmark Analysis
    hallmark_gene_sets = msigdbr(species = "Mus musculus", category = "H")
    hallmark <- enricher(gene = gene_ids,
                        TERM2GENE = dplyr::select(hallmark_gene_sets, gs_name, entrez_gene),
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05)
    
    return(list(
        go = go_bp,
        kegg = kegg,
        hallmark = hallmark
    ))
}

# Modified main analysis pipeline
analyze_myeloid_programs <- function(myeloid, usage_norm, top_genes_df) {
    # Calculate correlations between programs
    program_correlations <- cor(usage_norm, method = "spearman")
    
    # Run pathway analysis for each program
    pathway_results <- list()
    for(i in 1:ncol(top_genes_df)) {
        genes <- top_genes_df[,i]
        pathway_results[[i]] <- analyze_program_pathways(genes, paste0("Program_", i))
    }
    
    # Get top pathways and create annotations (using GO BP as primary annotation)
    pathway_annotations <- sapply(pathway_results, function(x) {
        if(nrow(x$go@result) > 0) {
            return(x$go@result$Description[1])
        } else {
            return("No significant pathways")
        }
    })
    
    # Create main correlation plot
    create_pathway_correlation_plot(
        correlation_matrix = program_correlations,
        pathway_annotations = pathway_annotations,
        output_file = "program_correlations.pdf"
    )
    
    # Create pathway enrichment plots for all three analyses
    pdf("pathway_enrichment_all.pdf", width = 12, height = 8)
    for(i in 1:length(pathway_results)) {
        # GO Plot
        if(nrow(pathway_results[[i]]$go@result) > 0) {
            print(dotplot(pathway_results[[i]]$go, 
                         showCategory = 10, 
                         title = paste0("Program ", i, " - GO Pathways"),
                         font.size = 10) +
                  theme_minimal() +
                  theme(axis.text.y = element_text(size = 8)))
        }
        
        # KEGG Plot
        if(nrow(pathway_results[[i]]$kegg@result) > 0) {
            print(dotplot(pathway_results[[i]]$kegg, 
                         showCategory = 10, 
                         title = paste0("Program ", i, " - KEGG Pathways"),
                         font.size = 10) +
                  theme_minimal() +
                  theme(axis.text.y = element_text(size = 8)))
        }
        
        # Hallmark Plot
        if(nrow(pathway_results[[i]]$hallmark@result) > 0) {
            print(dotplot(pathway_results[[i]]$hallmark, 
                         showCategory = 10, 
                         title = paste0("Program ", i, " - Hallmark Pathways"),
                         font.size = 10) +
                  theme_minimal() +
                  theme(axis.text.y = element_text(size = 8)))
        }
    }
    dev.off()
    
    # Create detailed report
    sink("pathway_analysis_report.txt")
    cat("Myeloid Program Pathway Analysis Report\n")
    cat("======================================\n\n")
    
    for(i in 1:length(pathway_results)) {
        cat(sprintf("\nProgram %d:\n", i))
        cat("-------------\n")
        
        cat("\nTop GO Pathways:\n")
        if(nrow(pathway_results[[i]]$go@result) > 0) {
            top_paths <- head(pathway_results[[i]]$go@result, 5)
            print(top_paths[, c("Description", "pvalue", "p.adjust", "Count")])
        } else {
            cat("No significant GO pathways found\n")
        }
        
        cat("\nTop KEGG Pathways:\n")
        if(nrow(pathway_results[[i]]$kegg@result) > 0) {
            top_paths <- head(pathway_results[[i]]$kegg@result, 5)
            print(top_paths[, c("Description", "pvalue", "p.adjust", "Count")])
        } else {
            cat("No significant KEGG pathways found\n")
        }
        
        cat("\nTop Hallmark Pathways:\n")
        if(nrow(pathway_results[[i]]$hallmark@result) > 0) {
            top_paths <- head(pathway_results[[i]]$hallmark@result, 5)
            print(top_paths[, c("Description", "pvalue", "p.adjust", "Count")])
        } else {
            cat("No significant Hallmark pathways found\n")
        }
    }
    sink()
    
    return(list(
        correlations = program_correlations,
        pathway_results = pathway_results,
        annotations = pathway_annotations
    ))
}

# Usage example:
# results <- analyze_myeloid_programs(myeloid, usage_norm, top_genes_df)

#####

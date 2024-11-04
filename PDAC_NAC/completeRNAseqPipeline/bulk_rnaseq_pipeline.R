#!/usr/bin/env Rscript

# ===============================================
# RNA-seq Data Processing and Analysis Pipeline
# Version: 1.0
# Created: 2024-11-04
# Description: Comprehensive pipeline for RNA-seq count data processing,
#              including QC, normalization, and visualization
# ===============================================

# Load required libraries with error checking
required_packages <- c("edgeR", "ggplot2", "pheatmap", "reshape2", "dplyr")

for (package in required_packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
        stop(paste("Package", package, "is required but not installed. Please install it using install.packages('", package, "')"))
    }
    library(package, character.only = TRUE)
}

# Function to create directories safely
create_dir <- function(path) {
    if (!dir.exists(path)) {
        dir.create(path, recursive = TRUE)
    }
}

# 1. Data Reading and Preparation
prepare_count_data <- function(file_path) {
    message("Reading count data...")
    
    # Check if file exists
    if (!file.exists(file_path)) {
        stop("Input file does not exist: ", file_path)
    }
    
    # Read the count matrix
    counts <- read.csv(file_path, row.names = 1, header = TRUE)
    
    # Basic data validation
    required_cols <- c("Gene_Name", "Gene_Type")
    if (!all(required_cols %in% colnames(counts))) {
        stop("Input file must contain Gene_Name and Gene_Type columns")
    }
    
    # Verify data structure
    cat("Initial data verification:\n")
    cat("Dimensions:", dim(counts), "\n")
    cat("First few column names:", head(colnames(counts)), "\n")
    cat("First few row names:", head(rownames(counts)), "\n")
    
    # Separate gene annotation and count data
    gene_info <- data.frame(
        Gene_Name = counts$Gene_Name,
        Gene_Type = counts$Gene_Type,
        row.names = rownames(counts)
    )
    
    # Extract only count columns (starting with M)
    count_columns <- grep("^M", colnames(counts))
    if (length(count_columns) == 0) {
        stop("No sample columns found (columns starting with 'M')")
    }
    
    count_matrix <- as.matrix(counts[, count_columns])
    mode(count_matrix) <- "numeric"
    
    return(list(
        counts = count_matrix,
        gene_info = gene_info
    ))
}

# 2. Create filtered DGEList
create_filtered_dge <- function(count_data) {
    message("Filtering low expression genes...")
    
    # Create DGEList object
    dge <- DGEList(counts = count_data$counts)
    
    # Calculate CPM
    cpm <- cpm(dge)
    
    # Filter low expression genes (CPM > 1 in at least 3 samples)
    keep <- rowSums(cpm > 1) >= 3
    dge_filtered <- dge[keep, ]
    
    # Print filtering summary
    cat("\nFiltering summary:\n")
    cat("Original genes:", nrow(count_data$counts), "\n")
    cat("Genes passing filter:", sum(keep), "\n")
    cat("Percentage retained:", round(sum(keep)/nrow(count_data$counts)*100, 2), "%\n")
    
    # Add gene info to the filtered DGE object
    dge_filtered$genes <- count_data$gene_info[rownames(dge_filtered), ]
    
    return(dge_filtered)
}

# 3. Process duplicates
process_duplicates <- function(count_matrix) {
    message("Processing duplicate genes...")
    
    # Create output directory for QC
    create_dir("duplicate_genes_qc")
    
    # Add Gene_ID column from rownames
    count_matrix$Gene_ID <- rownames(count_matrix)
    
    # Find duplicated genes
    dups <- duplicated(count_matrix$Gene_Name) | 
            duplicated(count_matrix$Gene_Name, fromLast = TRUE)
    dup_genes <- unique(count_matrix$Gene_Name[dups])
    
    if (length(dup_genes) == 0) {
        message("No duplicate genes found")
        return(list(
            details = NULL,
            resolution_strategy = NULL,
            apply_resolution = function(x) x
        ))
    }
    
    # Analyze duplicates
    sample_cols <- setdiff(colnames(count_matrix), 
                         c("Gene_Name", "Gene_Type", "Gene_ID"))
    
    # Process duplicate genes
    dup_details <- data.frame()
    for(gene in dup_genes) {
        gene_rows <- count_matrix[count_matrix$Gene_Name == gene, ]
        gene_metrics <- data.frame(
            Gene_Name = gene_rows$Gene_Name,
            Gene_ID = gene_rows$Gene_ID,
            Gene_Type = gene_rows$Gene_Type,
            Mean_Expression = rowMeans(gene_rows[, sample_cols]),
            Max_Expression = apply(gene_rows[, sample_cols], 1, max),
            Non_Zero_Samples = rowSums(gene_rows[, sample_cols] > 0)
        )
        dup_details <- rbind(dup_details, gene_metrics)
    }
    
    # Generate resolution strategy
    resolution_strategy <- data.frame()
    for(gene in unique(dup_details$Gene_Name)) {
        gene_rows <- dup_details[dup_details$Gene_Name == gene, ]
        if("protein_coding" %in% gene_rows$Gene_Type) {
            protein_coding_rows <- gene_rows[gene_rows$Gene_Type == "protein_coding", ]
            selected_row <- protein_coding_rows[which.max(protein_coding_rows$Mean_Expression), ]
            strategy <- "protein_coding"
        } else {
            selected_row <- gene_rows[which.max(gene_rows$Mean_Expression), ]
            strategy <- "highest_expression"
        }
        resolution_strategy <- rbind(resolution_strategy,
                                  data.frame(
                                      Gene_Name = gene,
                                      Selected_ID = selected_row$Gene_ID,
                                      Strategy = strategy,
                                      Mean_Expression = selected_row$Mean_Expression
                                  ))
    }
    
    # Save results
    write.csv(dup_details, "duplicate_genes_qc/duplicate_genes_details.csv", row.names = FALSE)
    write.csv(resolution_strategy, "duplicate_genes_qc/resolution_strategy.csv", row.names = FALSE)
    
    return(list(
        details = dup_details,
        resolution_strategy = resolution_strategy,
        apply_resolution = function(original_matrix) {
            original_matrix$Gene_ID <- rownames(original_matrix)
            resolved <- original_matrix[original_matrix$Gene_ID %in% resolution_strategy$Selected_ID, ]
            resolved$Gene_ID <- NULL
            return(resolved)
        }
    ))
}
# 4. Generate normalized counts and visualizations
create_visualizations <- function(counts_data, output_file) {
    message("Creating visualizations...")
    
    # Get sample columns
    sample_cols <- colnames(counts_data)[!colnames(counts_data) %in% c("Gene_Name", "Gene_Type", "Gene_ID")]
    counts_matrix <- as.matrix(counts_data[, sample_cols])
    mode(counts_matrix) <- "numeric"
    
    # Start PDF
    pdf(output_file, width = 12, height = 8)
    
    # 1. Library size plot
    lib_sizes <- colSums(counts_matrix)/1e6
    par(mar = c(10, 4, 4, 2))
    barplot(lib_sizes,
            main = "Library Sizes",
            ylab = "Millions of Reads",
            las = 2,
            cex.names = 0.7)
    abline(h = median(lib_sizes), col = "red", lty = 2)
    
    # 2. Expression distribution
    par(mar = c(10, 4, 4, 2))
    boxplot(log2(counts_matrix + 1),
            main = "Expression Distribution",
            ylab = "log2(counts + 1)",
            las = 2,
            cex.axis = 0.7)
    
    # 3. Sample correlation heatmap
    tryCatch({
        log_counts <- log2(counts_matrix + 1)
        cor_matrix <- cor(log_counts)
        pheatmap(cor_matrix,
                main = "Sample Correlations",
                display_numbers = TRUE,
                number_format = "%.2f",
                fontsize_number = 7)
    }, error = function(e) {
        message("Could not generate correlation heatmap: ", e$message)
    })
    
    # 4. PCA plot
    if(ncol(counts_matrix) > 2) {
        tryCatch({
            # Scale the data
            scaled_data <- t(scale(t(log2(counts_matrix + 1))))
            # Remove rows with NA values
            scaled_data <- scaled_data[complete.cases(scaled_data),]
            
            # Perform PCA
            pca_result <- prcomp(t(scaled_data))
            var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)
            
            # Plot
            plot(pca_result$x[,1], 
                 pca_result$x[,2],
                 main = "PCA Plot",
                 xlab = paste0("PC1 (", var_explained[1], "%)"),
                 ylab = paste0("PC2 (", var_explained[2], "%)"),
                 pch = 19)
            text(pca_result$x[,1], 
                 pca_result$x[,2],
                 labels = colnames(counts_matrix),
                 pos = 3,
                 cex = 0.8)
        }, error = function(e) {
            message("Could not generate PCA plot: ", e$message)
        })
    }
    
    # 5. Mean-variance plot
    par(mar = c(4, 4, 4, 2))
    mean_expr <- rowMeans(counts_matrix)
    var_expr <- apply(counts_matrix, 1, var)
    plot(log2(mean_expr + 1), log2(var_expr + 1),
         main = "Mean-Variance Relationship",
         xlab = "log2(mean + 1)",
         ylab = "log2(variance + 1)",
         pch = 16,
         cex = 0.5)
    
    dev.off()
}

# 5. Generate normalized counts (CPM)
get_normalized_counts <- function(counts_data) {
    message("Generating normalized counts...")
    
    sample_cols <- colnames(counts_data)[!colnames(counts_data) %in% c("Gene_Name", "Gene_Type", "Gene_ID")]
    counts_matrix <- as.matrix(counts_data[, sample_cols])
    mode(counts_matrix) <- "numeric"
    
    # Add small pseudocount to handle zeros
    counts_matrix <- counts_matrix + 1
    
    # Calculate CPM
    lib_sizes <- colSums(counts_matrix)
    cpm_matrix <- t(t(counts_matrix) / lib_sizes) * 1e6
    
    normalized_counts <- data.frame(
        Gene_Name = counts_data$Gene_Name,
        Gene_Type = counts_data$Gene_Type,
        cpm_matrix,
        check.names = FALSE
    )
    
    return(normalized_counts)
}

# 6. Generate comprehensive statistics
generate_comprehensive_stats <- function(counts_data) {
    message("Generating comprehensive statistics...")
    
    sample_cols <- colnames(counts_data)[!colnames(counts_data) %in% c("Gene_Name", "Gene_Type")]
    counts_matrix <- as.matrix(counts_data[, sample_cols])
    mode(counts_matrix) <- "numeric"
    
    # Sample Statistics
    sample_stats <- data.frame(
        Sample = sample_cols,
        Total_Counts = colSums(counts_matrix),
        Mean_Counts = colMeans(counts_matrix),
        Median_Counts = apply(counts_matrix, 2, median),
        SD_Counts = apply(counts_matrix, 2, sd),
        Detected_Genes = colSums(counts_matrix > 0),
        Zeros_Percent = (colSums(counts_matrix == 0) / nrow(counts_matrix)) * 100
    )
    
    # Add quantiles
    sample_stats$Q1 <- apply(counts_matrix, 2, quantile, probs = 0.25)
    sample_stats$Q3 <- apply(counts_matrix, 2, quantile, probs = 0.75)
    
    # Gene Statistics
    gene_stats <- data.frame(
        Gene_Name = counts_data$Gene_Name,
        Gene_Type = counts_data$Gene_Type,
        Mean_Expression = rowMeans(counts_matrix),
        SD_Expression = apply(counts_matrix, 1, sd),
        CV = apply(counts_matrix, 1, function(x) sd(x)/mean(x)),
        Zeros_Percent = rowSums(counts_matrix == 0) / ncol(counts_matrix) * 100,
        Samples_Detected = rowSums(counts_matrix > 0)
    )
    
    return(list(
        sample_stats = sample_stats,
        gene_stats = gene_stats
    ))
}

# 7. Generate README and documentation
generate_detailed_readme <- function(results, input_file) {
    message("Generating documentation...")
    
    # Calculate basic statistics
    sample_cols <- colnames(results)[!colnames(results) %in% c("Gene_Name", "Gene_Type")]
    counts_matrix <- as.matrix(results[, sample_cols])
    mode(counts_matrix) <- "numeric"
    
    readme_text <- c(
        "RNA-seq Data Processing and Analysis Pipeline",
        "===========================================",
        "",
        sprintf("Analysis Date: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        sprintf("Input File: %s", input_file),
        "",
        "Pipeline Overview",
        "----------------",
        "This pipeline processes RNA-seq count data through multiple steps of quality control,",
        "normalization, and analysis. It handles duplicate genes, performs expression-based",
        "filtering, and generates comprehensive quality metrics.",
        "",
        "Data Processing Steps",
        "--------------------",
        "1. Initial Data Processing:",
        sprintf("   - Input genes: %d", nrow(counts_matrix)),
        sprintf("   - Number of samples: %d", ncol(counts_matrix)),
        "",
        "2. Quality Control:",
        "   - Removal of lowly expressed genes (CPM > 1 in ≥3 samples)",
        "   - Handling of duplicate gene entries",
        "   - Library size normalization",
        "",
        "3. Normalization:",
        "   - CPM (Counts Per Million) normalization",
        "   - Addition of small pseudocount for zero handling",
        "",
        "Output Files",
        "------------",
        "1. Count Matrices:",
        "   - raw_counts.csv: Filtered and cleaned raw counts",
        "   - normalized_counts.csv: CPM-normalized counts",
        "",
        "2. Quality Control:",
        "   - sample_statistics.csv: Per-sample metrics",
        "   - gene_statistics.csv: Per-gene metrics",
        "   - analysis_plots.pdf: Visualization of key metrics",
        "",
        "3. Visualizations:",
        "   - Library size distribution",
        "   - Expression boxplots",
        "   - Sample correlation heatmap",
        "   - PCA plot",
        "   - Mean-variance relationship",
        "",
        "Quality Metrics",
        "--------------",
        sprintf("1. Library Sizes:"),
        sprintf("   - Median: %s", format(median(colSums(counts_matrix)), big.mark=",")),
        sprintf("   - Mean: %s", format(mean(colSums(counts_matrix)), big.mark=",")),
        "",
        sprintf("2. Gene Detection:"),
        sprintf("   - Mean genes detected per sample: %.1f", mean(colSums(counts_matrix > 0))),
        sprintf("   - Median genes detected per sample: %.1f", median(colSums(counts_matrix > 0))),
        "",
        "Usage Notes",
        "-----------",
        "1. For differential expression analysis:",
        "   - Use raw_counts.csv",
        "   - Suitable for DESeq2 or edgeR",
        "",
        "2. For visualization and clustering:",
        "   - Use normalized_counts.csv",
        "   - Data is CPM normalized",
        "",
        "3. Quality assessment:",
        "   - Review analysis_plots.pdf",
        "   - Check sample correlations",
        "   - Examine PCA plot for batch effects",
        "",
        "Methods Details",
        "---------------",
        "1. Expression Filtering:",
        "   - Criterion: CPM > 1 in at least 3 samples",
        "   - Purpose: Remove unreliably detected genes",
        "",
        "2. Duplicate Resolution:",
        "   - Primary: Keep protein-coding version",
        "   - Secondary: Highest expressed version",
        "",
        "3. Normalization Method:",
        "   - CPM normalization",
        "   - Formula: (count / library_size) * 1e6",
        "",
        "References",
        "----------",
        "1. Robinson MD, Oshlack A (2010)",
        "   'A scaling normalization method for differential expression analysis of RNA-seq data'",
        "   Genome Biology 11:R25",
        "",
        "Session Information",
        "------------------"
    )
    
    # Add session info
    session_info <- capture.output(sessionInfo())
    readme_text <- c(readme_text, session_info)
    
    # Write README
    writeLines(readme_text, "clean_counts_analysis/README.md")
    
    # Create quick reference
    quick_ref <- c(
        "Quick Reference Guide",
        "===================",
        "",
        "Key Files:",
        "1. raw_counts.csv - Filtered count matrix",
        "2. normalized_counts.csv - CPM normalized data",
        "3. analysis_plots.pdf - QC visualizations",
        "",
        "Important Notes:",
        "- Filtered for reliable expression",
        "- Duplicate genes resolved",
        "- CPM normalization applied",
        "",
        sprintf("Last updated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
    )
    
    writeLines(quick_ref, "clean_counts_analysis/QUICK_REFERENCE.md")
}

# Function to ensure directories exist
setup_directories <- function() {
    # Create main output directory
    dirs <- c(
        "clean_counts_analysis",
        "clean_counts_analysis/raw_counts",
        "clean_counts_analysis/normalized_counts",
        "clean_counts_analysis/statistics",
        "clean_counts_analysis/plots",
        "clean_counts_analysis/metadata",
        "duplicate_genes_qc"
    )
    
    # Create all directories
    for (dir in dirs) {
        dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
}

# Modified main execution function
run_complete_pipeline <- function(input_file) {
    message("Starting RNA-seq analysis pipeline...")
    
    # Create directory structure
    setup_directories()
    
    # 1. Read and prepare data
    count_data <- prepare_count_data(input_file)
    
    # 2. Create filtered DGEList
    filtered_dge <- create_filtered_dge(count_data)
    
    # 3. Process duplicates
    combined_data <- data.frame(
        Gene_Name = filtered_dge$genes$Gene_Name,
        Gene_Type = filtered_dge$genes$Gene_Type,
        filtered_dge$counts
    )
    rownames(combined_data) <- rownames(filtered_dge)
    
    duplicate_results <- process_duplicates(combined_data)
    resolved_counts <- duplicate_results$apply_resolution(combined_data)
    
    # 4. Generate normalized counts and save results
    message("Generating and saving normalized counts...")
    normalized_counts <- get_normalized_counts(resolved_counts)
    
    # Save raw counts
    write.csv(resolved_counts, 
              file.path("clean_counts_analysis", "raw_counts", "clean_counts.csv"))
    
    # Save normalized counts
    write.csv(normalized_counts, 
              file.path("clean_counts_analysis", "normalized_counts", "normalized_counts.csv"))
    
    # 5. Generate statistics
    message("Generating statistics...")
    stats <- generate_comprehensive_stats(resolved_counts)
    write.csv(stats$sample_stats, 
              file.path("clean_counts_analysis", "statistics", "sample_statistics.csv"))
    write.csv(stats$gene_stats, 
              file.path("clean_counts_analysis", "statistics", "gene_statistics.csv"))
    
    # 6. Create visualizations
    message("Creating visualizations...")
    create_visualizations(resolved_counts, 
                         file.path("clean_counts_analysis", "plots", "analysis_plots.pdf"))
    
    # 7. Generate documentation
    message("Generating documentation...")
    generate_detailed_readme(resolved_counts, input_file)
    
    # Save session info
    capture.output(sessionInfo(), 
                  file = file.path("clean_counts_analysis", "metadata", "session_info.txt"))
    
    message("\nProcessing complete! Results are available in the 'clean_counts_analysis' directory")
    message("Please check README.md for detailed information about the analysis")
    
    return(list(
        raw_counts = resolved_counts,
        normalized_counts = normalized_counts,
        statistics = stats
    ))
}

# Example usage
# Set your input file path
# input_file <- "/path/to/your/counts.csv"
# Run the pipeline
# results <- run_complete_pipeline(input_file)

####

# Source the updated script
#source("/Users/akhaliq/Desktop/asif/bulk_counts/bulk_rnaseq_pipeline.R")

# Set your input file path
#input_file <- "/Users/akhaliq/Desktop/asif/bulk_counts/run_20241103_093957/counts/counts_protein_coding.csv"

# Run the pipeline
#results <- run_complete_pipeline(input_file)

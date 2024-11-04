# Comprehensive RNA-seq QC function
perform_comprehensive_qc <- function(count_matrix) {
    library(edgeR)
    library(ggplot2)
    library(reshape2)
    library(pheatmap)
    
    # Create QC output directory
    dir.create("comprehensive_qc", showWarnings = FALSE)
    dir.create("comprehensive_qc/plots", showWarnings = FALSE)
    
    # Get sample columns
    sample_cols <- grep("^M", colnames(count_matrix), value = TRUE)
    counts <- count_matrix[, sample_cols]
    
    # Create DGEList object
    dge <- DGEList(counts = counts)
    
    # Calculate CPM
    cpm <- cpm(dge)
    log_cpm <- cpm(dge, log=TRUE)
    
    # 1. Basic Statistics
    basic_stats <- data.frame(
        Sample = colnames(counts),
        Total_Reads = colSums(counts),
        Genes_Detected = colSums(counts > 0),
        Genes_Over_1CPM = colSums(cpm > 1),
        Median_Counts = apply(counts, 2, median),
        Mean_Counts = colMeans(counts)
    )
    
    # Add percentage metrics
    basic_stats$Percent_Genes_Detected <- (basic_stats$Genes_Detected / nrow(counts)) * 100
    basic_stats$Percent_Genes_Over_1CPM <- (basic_stats$Genes_Over_1CPM / nrow(counts)) * 100
    
    # 2. Filtering low expressed genes
    keep <- rowSums(cpm > 1) >= 3
    dge_filtered <- dge[keep, ]
    counts_filtered <- counts[keep, ]
    
    # 3. Sample correlation analysis
    cor_matrix <- cor(log_cpm)
    
    # 4. Create visualizations
    # Library size plot
    pdf("comprehensive_qc/plots/1_library_sizes.pdf", width=10, height=6)
    par(mar=c(8, 4, 4, 2))
    barplot(basic_stats$Total_Reads/1e6, 
            names.arg=basic_stats$Sample,
            main="Library Sizes",
            ylab="Millions of Reads",
            las=2,
            cex.names=0.8)
    abline(h=median(basic_stats$Total_Reads/1e6), col="red", lty=2)
    dev.off()
    
    # Genes detected plot
    pdf("comprehensive_qc/plots/2_genes_detected.pdf", width=10, height=6)
    par(mar=c(8, 4, 4, 2))
    barplot(basic_stats$Genes_Detected,
            names.arg=basic_stats$Sample,
            main="Number of Genes Detected",
            ylab="Number of Genes",
            las=2,
            cex.names=0.8)
    abline(h=median(basic_stats$Genes_Detected), col="red", lty=2)
    dev.off()
    
    # Expression distribution plot
    pdf("comprehensive_qc/plots/3_expression_distribution.pdf", width=10, height=6)
    boxplot(log_cpm, 
            main="Expression Distribution (log CPM)",
            ylab="log2 CPM",
            las=2,
            cex.axis=0.7)
    dev.off()
    
    # Sample correlation heatmap
    pdf("comprehensive_qc/plots/4_sample_correlation.pdf", width=10, height=10)
    pheatmap(cor_matrix,
             main="Sample Correlation",
             display_numbers=TRUE,
             number_format="%.2f",
             fontsize_number=7)
    dev.off()
    
    # MDS plot
    pdf("comprehensive_qc/plots/5_mds_plot.pdf", width=8, height=8)
    plotMDS(dge_filtered, main="MDS Plot")
    dev.off()
    
    # MA plots for each sample
    pdf("comprehensive_qc/plots/6_MA_plots.pdf", width=12, height=8)
    par(mfrow=c(2,2))
    for(i in seq_along(sample_cols)) {
        A <- rowMeans(log_cpm)
        M <- log_cpm[,i] - A
        smoothScatter(A, M, 
                     main=paste("MA Plot -", sample_cols[i]),
                     xlab="Average log CPM",
                     ylab="M")
        abline(h=0, col="red")
    }
    dev.off()
    
    # 5. Gene statistics
    gene_stats <- data.frame(
        Gene_Name = count_matrix$Gene_Name,
        Gene_Type = count_matrix$Gene_Type,
        Mean_Expression = rowMeans(counts),
        SD_Expression = apply(counts, 1, sd),
        CV = apply(counts, 1, function(x) sd(x)/mean(x)),
        Zeros_Percent = rowSums(counts == 0) / ncol(counts) * 100,
        Samples_Detected = rowSums(counts > 0)
    )
    
    # 6. Generate comprehensive report
    sink("comprehensive_qc/QC_report.txt")
    
    cat("RNA-seq Quality Control Report\n")
    cat("============================\n\n")
    cat("Analysis Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
    
    cat("1. Basic Statistics\n")
    cat("------------------\n")
    cat("Total samples:", ncol(counts), "\n")
    cat("Total genes:", nrow(counts), "\n")
    cat("Genes passing filter (>1 CPM in ≥3 samples):", sum(keep), "\n\n")
    
    cat("Library Size Summary:\n")
    print(summary(basic_stats$Total_Reads))
    cat("\n")
    
    cat("Genes Detected Summary:\n")
    print(summary(basic_stats$Genes_Detected))
    cat("\n")
    
    cat("2. Sample Correlations\n")
    cat("--------------------\n")
    cat("Minimum correlation:", min(cor_matrix[upper.tri(cor_matrix)]), "\n")
    cat("Median correlation:", median(cor_matrix[upper.tri(cor_matrix)]), "\n")
    cat("Maximum correlation:", max(cor_matrix[upper.tri(cor_matrix)]), "\n\n")
    
    cat("3. Gene Type Summary\n")
    cat("------------------\n")
    print(table(gene_stats$Gene_Type))
    cat("\n")
    
    cat("4. Recommendations\n")
    cat("----------------\n")
    
    # Library size recommendations
    lib_size_cv <- sd(basic_stats$Total_Reads)/mean(basic_stats$Total_Reads)
    if(lib_size_cv > 0.5) {
        cat("WARNING: High variability in library sizes detected.\n")
    }
    
    # Correlation recommendations
    if(min(cor_matrix[upper.tri(cor_matrix)]) < 0.8) {
        cat("WARNING: Low sample correlations detected. Check for potential outliers.\n")
    }
    
    sink()
    
    # 7. Save results
    write.csv(basic_stats, "comprehensive_qc/sample_statistics.csv", row.names=FALSE)
    write.csv(gene_stats, "comprehensive_qc/gene_statistics.csv", row.names=FALSE)
    saveRDS(dge_filtered, "comprehensive_qc/filtered_dge.rds")
    
    return(list(
        basic_stats = basic_stats,
        gene_stats = gene_stats,
        correlation_matrix = cor_matrix,
        filtered_dge = dge_filtered,
        filter_criteria = keep
    ))
}

      # Run comprehensive QC
qc_results <- perform_comprehensive_qc(global_counts)

# Review results
file.show("comprehensive_qc/QC_report.txt")

# Check basic statistics
head(qc_results$basic_stats)

# Look at filtered object for downstream analysis
filtered_dge <- qc_results$filtered_dge

      

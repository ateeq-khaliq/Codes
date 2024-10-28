#!/bin/bash
#SBATCH --mail-user=akhaliq@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=bam_count
#SBATCH --error=bam_count_%j.error
#SBATCH --output=bam_count_%j.out
#SBATCH --time=40:00:00
#SBATCH --mem=500G
#SBATCH --partition=gpu
#SBATCH --account=r00583

# Set up environment
echo "Loading modules..."
module load miniconda
source activate spatial
set -e
set -o pipefail

# Verify GTF file exists
GTF_FILE="/N/project/akhaliq/Ateeq_dwd/reference/gencode.v44.annotation.gtf"
if [ ! -f "$GTF_FILE" ]; then
    echo "Error: GTF file not found at $GTF_FILE"
    exit 1
fi

# Create output directory with timestamp in the specified location
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
OUTPUT_DIR="/N/project/akhaliq/Ateeq_dwd/bulk_counts_${TIMESTAMP}"
mkdir -p ${OUTPUT_DIR}

echo "Created output directory: ${OUTPUT_DIR}"

# Rest of the script remains the same until the R script part
cat > ${OUTPUT_DIR}/count_bam.R << 'EOL'
# Install and load required packages
install_and_load_packages <- function() {
    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager", repos='http://cran.rstudio.com/')
    
    packages <- c("GenomicAlignments", "GenomicFeatures", 
                 "Rsamtools", "rtracklayer", "DESeq2")
    
    for(pkg in packages) {
        if (!require(pkg, character.only = TRUE)) {
            BiocManager::install(pkg)
            library(pkg, character.only = TRUE)
        }
    }
}

install_and_load_packages()

# Function to log messages
log_message <- function(message) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    message_with_timestamp <- sprintf("[%s] %s\n", timestamp, message)
    cat(message_with_timestamp)
    # Also write to log file
    write(message_with_timestamp, file = "bam_counting.log", append = TRUE)
}

# Function to find RNA-seq BAM files
find_rnaseq_bams <- function(base_dir) {
    log_message("Searching for RNA-seq BAM files...")
    
    sample_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
    bam_files <- character()
    sample_names <- character()
    
    for(dir in sample_dirs) {
        # Look for RNA-seq subdirectory
        rna_dirs <- list.files(dir, pattern = "RS\\.v2-RNA", full.names = TRUE)
        
        if(length(rna_dirs) > 0) {
            bam_file <- list.files(rna_dirs[1], 
                                 pattern = "_T_sorted\\.bam$", 
                                 full.names = TRUE)
            
            if(length(bam_file) > 0) {
                # Verify BAM index exists
                if(!file.exists(paste0(bam_file[1], ".bai"))) {
                    log_message(sprintf("Warning: Index missing for %s", bam_file[1]))
                    next
                }
                
                bam_files <- c(bam_files, bam_file[1])
                sample_name <- basename(dirname(dirname(bam_file[1])))
                sample_names <- c(sample_names, sample_name)
                log_message(sprintf("Found BAM file for sample %s", sample_name))
            }
        }
    }
    
    names(bam_files) <- sample_names
    return(bam_files)
}

create_count_matrix <- function(bam_files, gtf_file, output_prefix) {
    log_message(sprintf("Processing %d BAM files", length(bam_files)))
    
    log_message("Creating transcript database from GTF...")
    txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")
    
    log_message("Extracting exons by gene...")
    exons_by_gene <- exonsBy(txdb, by = "gene")
    
    log_message("Creating BAM file list...")
    bfl <- BamFileList(bam_files)
    
    log_message("Counting reads across all samples...")
    se <- summarizeOverlaps(
        features = exons_by_gene,
        reads = bfl,
        mode = "Union",
        singleEnd = TRUE,
        ignore.strand = TRUE,
        BPPARAM = MulticoreParam(workers = 16)
    )
    
    # Extract count matrix
    count_matrix <- assay(se)
    
    # Calculate QC metrics
    log_message("Calculating QC metrics...")
    total_reads <- colSums(count_matrix)
    detected_genes <- colSums(count_matrix > 0)
    median_counts <- apply(count_matrix, 2, median)
    
    qc_metrics <- data.frame(
        Sample = colnames(count_matrix),
        Total_Reads = total_reads,
        Detected_Genes = detected_genes,
        Detected_Genes_Ratio = detected_genes / nrow(count_matrix) * 100,
        Median_Counts = median_counts
    )
    
    # Save outputs
    log_message("Saving results...")
    
    # Raw counts
    write.csv(count_matrix, file = file.path(output_prefix, "raw_counts.csv"))
    saveRDS(count_matrix, file = file.path(output_prefix, "raw_counts.rds"))
    
    # Normalized counts (DESeq2)
    log_message("Performing DESeq2 normalization...")
    dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                                 colData = data.frame(row.names = colnames(count_matrix),
                                                    condition = factor(rep("A", ncol(count_matrix)))),
                                 design = ~ 1)
    dds <- estimateSizeFactors(dds)
    normalized_counts <- counts(dds, normalized = TRUE)
    write.csv(normalized_counts, file = file.path(output_prefix, "normalized_counts.csv"))
    
    # QC metrics
    write.csv(qc_metrics, file = file.path(output_prefix, "qc_metrics.csv"), 
              row.names = FALSE)
    
    # Generate QC plots
    log_message("Generating QC plots...")
    pdf(file.path(output_prefix, "qc_plots.pdf"))
    
    # Library sizes
    barplot(total_reads/1e6, 
            main="Library Sizes", 
            las=2, 
            cex.names=0.7, 
            ylab="Million Reads")
    
    # Detected genes
    barplot(detected_genes, 
            main="Detected Genes", 
            las=2, 
            cex.names=0.7, 
            ylab="Number of Genes")
    
    # Sample correlations
    heatmap(cor(count_matrix), 
            main="Sample Correlations",
            cex.main=0.8)
    
    # MA plot
    plotMA(dds, main="MA Plot")
    
    dev.off()
    
    return(list(counts = count_matrix, 
                normalized_counts = normalized_counts,
                qc = qc_metrics))
}

# Main execution
tryCatch({
    # Set working directory and output directory
    project_dir <- "/N/project/akhaliq/Ateeq_dwd"
    output_dir <- Sys.getenv("OUTPUT_DIR")  # Get the output directory from environment variable
    setwd(project_dir)
    log_message(sprintf("Working directory set to: %s", project_dir))
    log_message(sprintf("Output directory set to: %s", output_dir))
    
    gtf_file <- "/N/project/akhaliq/reference/gencode.v44.annotation.gtf"
    log_message(sprintf("Using GTF file: %s", gtf_file))
    
    # Find BAM files
    bam_files <- find_rnaseq_bams(project_dir)
    log_message(sprintf("Found %d BAM files", length(bam_files)))
    
    # Generate count matrix
    results <- create_count_matrix(
        bam_files = bam_files,
        gtf_file = gtf_file,
        output_prefix = output_dir
    )
    
    # Print summary
    log_message("Processing complete!")
    log_message(sprintf("Processed %d samples", ncol(results$counts)))
    log_message(sprintf("Quantified %d genes", nrow(results$counts)))
    
}, error = function(e) {
    log_message(sprintf("Error occurred: %s", e$message))
    quit(status = 1)
})
EOL

# Pass the output directory to the R script
export OUTPUT_DIR

# Run the R script
log_file="${OUTPUT_DIR}/bam_count.log"
Rscript ${OUTPUT_DIR}/count_bam.R 2>&1 | tee ${log_file}

# Check if the script completed successfully
if [ $? -eq 0 ]; then
    echo "BAM counting completed successfully. Check results in ${OUTPUT_DIR}"
    echo "Output files:"
    ls -l ${OUTPUT_DIR}
else
    echo "Error occurred during BAM counting. Check log file: ${log_file}"
    exit 1
fi

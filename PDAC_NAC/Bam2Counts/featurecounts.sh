#!/bin/bash
#SBATCH --mail-user=akhaliq@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=feat_count
#SBATCH --error=feat_count_%j.error
#SBATCH --output=feat_count_%j.out
#SBATCH --time=40:00:00
#SBATCH --mem=500G
#SBATCH --partition=gpu
#SBATCH --account=r00583

# Set up environment
echo "Loading modules..."
module load miniconda
module load samtools
source activate spatial
set -e
set -o pipefail

# Create main output directory
MAIN_DIR="/N/project/akhaliq/Ateeq_dwd/featureCounts"
mkdir -p ${MAIN_DIR}

# Create timestamped subdirectory for this run
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
OUTPUT_DIR="${MAIN_DIR}/run_${TIMESTAMP}"
mkdir -p ${OUTPUT_DIR}/{counts,qc,logs,summary}

echo "Created output directories:"
echo "- Main output: ${OUTPUT_DIR}"
echo "- Counts directory: ${OUTPUT_DIR}/counts"
echo "- QC directory: ${OUTPUT_DIR}/qc"
echo "- Logs directory: ${OUTPUT_DIR}/logs"
echo "- Summary directory: ${OUTPUT_DIR}/summary"

# Create the R script
cat > ${OUTPUT_DIR}/logs/run_featurecounts.R << 'EOL'
# Install and load required packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos='http://cran.rstudio.com/')

packages <- c("Rsubread", "DESeq2", "rtracklayer", 
             "org.Hs.eg.db", "AnnotationDbi", "GenomicRanges",
             "Rsamtools", "GenomicFeatures")

for(pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
        BiocManager::install(pkg, update=FALSE, ask=FALSE)
        library(pkg, character.only = TRUE)
    }
}

# Function to log messages
log_message <- function(message) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    msg <- sprintf("[%s] %s", timestamp, message)
    cat(msg, "\n")
    write(msg, file = file.path(output_dir, "logs", "featurecounts.log"), append = TRUE)
}

# Function to check BAM chromosomes
check_bam_chromosomes <- function(bam_file) {
    log_message(sprintf("Checking chromosome names in BAM file: %s", basename(bam_file)))
    header <- scanBamHeader(bam_file)
    chroms <- names(header[[1]]$targets)
    log_message(sprintf("Number of chromosomes: %d", length(chroms)))
    log_message("First few chromosome names:")
    print(head(chroms))
    return(chroms)
}

# Function to find RNA-seq BAM files
find_rnaseq_bams <- function(base_dir) {
    log_message("Searching for RNA-seq BAM files...")
    
    sample_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
    bam_files <- character()
    sample_names <- character()
    
    for(dir in sort(sample_dirs)) {
        if(grepl("^(M|MV)", basename(dir))) {
            rna_dirs <- list.files(dir, pattern = "RS\\.v2-RNA", full.names = TRUE)
            
            if(length(rna_dirs) > 0) {
                bam_file <- list.files(rna_dirs[1], 
                                     pattern = "_T_sorted\\.bam$", 
                                     full.names = TRUE)
                
                if(length(bam_file) > 0 && file.exists(bam_file[1])) {
                    # Check chromosome names
                    chroms <- check_bam_chromosomes(bam_file[1])
                    
                    bam_files <- c(bam_files, bam_file[1])
                    sample_name <- basename(dir)
                    sample_names <- c(sample_names, sample_name)
                    log_message(sprintf("Found BAM file for sample %s", sample_name))
                }
            }
        }
    }
    
    names(bam_files) <- sample_names
    return(bam_files)
}

# Function to run featureCounts
run_featurecounts <- function(bam_files, gtf_file, output_dir) {
    log_message("Starting featureCounts analysis...")
    
    # Run featureCounts with GTF
    log_message("Running featureCounts...")
    
	   fc <- featureCounts(
        files = as.character(bam_files),        # Convert to character vector
        GTF.featureType = "exon",
        GTF.attrType = "gene_id",
        annot.inbuilt = NULL,                   # Set to NULL
        annot.ext = gtf_file,                   # Use annot.ext instead of annot
        isGTFAnnotationFile = TRUE,
        useMetaFeatures = TRUE,
        allowMultiOverlap = TRUE,
        isPairedEnd = FALSE,
        nthreads = 16,
        strandSpecific = 0,
        verbose = TRUE
    )
     
    # Get count matrix
    count_matrix <- fc$counts
    colnames(count_matrix) <- names(bam_files)
    
    # Extract gene names from GTF
    log_message("Processing gene names from GTF...")
    gtf <- import(gtf_file)
    gene_id_to_name <- unique(data.frame(
        gene_id = gtf$gene_id,
        gene_name = gtf$gene_name
    ))
    gene_id_to_name <- gene_id_to_name[!duplicated(gene_id_to_name$gene_id),]
    
    # Match gene names to count matrix
    gene_names <- gene_id_to_name$gene_name[match(rownames(count_matrix), 
                                                 gene_id_to_name$gene_id)]
    
    # Save raw counts
    write.csv(count_matrix, 
              file = file.path(output_dir, "counts", "raw_counts_ensembl.csv"))
    
    # Save raw counts with gene names
    count_matrix_named <- count_matrix
    rownames(count_matrix_named) <- ifelse(is.na(gene_names), 
                                         rownames(count_matrix), 
                                         gene_names)
    write.csv(count_matrix_named, 
              file = file.path(output_dir, "counts", "raw_counts_gene_names.csv"))
    
    # Save featureCounts stats
    write.csv(fc$stat, 
              file = file.path(output_dir, "qc", "featurecounts_stats.csv"))
    
    # Calculate and save QC metrics
    qc_metrics <- data.frame(
        Sample = colnames(count_matrix),
        Total_Reads = colSums(count_matrix),
        Detected_Genes = colSums(count_matrix > 0),
        Percent_Detected = (colSums(count_matrix > 0) / nrow(count_matrix)) * 100,
        Mean_Counts = colMeans(count_matrix),
        Median_Counts = apply(count_matrix, 2, median)
    )
    write.csv(qc_metrics, 
              file = file.path(output_dir, "qc", "qc_metrics.csv"),
              row.names = FALSE)
    
    # DESeq2 normalization
    log_message("Performing DESeq2 normalization...")
    dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                                 colData = data.frame(row.names = colnames(count_matrix)),
                                 design = ~1)
    
    dds <- estimateSizeFactors(dds)
    normalized_counts <- counts(dds, normalized = TRUE)
    
    # Save normalized counts
    write.csv(normalized_counts, 
              file = file.path(output_dir, "counts", "normalized_counts_ensembl.csv"))
    
    normalized_counts_named <- normalized_counts
    rownames(normalized_counts_named) <- ifelse(is.na(gene_names), 
                                              rownames(normalized_counts), 
                                              gene_names)
    write.csv(normalized_counts_named, 
              file = file.path(output_dir, "counts", "normalized_counts_gene_names.csv"))
    
    # Create summary report
    sink(file.path(output_dir, "summary", "analysis_summary.txt"))
    cat("featureCounts Analysis Summary\n")
    cat("============================\n\n")
    cat(sprintf("Date: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    cat(sprintf("Total samples processed: %d\n", ncol(count_matrix)))
    cat(sprintf("Total genes quantified: %d\n", nrow(count_matrix)))
    cat("\nSample Summary:\n")
    print(qc_metrics)
    cat("\nFeatureCounts Statistics:\n")
    print(fc$stat)
    sink()
    
    # Create sample list
    write.table(data.frame(Sample = names(bam_files), 
                          BAM_File = bam_files),
                file = file.path(output_dir, "summary", "processed_samples.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    return(list(counts = count_matrix,
                normalized = normalized_counts,
                stats = fc$stat,
                qc = qc_metrics))
}

# Main execution
tryCatch({
    project_dir <- "/N/project/akhaliq/Ateeq_dwd"
    output_dir <- Sys.getenv("OUTPUT_DIR")
    setwd(project_dir)
    
    gtf_file <- file.path(project_dir, "reference/gencode.v44.annotation.gtf")
    
    # Find BAM files
    bam_files <- find_rnaseq_bams(project_dir)
    log_message(sprintf("Found %d BAM files", length(bam_files)))
    
    # Run featureCounts
    results <- run_featurecounts(
        bam_files = bam_files,
        gtf_file = gtf_file,
        output_dir = output_dir
    )
    
    log_message("Processing complete!")
    
}, error = function(e) {
    log_message(sprintf("Error occurred: %s", e$message))
    quit(status = 1)
})
EOL

# Export output directory
export OUTPUT_DIR

# Run the R script
log_file="${OUTPUT_DIR}/logs/featurecounts.log"
Rscript ${OUTPUT_DIR}/logs/run_featurecounts.R 2>&1 | tee ${log_file}

# Check completion and create summary links
if [ $? -eq 0 ]; then
    echo "featureCounts analysis completed successfully!"
    
    # Create symbolic links to latest results
    ln -sfn ${OUTPUT_DIR} ${MAIN_DIR}/latest
    
    # Print summary
    echo -e "\nOutput files are organized as follows:"
    echo "${OUTPUT_DIR}/"
    echo "├── counts/               # Raw and normalized count matrices"
    echo "│   ├── raw_counts_ensembl.csv"
    echo "│   ├── raw_counts_gene_names.csv"
    echo "│   ├── normalized_counts_ensembl.csv"
    echo "│   └── normalized_counts_gene_names.csv"
    echo "├── qc/                   # Quality control metrics"
    echo "│   ├── featurecounts_stats.csv"
    echo "│   └── qc_metrics.csv"
    echo "├── logs/                 # Log files"
    echo "│   ├── featurecounts.log"
    echo "│   └── run_featurecounts.R"
    echo "└── summary/              # Analysis summaries"
    echo "    ├── analysis_summary.txt"
    echo "    └── processed_samples.txt"
    echo -e "\nLatest results linked to: ${MAIN_DIR}/latest"
    
else
    echo "Error occurred during analysis. Check log file: ${log_file}"
    exit 1
fi

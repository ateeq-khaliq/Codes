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

# Function to process GTF and extract gene information
process_gene_info <- function(gtf_file) {
    log_message("Processing GTF file for gene information...")
    
    # Import GTF
    gtf <- import(gtf_file)
    
    # Extract gene-level information only
    genes <- gtf[gtf$type == "gene"]
    
    # Create gene info dataframe
    gene_info <- data.frame(
        gene_id = genes$gene_id,
        gene_name = genes$gene_name,
        gene_type = genes$gene_type,
        stringsAsFactors = FALSE
    )
    
    # Remove duplicates
    gene_info <- gene_info[!duplicated(gene_info$gene_id),]
    
    # Filter out non-standard gene names
    gene_info$gene_name <- ifelse(
        !grepl("^ENS", gene_info$gene_name) & 
        !is.na(gene_info$gene_name) & 
        gene_info$gene_name != "",
        gene_info$gene_name,
        gene_info$gene_id
    )
    
    return(gene_info)
}

# Function to run featureCounts
run_featurecounts <- function(bam_files, gtf_file, output_dir) {
    log_message("Starting featureCounts analysis...")
    
    # Get gene information first
    gene_info <- process_gene_info(gtf_file)
    
    # Run featureCounts
    log_message("Running featureCounts...")
    fc <- featureCounts(
        files = as.character(bam_files),
        annot.ext = gtf_file,
        isGTFAnnotationFile = TRUE,
        GTF.featureType = "exon",
        GTF.attrType = "gene_id",
        useMetaFeatures = TRUE,
        isPairedEnd = FALSE,
        nthreads = 16,
        countMultiMappingReads = FALSE,
        strandSpecific = 0,
        verbose = TRUE
    )
    
    # Get count matrix
    count_matrix <- fc$counts
    colnames(count_matrix) <- names(bam_files)
    
    # Match gene names using the processed gene info
    gene_names <- gene_info$gene_name[match(rownames(count_matrix), gene_info$gene_id)]
    gene_types <- gene_info$gene_type[match(rownames(count_matrix), gene_info$gene_id)]
    
    # Create annotated count matrix
    count_matrix_annotated <- data.frame(
        Gene_ID = rownames(count_matrix),
        Gene_Name = gene_names,
        Gene_Type = gene_types,
        count_matrix,
        stringsAsFactors = FALSE
    )
    
    # Save different versions of the count matrix
    # 1. Full annotated version
    write.csv(count_matrix_annotated, 
              file = file.path(output_dir, "counts", "counts_full_annotation.csv"),
              row.names = FALSE)
    
    # 2. Raw counts with gene names
    count_matrix_named <- count_matrix
    rownames(count_matrix_named) <- gene_names
    write.csv(count_matrix_named, 
              file = file.path(output_dir, "counts", "raw_counts_gene_names.csv"))
    
    # 3. Filtered version (protein coding genes only)
    protein_coding <- count_matrix_annotated[count_matrix_annotated$Gene_Type == "protein_coding",]
    write.csv(protein_coding, 
              file = file.path(output_dir, "counts", "counts_protein_coding.csv"),
              row.names = FALSE)
    
    # Save feature counts statistics
    write.csv(fc$stat, 
              file = file.path(output_dir, "qc", "featurecounts_stats.csv"))
    
    # Calculate QC metrics
    qc_metrics <- data.frame(
        Sample = colnames(count_matrix),
        Total_Reads = colSums(count_matrix),
        Detected_Genes = colSums(count_matrix > 0),
        Detected_Protein_Coding = colSums(count_matrix[gene_types == "protein_coding",] > 0),
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
    
    # Save normalized counts with annotation
    normalized_counts_annotated <- data.frame(
        Gene_ID = rownames(normalized_counts),
        Gene_Name = gene_names,
        Gene_Type = gene_types,
        normalized_counts,
        stringsAsFactors = FALSE
    )
    
    write.csv(normalized_counts_annotated, 
              file = file.path(output_dir, "counts", "normalized_counts_annotated.csv"),
              row.names = FALSE)
    
    # Create gene type summary
    gene_type_summary <- table(gene_types)
    gene_type_counts <- data.frame(
        Gene_Type = names(gene_type_summary),
        Count = as.numeric(gene_type_summary)
    )
    write.csv(gene_type_counts,
              file = file.path(output_dir, "summary", "gene_type_summary.csv"),
              row.names = FALSE)
    
    # Create comprehensive summary report
    sink(file.path(output_dir, "summary", "analysis_summary.txt"))
    cat("featureCounts Analysis Summary\n")
    cat("============================\n\n")
    cat(sprintf("Date: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    cat(sprintf("Total samples processed: %d\n", ncol(count_matrix)))
    cat(sprintf("Total genes quantified: %d\n", nrow(count_matrix)))
    cat(sprintf("Protein-coding genes: %d\n", sum(gene_types == "protein_coding")))
    cat("\nGene Type Summary:\n")
    print(gene_type_counts)
    cat("\nSample Summary:\n")
    print(qc_metrics)
    cat("\nFeatureCounts Statistics:\n")
    print(fc$stat)
    sink()
    
    return(list(counts = count_matrix,
                normalized = normalized_counts,
                stats = fc$stat,
                qc = qc_metrics,
                gene_info = gene_info))
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
    echo "├── counts/"
    echo "│   ├── counts_full_annotation.csv     # Complete counts with annotations"
    echo "│   ├── counts_protein_coding.csv      # Protein-coding genes only"
    echo "│   ├── raw_counts_gene_names.csv      # Raw counts with gene names"
    echo "│   └── normalized_counts_annotated.csv # Normalized counts with annotations"
    echo "├── qc/"
    echo "│   ├── featurecounts_stats.csv       # Counting statistics"
    echo "│   └── qc_metrics.csv                # Quality metrics"
    echo "├── logs/"
    echo "│   ├── featurecounts.log            # Processing log"
    echo "│   └── run_featurecounts.R          # R script"
    echo "└── summary/"
    echo "    ├── analysis_summary.txt         # Complete analysis summary"
    echo "    └── gene_type_summary.csv        # Gene type statistics"
    echo -e "\nLatest results linked to: ${MAIN_DIR}/latest"
    
else
    echo "Error occurred during analysis. Check log file: ${log_file}"
    exit 1
fi

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

# Create output directory
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
OUTPUT_DIR="/N/project/akhaliq/Ateeq_dwd/featurecounts_${TIMESTAMP}"
mkdir -p ${OUTPUT_DIR}

# Create the R script
cat > ${OUTPUT_DIR}/run_featurecounts.R << 'EOL'
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
    write(msg, file = "featurecounts.log", append = TRUE)
}

# Function to check chromosome names in BAM file
check_bam_chromosomes <- function(bam_file) {
    log_message(sprintf("Checking chromosome names in BAM file: %s", basename(bam_file)))
    header <- scanBamHeader(bam_file)
    chroms <- names(header[[1]]$targets)
    log_message(sprintf("Number of chromosomes: %d", length(chroms)))
    log_message("First few chromosome names:")
    print(head(chroms))
    return(chroms)
}

# Function to find RNA-seq BAM files and check their chromosome names
find_rnaseq_bams <- function(base_dir) {
    log_message("Searching for RNA-seq BAM files...")
    
    sample_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
    bam_files <- character()
    sample_names <- character()
    chrom_names <- list()
    
    for(dir in sort(sample_dirs)) {
        if(grepl("^(M|MV)", basename(dir))) {
            rna_dirs <- list.files(dir, pattern = "RS\\.v2-RNA", full.names = TRUE)
            
            if(length(rna_dirs) > 0) {
                bam_file <- list.files(rna_dirs[1], 
                                     pattern = "_T_sorted\\.bam$", 
                                     full.names = TRUE)
                
                if(length(bam_file) > 0 && file.exists(bam_file[1])) {
                    # Check if index exists
                    if(!file.exists(paste0(bam_file[1], ".bai"))) {
                        log_message(sprintf("Warning: Index missing for %s", bam_file[1]))
                        next
                    }
                    
                    # Check chromosome names
                    chroms <- check_bam_chromosomes(bam_file[1])
                    chrom_names[[basename(dir)]] <- chroms
                    
                    bam_files <- c(bam_files, bam_file[1])
                    sample_name <- basename(dir)
                    sample_names <- c(sample_names, sample_name)
                    log_message(sprintf("Found BAM file for sample %s", sample_name))
                }
            }
        }
    }
    
    # Check chromosome name consistency
    log_message("Checking chromosome name consistency across BAM files...")
    all_chroms <- unique(unlist(chrom_names))
    has_chr_prefix <- any(grepl("^chr", all_chroms))
    log_message(sprintf("Chromosomes have 'chr' prefix: %s", has_chr_prefix))
    
    # Save chromosome information
    chrom_info <- data.frame(
        Sample = rep(names(chrom_names), sapply(chrom_names, length)),
        Chromosome = unlist(chrom_names)
    )
    write.csv(chrom_info, "chromosome_names.csv", row.names = FALSE)
    
    names(bam_files) <- sample_names
    return(list(bam_files = bam_files, 
               chrom_names = chrom_names, 
               has_chr_prefix = has_chr_prefix))
}

# Function to process GTF and harmonize chromosome names
process_gtf <- function(gtf_file, has_chr_prefix) {
    log_message("Processing GTF file...")
    gtf <- import(gtf_file)
    
    # Check GTF chromosome names
    gtf_chroms <- unique(seqnames(gtf))
    log_message("Original GTF chromosome names:")
    print(gtf_chroms)
    
    # Harmonize chromosome names based on BAM files
    if(!has_chr_prefix) {
        log_message("Removing 'chr' prefix from GTF chromosome names...")
        seqlevels(gtf) <- sub("^chr", "", seqlevels(gtf))
    }
    
    log_message("Modified GTF chromosome names:")
    print(unique(seqnames(gtf)))
    
    # Create SAF format
    saf <- data.frame(
        GeneID = gtf$gene_id,
        Chr = seqnames(gtf),
        Start = start(gtf),
        End = end(gtf),
        Strand = strand(gtf),
        Gene_Name = gtf$gene_name
    )
    
    # Remove duplicates
    saf <- saf[!duplicated(saf$GeneID),]
    
    log_message(sprintf("Number of genes in SAF: %d", nrow(saf)))
    log_message("SAF format chromosome names:")
    print(unique(saf$Chr))
    
    return(saf)
}

# Main function to run featureCounts
run_featurecounts <- function(bam_files, gtf_file, has_chr_prefix, output_prefix) {
    log_message("Starting featureCounts analysis...")
    
    # Process GTF file
    saf <- process_gtf(gtf_file, has_chr_prefix)
    saf_file <- file.path(output_prefix, "genes.saf")
    write.table(saf, file=saf_file, sep="\t", quote=FALSE, row.names=FALSE)
    
    # Run featureCounts
    log_message("Running featureCounts...")
    fc <- featureCounts(files = bam_files,
                       annot.ext = saf_file,
                       annot.ext.format = "SAF",
                       isGTFAnnotationFile = FALSE,
                       isPairedEnd = FALSE,
                       nthreads = 16,
                       countMultiMappingReads = FALSE,
                       strandSpecific = 0,
                       verbose = TRUE)
    
    # Get count matrix
    count_matrix <- fc$counts
    colnames(count_matrix) <- names(bam_files)
    
    # Add gene names
    gene_names <- saf$Gene_Name[match(rownames(count_matrix), saf$GeneID)]
    
    # Save raw counts with Ensembl IDs
    write.csv(count_matrix, 
              file = file.path(output_prefix, "raw_counts_ensembl.csv"))
    
    # Save raw counts with gene names
    count_matrix_named <- count_matrix
    rownames(count_matrix_named) <- ifelse(is.na(gene_names), 
                                         rownames(count_matrix), 
                                         gene_names)
    write.csv(count_matrix_named, 
              file = file.path(output_prefix, "raw_counts_gene_names.csv"))
    
    # Save featureCounts stats
    write.csv(fc$stat, 
              file = file.path(output_prefix, "featurecounts_stats.csv"))
    
    # Calculate and save additional QC metrics
    qc_metrics <- data.frame(
        Sample = colnames(count_matrix),
        Total_Reads = colSums(count_matrix),
        Detected_Genes = colSums(count_matrix > 0),
        Percent_Detected = (colSums(count_matrix > 0) / nrow(count_matrix)) * 100,
        Mean_Counts = colMeans(count_matrix),
        Median_Counts = apply(count_matrix, 2, median)
    )
    write.csv(qc_metrics, 
              file = file.path(output_prefix, "qc_metrics.csv"),
              row.names = FALSE)
    
    # DESeq2 normalization
    log_message("Performing DESeq2 normalization...")
    dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                                 colData = data.frame(row.names = colnames(count_matrix)),
                                 design = ~1)
    
    dds <- estimateSizeFactors(dds)
    normalized_counts <- counts(dds, normalized = TRUE)
    
    # Save normalized counts with Ensembl IDs
    write.csv(normalized_counts, 
              file = file.path(output_prefix, "normalized_counts_ensembl.csv"))
    
    # Save normalized counts with gene names
    normalized_counts_named <- normalized_counts
    rownames(normalized_counts_named) <- ifelse(is.na(gene_names), 
                                              rownames(normalized_counts), 
                                              gene_names)
    write.csv(normalized_counts_named, 
              file = file.path(output_prefix, "normalized_counts_gene_names.csv"))
    
    # Create summary report
    sink(file.path(output_prefix, "analysis_summary.txt"))
    cat("featureCounts Analysis Summary\n")
    cat("============================\n\n")
    cat(sprintf("Total samples processed: %d\n", ncol(count_matrix)))
    cat(sprintf("Total genes quantified: %d\n", nrow(count_matrix)))
    cat("\nSample Summary:\n")
    print(qc_metrics)
    cat("\nFeatureCounts Statistics:\n")
    print(fc$stat)
    sink()
    
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
    
    # Find BAM files and check chromosome names
    bam_info <- find_rnaseq_bams(project_dir)
    log_message(sprintf("Found %d BAM files", length(bam_info$bam_files)))
    
    # Run featureCounts with chromosome name handling
    results <- run_featurecounts(
        bam_files = bam_info$bam_files,
        gtf_file = gtf_file,
        has_chr_prefix = bam_info$has_chr_prefix,
        output_prefix = output_dir
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
log_file="${OUTPUT_DIR}/featurecounts.log"
Rscript ${OUTPUT_DIR}/run_featurecounts.R 2>&1 | tee ${log_file}

# Check completion
if [ $? -eq 0 ]; then
    echo "featureCounts analysis completed successfully!"
    echo "Results are in: ${OUTPUT_DIR}"
    echo "Please check the following files:"
    echo "1. raw_counts_gene_names.csv (Raw counts with gene symbols)"
    echo "2. normalized_counts_gene_names.csv (Normalized counts)"
    echo "3. qc_metrics.csv (Quality control metrics)"
    echo "4. analysis_summary.txt (Complete analysis summary)"
    echo "5. featurecounts.log (Detailed log file)"
else
    echo "Error occurred during analysis. Check log file: ${log_file}"
    exit 1
fi

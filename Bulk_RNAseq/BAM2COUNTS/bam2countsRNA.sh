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
module load samtools
source activate spatial
set -e
set -o pipefail

# Create output directory
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
OUTPUT_DIR="/N/project/akhaliq/Ateeq_dwd/bulk_counts_${TIMESTAMP}"
mkdir -p ${OUTPUT_DIR}

# Reindex BAM files that need it
echo "Checking and updating BAM indices..."
for bam in $(find /N/project/akhaliq/Ateeq_dwd -name "*_T_sorted.bam"); do
    if [ -f "${bam}" ]; then
        if [ -f "${bam}.bai" ]; then
            if [ "${bam}" -nt "${bam}.bai" ]; then
                echo "Reindexing ${bam}..."
                samtools index -b ${bam}
            fi
        else
            echo "Creating index for ${bam}..."
            samtools index -b ${bam}
        fi
    fi
done

# Create the R script
cat > ${OUTPUT_DIR}/count_bam.R << 'EOL'
# Install and load required packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos='http://cran.rstudio.com/')

packages <- c("GenomicAlignments", "GenomicFeatures", "Rsamtools", 
              "rtracklayer", "DESeq2", "BiocParallel", 
              "GenomicRanges", "org.Hs.eg.db", "AnnotationDbi")

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
    write(msg, file = "bam_counting.log", append = TRUE)
}

# Function to find RNA-seq BAM files
find_rnaseq_bams <- function(base_dir) {
    log_message("Searching for RNA-seq BAM files...")
    
    sample_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
    bam_files <- character()
    sample_names <- character()
    
    for(dir in sort(sample_dirs)) {
        if(grepl("^(M|MV)", basename(dir))) {  # Only process M* and MV* directories
            rna_dirs <- list.files(dir, pattern = "RS\\.v2-RNA", full.names = TRUE)
            
            if(length(rna_dirs) > 0) {
                bam_file <- list.files(rna_dirs[1], 
                                     pattern = "_T_sorted\\.bam$", 
                                     full.names = TRUE)
                
                if(length(bam_file) > 0) {
                    tryCatch({
                        check_bam(bam_file[1])
                        bam_files <- c(bam_files, bam_file[1])
                        sample_name <- basename(dir)
                        sample_names <- c(sample_names, sample_name)
                        log_message(sprintf("Found valid BAM file for sample %s", sample_name))
                    }, error = function(e) {
                        log_message(sprintf("Error with BAM file for %s: %s", basename(dir), e$message))
                    })
                }
            }
        }
    }
    
    names(bam_files) <- sample_names
    return(bam_files)
}

# Function to check BAM file
check_bam <- function(bam_file) {
    log_message(sprintf("Checking BAM file: %s", basename(bam_file)))
    
    if(!file.exists(bam_file)) {
        stop(sprintf("BAM file not found: %s", bam_file))
    }
    
    bam <- scanBam(bam_file, param=ScanBamParam(what=c("qname"), tag=character(0)))[[1]]
    if(length(bam$qname) == 0) {
        stop(sprintf("No reads found in BAM file: %s", bam_file))
    }
    
    log_message(sprintf("BAM file ok, found %d reads in quick check", length(bam$qname)))
}

# Function to load GTF file with modified chromosome names and create TxDb object
create_txdb <- function(gtf_file) {
    log_message("Loading GTF file and creating TxDb object with modified chromosome names...")
    gtf <- import(gtf_file)
    seqlevels(gtf) <- sub("^chr", "", seqlevels(gtf))
    txdb <- makeTxDbFromGRanges(gtf)
    return(txdb)
}

# Main function to process BAM files
create_count_matrix <- function(bam_files, gtf_file, output_prefix) {
    log_message(sprintf("Processing %d BAM files", length(bam_files)))
    
    txdb <- create_txdb(gtf_file)
    log_message("Extracting exons by gene...")
    exons_by_gene <- exonsBy(txdb, by = "gene")
    
    param <- ScanBamParam(what = c("qname", "flag", "mapq"))
    bfl <- BamFileList(bam_files, yieldSize = 1000000)
    register(MulticoreParam(workers = 16))
    
    log_message("Counting reads across all samples...")
    se <- summarizeOverlaps(
        features = exons_by_gene,
        reads = bfl,
        mode = "Union",
        singleEnd = TRUE,
        ignore.strand = TRUE,
        param = param
    )
    
    count_matrix <- assay(se)
    write.csv(count_matrix, file = file.path(output_prefix, "raw_counts.csv"))
    
    log_message("Converting Ensembl IDs to gene names for raw counts...")
    gene_names <- mapIds(org.Hs.eg.db, keys=row.names(count_matrix), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
    rownames(count_matrix) <- ifelse(is.na(gene_names), row.names(count_matrix), gene_names)
    write.csv(count_matrix, file = file.path(output_prefix, "raw_counts_with_gene_names.csv"))
    
    log_message("Normalizing counts with DESeq2...")
    dds <- DESeqDataSetFromMatrix(countData = assay(se), colData = DataFrame(row.names = colnames(count_matrix)), design = ~1)
    dds <- DESeq(dds)
    normalized_counts <- counts(dds, normalized=TRUE)
    write.csv(normalized_counts, file = file.path(output_prefix, "normalized_counts.csv"))
    
    log_message("Converting Ensembl IDs to gene names for normalized counts...")
    rownames(normalized_counts) <- ifelse(is.na(gene_names), row.names(normalized_counts), gene_names)
    write.csv(normalized_counts, file = file.path(output_prefix, "normalized_counts_with_gene_names.csv"))
}

# Main execution
tryCatch({
    project_dir <- "/N/project/akhaliq/Ateeq_dwd"
    output_dir <- Sys.getenv("OUTPUT_DIR")
    setwd(project_dir)
    
    gtf_file <- file.path(project_dir, "reference/gencode.v44.annotation.gtf")
    bam_files <- find_rnaseq_bams(project_dir)
    
    log_message(sprintf("Found %d BAM files", length(bam_files)))
    create_count_matrix(
        bam_files = bam_files,
        gtf_file = gtf_file,
        output_prefix = output_dir
    )
}, error = function(e) {
    log_message(sprintf("Error occurred: %s", e$message))
    quit(status = 1)
})
EOL

# Export output directory
export OUTPUT_DIR

# Run the R script
log_file="${OUTPUT_DIR}/bam_count.log"
Rscript ${OUTPUT_DIR}/count_bam.R 2>&1 | tee ${log_file}

# Check completion
if [ $? -eq 0 ]; then
    echo "BAM counting completed. Check results in ${OUTPUT_DIR}"
else
    echo "Error occurred during BAM counting. Check log file: ${log_file}"
    exit 1
fi

#!/bin/bash
#SBATCH --mail-user=akhaliq@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=gsea_analysis
#SBATCH --error=gsea_analysis.error
#SBATCH --output=gsea_analysis.out
#SBATCH --time=24:00:00
#SBATCH --mem=500G
#SBATCH --partition=gpu
#SBATCH --account=r00583

# Set up environment
module load java/17.0.7
module load miniconda
source activate spatial
set -e
set -o pipefail

# Check Java version
java -version


# Set paths
BASE_DIR="/N/project/akhaliq/gsea/gsea_new"
GSEA_DIR="/N/project/akhaliq/gsea/GSEA_4.3.3"
GSEA_JAR="$GSEA_DIR/gsea-cli.sh"
GMT_DIR="/N/project/akhaliq/gsea/msigdb_v2024.1.Hs_files_to_download_locally/msigdb_v2024.1.Hs_GMTs"
CHIP_FILE="/N/project/akhaliq/gsea/msigdb_v2024.1.Hs_chip_files_to_download_locally/Human_Gene_Symbol_with_Remapping_MSigDB.v2024.1.Hs.chip"
GCT_FILE="$BASE_DIR/tumor_core_cc7_10.gct"
CLS_FILE="$BASE_DIR/tumor_core_cc7_10.cls"

# Set up output directory
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
OUTPUT_DIR="$BASE_DIR/test/gsea_results_$TIMESTAMP"
mkdir -p "$OUTPUT_DIR"

# Function to run GSEA
run_gsea() {
    local gmt_file="$1"
    local output_name="$2"
    
    echo "Starting GSEA for $output_name"
    echo "GMT file: $gmt_file"
    echo "Output directory: $OUTPUT_DIR/$output_name"
    
    # Change to GSEA directory
    cd "$GSEA_DIR"
    bash "$GSEA_JAR" GSEA \
        -gmx "$gmt_file" \
        -chip "$CHIP_FILE" \
        -res "$GCT_FILE" \
        -cls "$CLS_FILE#Treated_versus_Untreated" \
        -out "$OUTPUT_DIR/$output_name" \
        -rpt_label "$output_name" \
        -collapse false \
        -mode Max_probe \
        -norm meandiv \
        -nperm 1000 \
        -permute phenotype \
        -scoring_scheme weighted \
        -metric Signal2Noise \
        -rnd_seed 12345 \
        -set_max 500 \
        -set_min 15 \
        -zip_report false
    
    # Check if the command was successful
    if [ $? -ne 0 ]; then
        echo "Error running GSEA for $output_name" >&2
        return 1
    fi
    
    # Change back to the original directory
    cd -
    
    echo "GSEA completed for $output_name"
    echo "Contents of output directory:"
    ls -R "$OUTPUT_DIR/$output_name"
}

# Run GSEA for multiple gene set collections
collections=("h.all.v2024.1.Hs.symbols.gmt" "c2.cp.v2024.1.Hs.symbols.gmt" "c5.go.bp.v2024.1.Hs.symbols.gmt")

for collection in "${collections[@]}"; do
    gmt_file="$GMT_DIR/$collection"
    output_name="${collection%.gmt}"
    
    echo "Running GSEA for $output_name"
    run_gsea "$gmt_file" "$output_name"
    
    if [ $? -ne 0 ]; then
        echo "Error running GSEA for $output_name" >&2
        exit 1
    fi
done

# Organize results for Enrichment Map
ENRICHMENT_MAP_DIR="$OUTPUT_DIR/enrichment_map_input"
mkdir -p "$ENRICHMENT_MAP_DIR"

for collection in "${collections[@]}"; do
    output_name="${collection%.gmt}"
    result_dir="$OUTPUT_DIR/$output_name"
    
    echo "Checking results for $output_name"
    echo "Looking in directory: $result_dir"
    
    if [ -d "$result_dir" ]; then
        echo "Result directory found for $output_name"
        # Copy relevant files
        find "$result_dir" -name "gsea_report_for_Treated_*.xls" -exec cp {} "$ENRICHMENT_MAP_DIR/${output_name}_Treated.xls" \;
        find "$result_dir" -name "gsea_report_for_Untreated_*.xls" -exec cp {} "$ENRICHMENT_MAP_DIR/${output_name}_Untreated.xls" \;
        find "$result_dir" -name "ranked_gene_list*.xls" -exec cp {} "$ENRICHMENT_MAP_DIR/${output_name}_ranked_gene_list.xls" \;
        
        echo "Copied files for $output_name:"
        ls -l "$ENRICHMENT_MAP_DIR"
    else
        echo "Results directory not found for $output_name" >&2
    fi
done

echo "GSEA analysis complete. Results are in $OUTPUT_DIR"
echo "Files for Enrichment Map are in $ENRICHMENT_MAP_DIR"
echo "Contents of Enrichment Map directory:"
ls -l "$ENRICHMENT_MAP_DIR"


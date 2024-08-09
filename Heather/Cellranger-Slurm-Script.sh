#!/bin/bash
#SBATCH --mail-user=akhaliq@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=cellranger_multi
#SBATCH --error=cellranger_multi.error
#SBATCH --output=cellranger_multi.out
#SBATCH --time=1-00:00:00
#SBATCH --mem=400G
#SBATCH --partition=general
#SBATCH --account=r00583

# Change to the project directory
cd /N/project/cytassist/heather/Chrm_499_OHagan_IUB_FLEX15_Mouse_May2024

# Load necessary modules
module load cellranger/7.2.0

# Create the configuration file
echo "Creating CellRanger configuration file..."
cat << EOF > cellranger_config.csv
[gene-expression]
reference,/N/project/cytassist/heather/Chrm_499_OHagan_IUB_FLEX15_Mouse_May2024/data/external/reference/refdata-gex-mm10-2020-A
probe-set,/N/project/cytassist/heather/Chrm_499_OHagan_IUB_FLEX15_Mouse_May2024/data/external/probe_set/Chromium_Mouse_Transcriptome_Probe_Set_v1.0.1_mm10-2020-A.csv

[libraries]
fastq_id,fastqs,feature_types
Chrm_499_OHagan_IUB_FLEX15,/N/project/cytassist/heather/Chrm_499_OHagan_IUB_FLEX15_Mouse_May2024,Gene Expression

[samples]
sample_id,probe_barcode_ids,description
968,BC001,Vehicle treated female 1
966,BC002,Vehicle treated female 2
979,BC003,Vehicle treated female 3
963,BC004,Inhibitor treated female 1
975,BC005,Inhibitor treated female 2
964,BC006,Inhibitor treated female 3
960,BC007,Vehicle treated male 1
978,BC008,Vehicle treated male 2
969,BC009,Vehicle treated male 3
983,BC010,Vehicle treated male 4
0,BC011,Inhibitor treated male 1
959,BC012,Inhibitor treated male 2
965,BC013,Inhibitor treated male 3
967,BC014,Inhibitor treated male 4
977,BC015,Inhibitor treated male 5
EOF

echo "Configuration file created: cellranger_config.csv"

# Run CellRanger with the generated configuration
echo "Running CellRanger..."
cellranger multi --id=Chrm_499_OHagan_IUB_FLEX15_run6 --csv=cellranger_config.csv

echo "CellRanger multi job completed."

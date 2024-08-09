# CellRanger Analysis Script for Visium FLEX Mouse Data

## Project Overview

This project, conducted by Ateeq Khaliq for Heather's Syngenic Mouse Models, focuses on analyzing spatial transcriptomics data from mouse samples using the Visium FLEX platform. The study aims to investigate gene expression patterns in various conditions, including vehicle-treated and inhibitor-treated samples from both male and female mice.

## Data Description

- **Organism**: Mouse (Mus musculus)
- **Experiment Type**: Visium FLEX spatial transcriptomics
- **Sample Types**: 
  - Vehicle-treated female mice
  - Inhibitor-treated female mice
  - Vehicle-treated male mice
  - Inhibitor-treated male mice
- **Total Samples**: 15 (8 vehicle-treated, 7 inhibitor-treated)
- **Data Format**: FASTQ files from Illumina sequencing

## Script Overview

This SLURM script automates the process of running CellRanger multi on the Visium FLEX spatial transcriptomics data. It creates a configuration file with the correct probe barcode IDs and runs CellRanger using this configuration.

## Requirements

- Access to a SLURM-based HPC system
- CellRanger (version 7.2.0) installed and accessible as a module
- Sufficient storage space for CellRanger output (>500GB recommended)
- Input data:
  - FASTQ files from your Visium FLEX experiment
  - Mouse reference genome (mm10)
  - Visium FLEX probe set file for mouse

## Directory Structure

Ensure your project directory is structured as follows:

```
/N/project/cytassist/heather/Chrm_499_OHagan_IUB_FLEX15_Mouse_May2024/
├── cellranger-slurm-script.sh
├── Chrm_499_OHagan_IUB_FLEX15_S1_L00[7-8]_[R1,R2,I1,I2]_001.fastq.gz
└── data/
    └── external/
        ├── reference/
        │   └── refdata-gex-mm10-2020-A/
        └── probe_set/
            └── Chromium_Mouse_Transcriptome_Probe_Set_v1.0.1_mm10-2020-A.csv
```

## Usage

1. Save the script as `cellranger-slurm-script.sh` in your project directory.

2. Make the script executable:
   ```
   chmod +x cellranger-slurm-script.sh
   ```

3. Submit the job to SLURM:
   ```
   sbatch cellranger-slurm-script.sh
   ```

## Script Details

The script performs the following steps:

1. Changes to the project directory
2. Loads the CellRanger module (version 7.2.0)
3. Creates a configuration file (`cellranger_config.csv`) with:
   - Paths to mouse reference genome and probe set
   - FASTQ file location
   - Sample information and probe barcode IDs for all 15 samples
4. Runs CellRanger multi using the generated configuration file

## Sample Information

The script includes information for 15 samples:

1. Vehicle treated female 1 (Sample ID: 968, Barcode: BC001)
2. Vehicle treated female 2 (Sample ID: 966, Barcode: BC002)
3. Vehicle treated female 3 (Sample ID: 979, Barcode: BC003)
4. Inhibitor treated female 1 (Sample ID: 963, Barcode: BC004)
5. Inhibitor treated female 2 (Sample ID: 975, Barcode: BC005)
6. Inhibitor treated female 3 (Sample ID: 964, Barcode: BC006)
7. Vehicle treated male 1 (Sample ID: 960, Barcode: BC007)
8. Vehicle treated male 2 (Sample ID: 978, Barcode: BC008)
9. Vehicle treated male 3 (Sample ID: 969, Barcode: BC009)
10. Vehicle treated male 4 (Sample ID: 983, Barcode: BC010)
11. Inhibitor treated male 1 (Sample ID: 0, Barcode: BC011)
12. Inhibitor treated male 2 (Sample ID: 959, Barcode: BC012)
13. Inhibitor treated male 3 (Sample ID: 965, Barcode: BC013)
14. Inhibitor treated male 4 (Sample ID: 967, Barcode: BC014)
15. Inhibitor treated male 5 (Sample ID: 977, Barcode: BC015)

## Output

- The main output directory will be named `Chrm_499_OHagan_IUB_FLEX15_run6` (or similar, depending on the run number)
- SLURM output and error logs: `cellranger_multi.out` and `cellranger_multi.error`
- The configuration file: `cellranger_config.csv`

## Customization

- Adjust the SLURM directives (e.g., `--time`, `--mem`, `--cpus-per-task`) based on your job requirements and system policies
- Modify the sample information in the script if your experiment setup differs

## Troubleshooting

- If the job fails, check the error log (`cellranger_multi.error`)
- Ensure all paths in the script are correct for your system
- Verify that the probe barcode IDs (BC001, BC002, etc.) match your experimental setup

## Support

For issues related to this script or the project, please contact:
- Ateeq Khaliq (Project Lead)
- Your local bioinformatics support team

For CellRanger-specific issues, refer to the [10x Genomics support page](https://support.10xgenomics.com/).

## Project Affiliation

This work is part of Heather's Syngenic Mouse Models project, focusing on spatial transcriptomics analysis of mouse tissue samples to investigate gene expression patterns in various experimental conditions.


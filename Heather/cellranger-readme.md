# CellRanger Analysis Script for Visium FLEX Data

## Overview

This SLURM script automates the process of running CellRanger multi on Visium FLEX spatial transcriptomics data. It creates a configuration file with the correct probe barcode IDs and runs CellRanger using this configuration.

## Requirements

- Access to a SLURM-based HPC system
- CellRanger (version 7.2.0) installed and accessible as a module
- Sufficient storage space for CellRanger output (varies based on your data size, but typically >500GB)
- Input data:
  - FASTQ files from your Visium FLEX experiment
  - Reference genome
  - Probe set file

## Directory Structure

Ensure your project directory is structured as follows:

```
/N/project/cytassist/heather/Chrm_499_OHagan_IUB_FLEX15_Mouse_May2024/
├── cellranger-slurm-script.sh
├── *.fastq.gz (Your FASTQ files)
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
2. Loads the CellRanger module
3. Creates a configuration file (`cellranger_config.csv`) with:
   - Paths to reference genome and probe set
   - FASTQ file location
   - Sample information and probe barcode IDs
4. Runs CellRanger multi using the generated configuration file

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

For issues related to this script, please contact your local bioinformatics support team. For CellRanger-specific issues, refer to the [10x Genomics support page](https://support.10xgenomics.com/).


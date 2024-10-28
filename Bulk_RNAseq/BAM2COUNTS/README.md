# BAM Counting and Normalization Pipeline

This repository contains a pipeline to count reads from RNA-seq BAM files, convert Ensembl IDs to gene names, and normalize counts using DESeq2. The pipeline is designed for HPC environments using SLURM and can be customized for other setups.

## Table of Contents

- [Overview](#overview)
- [Prerequisites](#prerequisites)
- [Pipeline Details](#pipeline-details)
- [Files Generated](#files-generated)
- [Usage](#usage)
- [Output Files](#output-files)
- [Contributing](#contributing)
- [License](#license)

## Overview

This pipeline processes RNA-seq BAM files by performing the following steps:

1. **Indexing BAM files**: Reindexes BAM files if required.
2. **Counting Reads**: Uses R’s `GenomicAlignments` package to count reads that map to exons.
3. **Ensembl ID Conversion**: Converts Ensembl gene IDs to gene names using `org.Hs.eg.db`.
4. **Normalization**: Normalizes the raw counts using DESeq2, yielding normalized counts.
5. **Logging**: Generates logs that detail the workflow’s steps and any encountered errors.

## Prerequisites

Ensure that your HPC system has the following tools and libraries installed:

- **SLURM**: To schedule jobs on the HPC.
- **Miniconda**: To manage dependencies.
- **SAMtools**: For BAM file indexing.

The required R packages are installed by the script if they aren’t already available. These include:

- `GenomicAlignments`
- `GenomicFeatures`
- `Rsamtools`
- `rtracklayer`
- `DESeq2`
- `BiocParallel`
- `GenomicRanges`
- `org.Hs.eg.db`
- `AnnotationDbi`

## Pipeline Details

### Script Workflow

The pipeline script follows these steps:

- **SLURM Job Setup**: Configures job parameters like memory, CPUs, partition, and email notifications.
- **Environment Setup**: Loads necessary modules and activates the Conda environment.
- **Output Directory Creation**: Creates a timestamped output directory to store results.
- **BAM Reindexing**: Checks and updates BAM index files if they are outdated or missing.
- **R Script Generation**: Creates and executes an R script to:
  - Identify BAM files.
  - Count reads mapped to exons.
  - Convert Ensembl IDs to gene names.
  - Normalize counts using DESeq2.
- **Logging**: Saves process information in a log file.

### SLURM Job Options

The SLURM job uses the following options:

- **Job name**: `bam_count`
- **Nodes**: 1
- **CPUs per task**: 16
- **Memory**: 500GB
- **Partition**: GPU
- **Time limit**: 40 hours
- **Notifications**: Sent for job start, end, and fail events

## Files Generated

- **Raw Counts**: A CSV file with raw counts and Ensembl IDs.
- **Normalized Counts**: A CSV file with normalized counts and gene names.
- **Log File**: Captures details of the run, including errors.

## Usage

### Step 1: Clone the Repository

```bash

cd yourrepository

### Step 2: Submit the SLURM Job
Make the script executable and submit it to SLURM:

```bash
chmod +x bam_count_pipeline.sh
sbatch bam_count_pipeline.sh

```bash
### Step 3: Check the Output
The results will be saved in a timestamped directory under /N/project/akhaliq/Ateeq_dwd.

Output Files
raw_counts.csv: Contains raw counts with Ensembl IDs.
normalized_counts.csv: Contains normalized counts with gene names.
bam_counting.log: Log file with details of the run.
Contributing
Contributions are welcome! Please open an issue or submit a pull request.

License
This project is licensed under the MIT License - see the LICENSE file for details.

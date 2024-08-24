# Complete Guide: Single-cell and Spatial Transcriptomics in Oncology

## Table of Contents
1. [Introduction](#1-introduction)
2. [Single-cell RNA Sequencing (scRNA-seq)](#2-single-cell-rna-sequencing-scrna-seq)
3. [Spatial Transcriptomics](#3-spatial-transcriptomics)
4. [Data Analysis Techniques](#4-data-analysis-techniques)
5. [Applications in Oncology](#5-applications-in-oncology)
6. [Integration of scRNA-seq and Spatial Transcriptomics](#6-integration-of-scrna-seq-and-spatial-transcriptomics)
7. [Best Practices and Considerations](#7-best-practices-and-considerations)
8. [Future Directions](#8-future-directions)
9. [Resources and Tools](#9-resources-and-tools)
10. [Interview Preparation Tips](#10-interview-preparation-tips)
11. [Code Examples and Comparative Analysis](#11-code-examples-and-comparative-analysis)

## 1. Introduction

Single-cell and spatial transcriptomics have revolutionized our understanding of cellular heterogeneity and tissue organization in cancer. These technologies allow researchers to study gene expression patterns at unprecedented resolution, revealing insights into tumor complexity, microenvironment interactions, and treatment responses.

This guide is designed to prepare you for an interview for a senior scientist role in bioinformatics, focusing on single-cell and spatial transcriptomics in oncology. It covers key concepts, methodologies, analysis techniques, and applications specific to cancer research.

## 2. Single-cell RNA Sequencing (scRNA-seq)

### 2.1 Overview

scRNA-seq allows the measurement of gene expression in individual cells, providing a high-resolution view of cellular heterogeneity within tumors.

### 2.2 Key Technologies

1. **10x Genomics Chromium**
   - Uses gel beads in emulsion (GEM) to capture individual cells
   - High throughput (thousands of cells)
   - Relatively shallow sequencing depth per cell

2. **Smart-seq2**
   - Full-length transcript coverage
   - Lower throughput but higher sensitivity
   - Useful for detecting splice variants and fusion transcripts

3. **Drop-seq**
   - Uses microfluidic droplets to capture cells
   - Cost-effective for large numbers of cells
   - Lower sensitivity compared to Smart-seq2

### 2.3 Workflow

1. **Sample Preparation**: Dissociation of tissue into single-cell suspension
2. **Cell Isolation**: Encapsulation of individual cells (e.g., in droplets or wells)
3. **RNA Capture**: Reverse transcription and barcoding of mRNA
4. **Library Preparation**: Amplification and addition of sequencing adapters
5. **Sequencing**: Usually performed on Illumina platforms
6. **Data Analysis**: Computational processing and interpretation of results

### 2.4 Key Concepts

- **Unique Molecular Identifiers (UMIs)**: Random sequences used to tag individual mRNA molecules, allowing for the correction of amplification bias
- **Cell Barcodes**: Sequences used to identify reads originating from the same cell
- **Dropout Events**: When a gene is not detected in a cell due to technical limitations, even though it's actually expressed

## 3. Spatial Transcriptomics

### 3.1 Overview

Spatial transcriptomics techniques allow the measurement of gene expression while preserving spatial information within tissue sections.

### 3.2 Key Technologies

1. **10x Genomics Visium**
   - Uses a slide with pre-printed oligonucleotides to capture mRNA from tissue sections
   - Provides spatial resolution at ~55µm (spots containing 1-10 cells)

2. **MERFISH (Multiplexed Error-Robust Fluorescence In Situ Hybridization)**
   - Uses sequential rounds of RNA FISH to detect hundreds to thousands of genes
   - Provides single-cell resolution and subcellular localization of transcripts

3. **Slide-seq**
   - Uses DNA-barcoded beads on a slide to capture mRNA from tissue sections
   - Achieves near-single-cell resolution (~10µm)

### 3.3 Workflow

1. **Tissue Preparation**: Cryosectioning or FFPE tissue processing
2. **mRNA Capture**: Hybridization of mRNA to spatially barcoded probes
3. **Library Preparation**: Reverse transcription, amplification, and sequencing adapter addition
4. **Sequencing**: Usually performed on Illumina platforms
5. **Image Analysis**: Alignment of sequencing data with tissue images
6. **Data Analysis**: Spatial statistics and visualization of gene expression patterns

## 4. Data Analysis Techniques

### 4.1 Quality Control

- **Cell Filtering**: Remove low-quality cells based on:
  - Number of genes detected (e.g., 200-5000 genes per cell)
  - Number of UMIs (e.g., 500-20000 UMIs per cell)
  - Percentage of mitochondrial genes (e.g., <10-20%)

- **Gene Filtering**: Remove genes detected in very few cells

- **Doublet Detection**: Use algorithms like DoubletFinder or Scrublet to identify and remove cell doublets

### 4.2 Normalization

- **CPM (Counts Per Million)**: Scale counts to a common library size
- **SCTransform**: Regress out technical factors while preserving biological variability
- **Log-normalization**: Log-transform counts after adding a pseudocount

### 4.3 Feature Selection

- **Highly Variable Genes (HVGs)**: Identify genes with high cell-to-cell variability
- **Principal Component Analysis (PCA)**: Select top PCs for downstream analysis

### 4.4 Dimensionality Reduction

- **PCA**: Linear dimensionality reduction
- **t-SNE**: Non-linear dimensionality reduction, good for visualizing clusters
- **UMAP**: Non-linear dimensionality reduction, better at preserving global structure

### 4.5 Clustering

- **Graph-based Clustering**: Construct a K-nearest neighbor (KNN) graph and apply community detection algorithms
- **K-means Clustering**: Partition cells into K clusters based on their gene expression profiles

### 4.6 Differential Expression Analysis

- **MAST**: Zero-inflated regression model, accounts for bimodal expression in scRNA-seq
- **Wilcoxon Rank-Sum Test**: Non-parametric test for identifying differentially expressed genes between groups

### 4.7 Trajectory Inference

- **Monocle3**: Constructs a principal graph to model cellular differentiation trajectories
- **RNA Velocity**: Uses spliced and unspliced mRNA ratios to predict future cell states

## 5. Applications in Oncology

### 5.1 Tumor Heterogeneity Analysis

- Identify and characterize distinct cell populations within tumors
- Study clonal evolution and tumor progression

Example Research Question: "How does intratumoral heterogeneity change during treatment response in triple-negative breast cancer?"

### 5.2 Tumor Microenvironment

- Profile immune cell infiltration and activation states
- Characterize stromal cell populations and their interactions with tumor cells

Example Research Question: "What is the spatial distribution of different T cell subsets in relation to tumor cells in colorectal cancer, and how does this correlate with patient outcomes?"

### 5.3 Cancer Stem Cells

- Identify and characterize stem-like cell populations
- Study lineage dynamics and differentiation trajectories

Example Research Question: "Can we identify a gene signature for glioblastoma stem cells using scRNA-seq, and how does this signature change spatially within the tumor?"

### 5.4 Drug Response and Resistance

- Analyze single-cell drug perturbation experiments
- Identify resistant cell populations and mechanisms of resistance

Example Research Question: "Using scRNA-seq, can we identify pre-existing resistant cell populations in melanoma before BRAF inhibitor treatment, and what are their transcriptional characteristics?"

### 5.5 Biomarker Discovery

- Identify cell type-specific markers for improved diagnostics
- Discover spatial biomarkers indicative of cancer progression or treatment response

Example Research Question: "Can we use spatial transcriptomics to identify region-specific biomarkers in pancreatic cancer that correlate with metastatic potential?"

## 6. Integration of scRNA-seq and Spatial Transcriptomics

### 6.1 Data Integration Methods

- **Anchor-based Integration (Seurat)**: Uses shared features between datasets to align them
- **Harmony**: Iterative algorithm that projects cells into a shared embedding
- **LIGER**: Uses integrative non-negative matrix factorization to align datasets

### 6.2 Mapping scRNA-seq Data to Spatial Coordinates

- **Cell2location**: Probabilistic model to map scRNA-seq data to spatial transcriptomics data
- **SPOTlight**: Uses NMF to decompose spatial transcriptomics spots into cell type proportions
- **Tangram**: Uses optimal transport to map single-cell data to spatial data

## 7. Best Practices and Considerations

### 7.1 Experimental Design

- Include biological replicates to account for variability
- Consider time points for capturing dynamic processes
- Use appropriate controls (e.g., healthy tissue, treatment-naive samples)

### 7.2 Batch Effect Mitigation

- Use well-designed experimental batches
- Implement computational correction methods (e.g., ComBat, Harmony)

### 7.3 Data Management

- Use standardized data formats (e.g., AnnData, Seurat objects)
- Implement version control for analysis scripts (e.g., Git)
- Document analysis pipelines thoroughly

### 7.4 Reproducibility

- Use containerization (e.g., Docker) to ensure consistent environments
- Implement workflow management systems (e.g., Snakemake, Nextflow)
- Provide detailed methods sections in publications

## 8. Future Directions

### 8.1 Multi-omics Integration

- Single-cell multi-omics technologies (e.g., CITE-seq for simultaneous protein and RNA profiling)
- Integration of scRNA-seq with scATAC-seq for epigenomic insights

### 8.2 Spatial Multi-omics

- Combining spatial transcriptomics with spatial proteomics or spatial epigenomics
- Higher resolution spatial transcriptomics technologies

### 8.3 Machine Learning and AI

- Deep learning for cell type annotation and trajectory inference
- Graph neural networks for modeling cellular interaction networks
- AI-driven experimental design and analysis pipelines

## 9. Resources and Tools

### 9.1 R Packages

- Seurat: Comprehensive toolkit for scRNA-seq analysis
- Monocle3: Trajectory inference and differential expression
- Giotto: Spatial transcriptomics analysis

### 9.2 Python Packages

- Scanpy: Analysis of single-cell gene expression data
- Squidpy: Analysis and visualization of spatial molecular data
- scvi-tools: Deep probabilistic analysis of single-cell omics data

### 9.3 Databases and Atlases

- Human Cell Atlas: Comprehensive map of human cells
- The Cancer Genome Atlas (TCGA): Multi-omics data for various cancer types
- GEO (Gene Expression Omnibus): Repository for high-throughput sequencing data

## 10. Interview Preparation Tips

1. **Review Key Concepts**: Ensure you have a solid understanding of scRNA-seq and spatial transcriptomics technologies and analysis methods.

2. **Practice Data Analysis**: Work through publicly available datasets to gain hands-on experience with analysis pipelines.

3. **Stay Current**: Read recent publications in the field to familiarize yourself with cutting-edge applications and methodologies.

4. **Prepare Case Studies**: Be ready to discuss specific examples of how you've applied these technologies in oncology research.

5. **Understand Limitations**: Be prepared to discuss the limitations of these technologies and how to address them.

6. **Interdisciplinary Knowledge**: Brush up on relevant areas of cancer biology, immunology, and statistics.

7. **Problem-Solving Skills**: Practice explaining your approach to complex analytical challenges.

8. **Communication Skills**: Prepare to explain complex concepts to both technical and non-technical audiences.

[Previous content remains the same]

## 11. Code Examples and Comparative Analysis

### 11.1 Seurat Workflow for scRNA-seq Analysis

Here's a basic Seurat workflow for analyzing single-cell RNA-seq data:

```r
# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)

# Load the data (replace with your data file)
data <- Read10X(data.dir = "path/to/your/10x/data")
seurat_object <- CreateSeuratObject(counts = data, project = "MyProject", min.cells = 3, min.features = 200)

# Perform quality control
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter cells
seurat_object <- subset(seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalize data
seurat_object <- NormalizeData(seurat_object)

# Identify highly variable features
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)

# Scale data
all_genes <- rownames(seurat_object)
seurat_object <- ScaleData(seurat_object, features = all_genes)

# Perform linear dimensional reduction
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))

# Determine dimensionality
ElbowPlot(seurat_object)

# Cluster cells
seurat_object <- FindNeighbors(seurat_object, dims = 1:15)
seurat_object <- FindClusters(seurat_object, resolution = 0.5)

# Run non-linear dimensional reduction (UMAP)
seurat_object <- RunUMAP(seurat_object, dims = 1:15)

# Visualize results
DimPlot(seurat_object, reduction = "umap")

# Find markers for each cluster
markers <- FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top_markers <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

# Visualize top markers
VlnPlot(seurat_object, features = top_markers$gene[1:5], ncol = 5)
```

### 11.2 Spatial Transcriptomics Analysis with Seurat

Here's a basic workflow for analyzing spatial transcriptomics data using Seurat:

```r
library(Seurat)
library(ggplot2)
library(patchwork)

# Load the spatial data (replace with your data file)
spatial_data <- Load10X_Spatial(data.dir = "path/to/your/spatial/data")

# Normalize the data
spatial_data <- SCTransform(spatial_data, assay = "Spatial", verbose = FALSE)

# Visualize the spatial data
SpatialFeaturePlot(spatial_data, features = c("gene1", "gene2"))

# Identify spatially variable features
spatial_data <- FindSpatiallyVariableFeatures(spatial_data, assay = "SCT", features = VariableFeatures(spatial_data)[1:1000], selection.method = "markvariogram")

top_features <- head(SpatiallyVariableFeatures(spatial_data, selection.method = "markvariogram"), 6)
SpatialFeaturePlot(spatial_data, features = top_features, ncol = 3)

# Perform dimensionality reduction and clustering
spatial_data <- RunPCA(spatial_data, assay = "SCT", verbose = FALSE)
spatial_data <- FindNeighbors(spatial_data, reduction = "pca", dims = 1:30)
spatial_data <- FindClusters(spatial_data, verbose = FALSE)
spatial_data <- RunUMAP(spatial_data, reduction = "pca", dims = 1:30)

# Visualize clusters
DimPlot(spatial_data, reduction = "umap", label = TRUE)
SpatialDimPlot(spatial_data, label = TRUE, label.size = 3)
```

### 11.3 Comparative Analysis: Single-cell vs. Bulk vs. Spatial Transcriptomics

#### Data Structure and Resolution

1. **Bulk RNA-seq**:
   - Data structure: Gene expression matrix (genes x samples)
   - Resolution: Average expression across all cells in a sample
   - Computational complexity: Lowest

2. **Single-cell RNA-seq**:
   - Data structure: Gene expression matrix (genes x cells)
   - Resolution: Individual cell level
   - Computational complexity: High due to the large number of cells

3. **Spatial Transcriptomics**:
   - Data structure: Gene expression matrix (genes x spatial locations) + spatial coordinates
   - Resolution: Varies (from near single-cell to multi-cell resolution)
   - Computational complexity: High, with additional spatial components

#### Key Computational Differences

1. **Dimensionality**:
   - Bulk: Low-dimensional (number of samples)
   - Single-cell: High-dimensional (number of cells)
   - Spatial: High-dimensional with additional spatial dimensions

2. **Sparsity**:
   - Bulk: Dense data
   - Single-cell: Sparse data (many zero counts)
   - Spatial: Varies, often less sparse than scRNA-seq but sparser than bulk

3. **Normalization**:
   - Bulk: Global normalization methods (e.g., TPM, RPKM)
   - Single-cell: Cell-specific normalization (e.g., SCTransform)
   - Spatial: Similar to single-cell, with additional spatial normalization

4. **Batch Effect Correction**:
   - Bulk: Standard methods (e.g., ComBat)
   - Single-cell: Specialized methods (e.g., Harmony, LIGER)
   - Spatial: Similar to single-cell, with spatial batch effects consideration

5. **Dimensionality Reduction**:
   - Bulk: PCA, t-SNE
   - Single-cell: PCA followed by t-SNE or UMAP
   - Spatial: Similar to single-cell, with additional spatial dimensionality reduction techniques

6. **Clustering**:
   - Bulk: Hierarchical clustering, k-means
   - Single-cell: Graph-based clustering (e.g., Louvain, Leiden)
   - Spatial: Graph-based clustering with spatial constraints

7. **Differential Expression**:
   - Bulk: DESeq2, edgeR
   - Single-cell: MAST, Wilcoxon rank-sum test
   - Spatial: Spatial variance component analysis, SpatialDE

8. **Trajectory Inference**:
   - Bulk: Not applicable
   - Single-cell: Monocle, Slingshot
   - Spatial: Spatial trajectory methods (e.g., Giotto)

9. **Visualization**:
   - Bulk: Heatmaps, PCA plots
   - Single-cell: t-SNE/UMAP plots, feature plots
   - Spatial: Spatial feature plots, spatial dimension plots

10. **Integration with Other Data Types**:
    - Bulk: Straightforward integration with other bulk omics data
    - Single-cell: Requires specialized integration methods (e.g., MOFA+)
    - Spatial: Integration of spatial data with scRNA-seq (e.g., cell2location)

### 11.4 Computational Considerations

1. **Memory Usage**:
   - Bulk: Low to moderate
   - Single-cell: High (can require 100s of GB for large datasets)
   - Spatial: High, with additional memory for spatial data structures

2. **Processing Time**:
   - Bulk: Fast (minutes to hours)
   - Single-cell: Slow (hours to days for large datasets)
   - Spatial: Slow, with additional time for spatial calculations

3. **Scalability**:
   - Bulk: Easily scalable to many samples
   - Single-cell: Requires optimization for millions of cells (e.g., sketching techniques)
   - Spatial: Requires optimization for high-resolution or large tissue areas

4. **Data Storage**:
   - Bulk: Gigabytes
   - Single-cell: Terabytes for large-scale studies
   - Spatial: Terabytes, with additional storage for high-resolution images

5. **Quality Control**:
   - Bulk: Focus on sample-level QC
   - Single-cell: Cell-level and gene-level QC (e.g., doublet detection, ambient RNA removal)
   - Spatial: Additional QC for spatial artifacts and tissue quality

By understanding these computational differences, you can better appreciate the unique challenges and opportunities presented by each type of transcriptomic data in oncology research.

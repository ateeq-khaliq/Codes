# Single-cell and Spatial Transcriptomics in Oncology: A Comprehensive Guide

## Table of Contents
1. [Introduction](#introduction)
2. [Single-cell RNA Sequencing (scRNA-seq)](#single-cell-rna-sequencing-scrna-seq)
3. [Spatial Transcriptomics](#spatial-transcriptomics)
4. [Data Analysis Techniques](#data-analysis-techniques)
5. [Applications in Oncology](#applications-in-oncology)
6. [Integration of scRNA-seq and Spatial Transcriptomics](#integration-of-scrna-seq-and-spatial-transcriptomics)
7. [Best Practices and Considerations](#best-practices-and-considerations)
8. [Future Directions](#future-directions)
9. [Resources and Tools](#resources-and-tools)
10. [Interview Preparation Tips](#interview-preparation-tips)

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

Example (R with Seurat):
```r
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
```

### 4.3 Feature Selection

- **Highly Variable Genes (HVGs)**: Identify genes with high cell-to-cell variability
- **Principal Component Analysis (PCA)**: Select top PCs for downstream analysis

Example (R with Seurat):
```r
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
```

### 4.4 Dimensionality Reduction

- **PCA**: Linear dimensionality reduction
- **t-SNE**: Non-linear dimensionality reduction, good for visualizing clusters
- **UMAP**: Non-linear dimensionality reduction, better at preserving global structure

Example (R with Seurat):
```r
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
seurat_object <- RunUMAP(seurat_object, dims = 1:30)
```

### 4.5 Clustering

- **Graph-based Clustering**: Construct a K-nearest neighbor (KNN) graph and apply community detection algorithms
- **K-means Clustering**: Partition cells into K clusters based on their gene expression profiles

Example (R with Seurat):
```r
seurat_object <- FindNeighbors(seurat_object, dims = 1:30)
seurat_object <- FindClusters(seurat_object, resolution = 0.8)
```

### 4.6 Differential Expression Analysis

- **MAST**: Zero-inflated regression model, accounts for bimodal expression in scRNA-seq
- **Wilcoxon Rank-Sum Test**: Non-parametric test for identifying differentially expressed genes between groups

Example (R with Seurat):
```r
markers <- FindMarkers(seurat_object, ident.1 = "Cluster1", ident.2 = "Cluster2", test.use = "MAST")
```

### 4.7 Trajectory Inference

- **Monocle3**: Constructs a principal graph to model cellular differentiation trajectories
- **RNA Velocity**: Uses spliced and unspliced mRNA ratios to predict future cell states

Example (R with Monocle3):
```r
cds <- new_cell_data_set(expression_matrix, cell_metadata = cell_metadata, gene_metadata = gene_metadata)
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
```

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

Example (R with Seurat):
```r
integrated_data <- FindIntegrationAnchors(object.list = list(scRNA_data, spatial_data), anchor.features = 2000)
integrated_data <- IntegrateData(anchorset = integrated_data)
```

### 6.2 Mapping scRNA-seq Data to Spatial Coordinates

- **Cell2location**: Probabilistic model to map scRNA-seq data to spatial transcriptomics data
- **SPOTlight**: Uses NMF to decompose spatial transcriptomics spots into cell type proportions
- **Tangram**: Uses optimal transport to map single-cell data to spatial data

Example (Python with cell2location):
```python
import cell2location
model = cell2location.models.Cell2location(
    adata_vis, sc_data, N_cells_per_location=30, detection_alpha=20
)
model.train(max_epochs=1000, batch_size=2048, train_size=1)
```

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

Remember, the field of single-cell and spatial transcriptomics is rapidly evolving. Demonstrating your ability to learn and adapt to new technologies and methods is just as important as your current knowledge base.

Good luck with your interview preparation!

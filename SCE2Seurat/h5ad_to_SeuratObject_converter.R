library(SeuratObject)
library(Seurat)
library(zellkonverter)
sce <- readH5AD("/Users/akhaliq/Downloads/fetal_RAWCOUNTS_cellxgene.h5ad")

# Assuming `sce` is your SingleCellExperiment object
seurat_obj <- as.Seurat(sce, counts = "X", data = "X")

saveRDS(seurat_obj,"seurat_obj.rds")

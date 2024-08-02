***************************************************************************************************************************************
***************************************************************************************************************************************
***************************************************************************************************************************************

library(SingleR)
library(Seurat)


load("RData/stroma_singler.RData")

#Adding Metadata to the Suret Obj

surat<-singler_sub_epithelial$seurat

metadata <- read.table("Metadata.txt") # this is our Super metadata

surat@meta.data$orig.ident <- NULL # delete the available data here

epi_row<-row.names(surat@meta.data) # take the row names available with with epi suaret object

metadata$X <- rownames(metadata) # Bring row names as first column and name it as X

write.csv(epi_row,"epi_meta.csv",row.names=F)


epi_row<-read.table("epi_meta.csv",col.names="X")

new_meta <- merge(epi_row,metadata, by="X")


surat@meta.data$orig.ident <- new_meta$orig.ident # add sample Identification
surat@meta.data$Condition <- new_meta$Condition
surat@meta.data$Location <- new_meta$Location
surat@meta.data$MSI_Status <- new_meta$MSI_Status
surat@meta.data$Phase <- new_meta$Phase

#now remove the old surat obj from singular obj and replace it with the modified one
singler_sub_epithelial$seurat <- NULL

singler_sub_epithelial$seurat <- surat


out = SingleR.PlotTsne(singler.T.cells$singler[[1]]$SingleR.single,
singler.T.cells$meta.data$xy,do.label=FALSE,
do.letters = F,labels = singler.T.cells$singler[[1]]$SingleR.single$labels,
label.size = 4, dot.size = 3)




# now change in SinglerR also

singler_sub_epithelial$meta.data$orig.ident <- NULL
singler_sub_epithelial$meta.data$orig.ident <- surat$orig.ident
singler_sub_epithelial$meta.data$Condition <- surat$Condition
singler_sub_epithelial$meta.data$Location <- surat$Location
singler_sub_epithelial$meta.data$MSI_Status <- surat$MSI_Status
singler_sub_epithelial$meta.data$Phase <- surat$Phase

***************************************************************************************************************************************
##################### ENCODE
***************************************************************************************************************************************


# plot by orig.ident
out = SingleR.PlotTsne(singler$singler[[2]]$SingleR.single,
singler$meta.data$xy, do.label = TRUE, do.letters = F,
labels=singler$meta.data$orig.ident,label.size = 4,
dot.size = 3)

out = SingleR.PlotTsne(singler.T.cells$singler[[2]]$SingleR.single, do.label = TRUE, do.letters = F,
labels=singler.T.cells$meta.data$orig.ident,label.size = 4,
dot.size = 3)

png("TSNE_singler_by_sample.png", units="in", width=10, height=10, res=300)
out$p
dev.off()

# plot by cluster

out = SingleR.PlotTsne(singler_sub_epithelial$singler[[2]]$SingleR.single,
singler_sub_epithelial$meta.data$xy, do.label = TRUE, do.letters = F,
labels=singler_sub_epithelial$meta.data$cluster,label.size = 4,
dot.size = 3)

png("TSNE_singler_by_cluster.png", units="in", width=10, height=10, res=300)
out$p
dev.off()



png("heatmap_clas_before_fine_tune.png", units="in", width=10, height=10, res=300)
SingleR.DrawHeatmap(singler$singler[[2]]$SingleR.single.main, top.n = Inf,
clusters = singler$meta.data$orig.ident)
dev.off()

png("heatmap_clas_before_fine_tune_by_cluster.png", units="in", width=10, height=10, res=300)
SingleR.DrawHeatmap(singler_sub_epithelial$singler[[2]]$SingleR.single.main, top.n = Inf,
clusters = singler_sub_epithelial$meta.data$cluster)
dev.off()


png("heatmap_clas_before_fine_tune_by_cell_types.png", units="in", width=10, height=10, res=300)

SingleR.DrawHeatmap(singler_sub_epithelial$singler[[2]]$SingleR.single, top.n = 50,
clusters = singler_sub_epithelial$meta.data$orig.ident)
dev.off()

png("heatmap_clas_before_fine_tune_by_cell_types_clusters.png", units="in", width=10, height=10, res=300)

SingleR.DrawHeatmap(singler.cd8.cells$singler[[2]]$SingleR.single, top.n = 50,clusters = singler.cd8.cells$meta.data$clusters)
dev.off()



out = SingleR.PlotTsne(singler_sub_epithelial$singler[[2]]$SingleR.single,
singler_sub_epithelial$meta.data$xy,do.label=TRUE,
do.letters = F, labels = singler_sub_epithelial$singler[[2]]$SingleR.single$labels,
label.size = 4, dot.size = 3)

png("TSNE_singler_fine_tuned.png", units="in", width=10, height=10, res=300)
out$p
dev.off()


***************************************************************************************************************************************


#################. HPCA

***************************************************************************************************************************************



# plot by orig.ident
out = SingleR.PlotTsne(singler_sub_epithelial$singler[[1]]$SingleR.single,
singler_sub_epithelial$meta.data$xy, do.label = TRUE, do.letters = F,
labels=singler_sub_epithelial$meta.data$orig.ident,label.size = 4,
dot.size = 3)

png("TSNE_singler_by_sample_HPCA.png", units="in", width=10, height=10, res=300)
out$p
dev.off()

# plot by cluster

out = SingleR.PlotTsne(singler_sub_epithelial$singler[[1]]$SingleR.single,
singler_sub_epithelial$meta.data$xy, do.label = TRUE, do.letters = F,
labels=singler_sub_epithelial$meta.data$cluster,label.size = 4,
dot.size = 3)

png("TSNE_singler_by_cluster_HPCA.png", units="in", width=10, height=10, res=300)
out$p
dev.off()
dev.off()


png("heatmap_clas_before_fine_tune_HPCA.png", units="in", width=10, height=10, res=300)
SingleR.DrawHeatmap(singler_sub_epithelial$singler[[1]]$SingleR.single.main, top.n = Inf,
clusters = singler_sub_epithelial$meta.data$orig.ident)
dev.off()

png("heatmap_clas_before_fine_tune_by_cluster_HPCA.png", units="in", width=10, height=10, res=300)
SingleR.DrawHeatmap(singler_sub_epithelial$singler[[1]]$SingleR.single.main, top.n = Inf,
clusters = singler_sub_epithelial$meta.data$cluster)
dev.off()


png("heatmap_clas_before_fine_tune_by_cell_types_HPCA.png", units="in", width=10, height=10, res=300)

SingleR.DrawHeatmap(singler_sub_epithelial$singler[[1]]$SingleR.single, top.n = 50,
clusters = singler_sub_epithelial$meta.data$orig.ident)
dev.off()

png("heatmap_clas_before_fine_tune_by_cell_types_clusters_HPCA.png", units="in", width=10, height=10, res=300)

SingleR.DrawHeatmap(singler_sub_epithelial$singler[[1]]$SingleR.single, top.n = 50,
clusters = singler_sub_epithelial$meta.data$clusters)
dev.off()



out = SingleR.PlotTsne(singler_sub_epithelial$singler[[1]]$SingleR.single,
singler_sub_epithelial$meta.data$xy,do.label=TRUE,
do.letters = F, labels = singler_sub_epithelial$singler[[1]]$SingleR.single$labels,
label.size = 4, dot.size = 3)

png("TSNE_singler_fine_tuned_HPCA.png", units="in", width=50, height=30, res=300)
out$p
dev.off()

suppressPackageStartupMessages({
    library(Seurat)
    #library(venn)
    library(dplyr)
    library(cowplot)
    library(ggplot2)
    library(pheatmap)
    library(rafalib)
})


write.csv(table(singler_sub_epithelial$singler[[1]]$SingleR.single$labels,singler_sub_epithelial$meta.data$orig.ident), "labeling_samples.csv")


table(singler_sub_epithelial$meta.data$orig.ident,singler_sub_epithelial$singler[[1]]$SingleR.single.main$labels)



#################




suppressPackageStartupMessages({
    library(Seurat)
    library(venn)
    library(dplyr)
    library(cowplot)
    library(ggplot2)
    library(pheatmap)
    library(rafalib)
})


write.csv(table(singler_sub_epithelial$singler[[2]]$SingleR.single$labels,singler_sub_epithelial$meta.data$orig.ident), "labeling_samples.csv")


table(singler_sub_epithelial$meta.data$orig.ident,singler_sub_epithelial$singler[[2]]$SingleR.single.main$labels)


***************************************************************************************************************************************
To Subcluster based on the Conditions/Locations
***************************************************************************************************************************************
Ref: https://www.biostars.org/p/420831/


test.subset <- subset(x = alldata, subset = Condition == "Normal")

library(Seurat)
library(dbplyr)

data <- Read10X(data.dir = "C:/Users/haan1/Downloads/vdj_v1_hs_nsclc_5gex_filtered_gene_bc_matrices/")
alldata <- CreateSeuratObject(counts = data, project = "Ha An", min.cells = 3, min.features = 200)
alldata[["percent.mt"]] <- PercentageFeatureSet(alldata, pattern = "^MT-")
alldata <- subset(alldata, subset = nFeature_RNA>200 & nFeature_RNA<2500 & percent.mt <5)
alldata <- NormalizeData(alldata, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
alldata <- FindVariableFeatures(alldata, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
alldata <- ScaleData(alldata, verbose = FALSE)
alldata <- RunPCA(alldata, verbose = FALSE, features = VariableFeatures(object = alldata))
alldata <- RunTSNE(object = alldata, dims.use = 1:30)
alldata <- RunUMAP(alldata, dims = 1:30)

alldata <- FindNeighbors(alldata, dims = 1:30)
alldata <- FindClusters(alldata)


## Set Ident for a particular resolution

# Find variable genes based on the mean-dispersion relationship based on z-score for dispersion. 
pre_regressed_seurat <-  pre_regressed_seurat %>%
                          FindVariableFeatures(
                            mean.function = ExpMean,
                            dispersion.function = LogVMR,
                            do.plot = FALSE)


markers_genes <- FindAllMarkers(alldata, logfc.threshold = log(2), test.use = "wilcox", only.pos = TRUE, assay = "RNA")

write.csv(markers_genes, "Diff_exp_between_Normal_Tumor.csv")


 
write.table(markers_genes_0.8, file="Epi_subClust_Res_0.8.txt", sep="\t",append=F, row.names = F)
write.table(markers_genes_1, file="Epi_subClust_Res_1.txt", sep="\t",append=F, row.names = F)
write.table(markers_genes_1.5, file="Epi_subClust_Res_1.5.txt", sep="\t",append=F, row.names = F)
write.table(markers_genes_2, file="Epi_subClust_Res_2.txt", sep="\t",append=F, row.names = F)
write.table(markers_genes_2.5, file="Epi_subClust_Res_2.5.txt", sep="\t",append=F, row.names = F)
write.table(markers_genes_3, file="Epi_subClust_Res_3.txt", sep="\t",append=F, row.names = F)
write.table(markers_genes_3.5, file="Epi_subClust_Res_3.5.txt", sep="\t",append=F, row.names = F)
write.table(markers_genes_4, file="Epi_subClust_Res_4.txt", sep="\t",append=F, row.names = F)


top25 <- markers_genes %>% group_by(cluster) %>% top_n(-25, p_val_adj)



png("Diff_exp_heatmap_top25.png", units="in", width=15, height=20, res=300)

DoHeatmap(alldata, features = as.character(unique(top25$gene)), assay = "RNA")
dev.off()


# for all genes in Heatmap

png("Diff_exp_Epithelial_heatmap_ALL.png", units="in", width=40, height=50, res=300)

DoHeatmap(alldata, features = as.character(unique(top25$gene)), assay = "RNA")
dev.off()



png("Diff_exp_dotplot_top25.png", units="in", width=40, height=40, res=300)

DotPlot(alldata, features = as.character(unique(top5$gene)), assay = "RNA") + coord_flip()

dev.off()

png("Diff_exp_vinplot_top5.png", units="in", width=20, height=20, res=300)
VlnPlot(singler_sub_epithelial$seurat, features = as.character(unique(top5$gene)), ncol = 5, assay = "RNA")
dev.off()

saveRDS(alldata,"Epithelial_Normal_Vs_Tumor_052720.rds")
save.image("Epithelial_Normal_Vs_Tumor_052720.RData")
savehistory("Epithelial_Normal_Vs_Tumor_052720.RHistory")


png("Diff_exp_vinplot_top5.png", units="in", width=20, height=20, res=300)
VlnPlot(singler_sub_epithelial$seurat, features = as.character(unique(top5$gene)), ncol = 5, assay = "RNA")
dev.off()


c("tsne", "pheatmap", "MASS", "cluster", "mclust","flexmix", "lattice", "fpc", "RColorBrewer", "permute", "amap", "locfit","vegan", "Rtsne", "scran","randomForest","rgl")
*****************************************************************************************************************************************
############ 					Subsetting clusters
*****************************************************************************************************************************************

#load the res 3 

alldata <- readRDS("1_all_samples_Clustring_Diff_exp/RDS_Various_res/1_allsamples_Res_3_053120.rds")

library(Seurat)

# Subclustring Epithelial cells into a new cluster

a.epi <- SubsetData(alldata, subset.name = "seurat_clusters", accept.value = c(6,11,12,17,18,20,29,22))

table(a.epi$orig.ident)
table(a.epi$seurat_clusters)

a.epi = ae
'
DefaultAssay(a.epi) <- "RNA"

epithelial.list <- SplitObject(a.epi, split.by = "orig.ident")
epithelial.list
'
for (i in 1:length(x = epithelial.list)){
epithelial.list[[i]] <- NormalizeData(epithelial.list[[i]], verbose = FALSE)
epithelial.list[[i]] <- SCTransform(object = epithelial.list[[i]],verbose = FALSE,return.only.var.genes = FALSE)}

epithelial.features <- SelectIntegrationFeatures(object.list = epithelial.list, nfeatures = 3000)
epithelial.features

epithelial.list <- PrepSCTIntegration(object.list = epithelial.list, anchor.features = epithelial.features)


epithelial.anchors <- FindIntegrationAnchors(object.list = epithelial.list, normalization.method = "SCT", k.filter = 13, k.score = 12,dims = 1:25, anchor.features = epithelial.features)

epithelial.integrated <- IntegrateData(anchorset = epithelial.anchors, dims = 1:30)

epithelial.integrated<- ScaleData(epithelial.integrated, verbose = FALSE)

epithelial.integrated <- RunPCA(epithelial.integrated, npcs = 30, verbose = FALSE, reduction.name = "PCA_on_CCA")

epithelial.integrated<- RunUMAP(epithelial.integrated, reduction = "PCA_on_CCA", dims = 1:30,reduction.name = "UMAP_on_CCA")

epithelial.integrated <- RunTSNE(epithelial.integrated, reduction = "PCA_on_CCA", dims = 1:30,reduction.name = "TSNE_on_CCA")

png("test.png", units="in", width=15, height=15, res=300)
DimPlot(epithelial.integrated, reduction = "UMAP_on_CCA", group.by = "orig.ident")
dev.off()

png("test_pca.png", units="in", width=15, height=15, res=300)

DimPlot(epithelial.integrated, reduction = "TSNE_on_CCA", group.by = "orig.ident")
dev.off()

epithelial.integrated
alldata.conta.minus <- FindNeighbors(alldata.conta.minus, reduction = "PCA_on_CCA", dims = 1:30, k.param = 60, force.recalc = T, prune.SNN = 1/15)


names(epithelial.integrated@graphs)
epithelial.integrated_test_0.8 <- FindClusters(epithelial.integrated, graph.name = "integrated_snn", resolution = 0.8 , algorithm = 1)
epi_1.5 <- FindClusters(epithelial.integrated, graph.name = "integrated_snn", resolution = 1.5 , algorithm = 1)
epi_0.8 <- FindClusters(epithelial.integrated, graph.name = "integrated_snn", resolution = 0.8 , algorithm = 1)
epi_2 <- FindClusters(epithelial.integrated, graph.name = "integrated_snn", resolution = 2 , algorithm = 1)
epi_2.5 <- FindClusters(epithelial.integrated, graph.name = "integrated_snn", resolution = 2.5 , algorithm = 1)
epi_3 <- FindClusters(epithelial.integrated, graph.name = "integrated_snn", resolution = 3 , algorithm = 1)
epi_3.5 <- FindClusters(epithelial.integrated, graph.name = "integrated_snn", resolution = 3.5 , algorithm = 1)
epi_4 <- FindClusters(epithelial.integrated, graph.name = "integrated_snn", resolution = 4 , algorithm = 1)
epi_0.5 <- FindClusters(epithelial.integrated, graph.name = "integrated_snn", resolution = 0.5 , algorithm = 1)



suppressPackageStartupMessages({
  library(Seurat)
  library(cowplot)
  library(ggplot2)
  library(pheatmap)
  library(rafalib)
})

plot_grid(ncol = 3,
  DimPlot(epi_0.5, reduction = "UMAP_on_CCA", label = T)+ggtitle("louvain_0.5"),
  DimPlot(epi_0.8, reduction = "UMAP_on_CCA", label = T)+ggtitle("louvain_0.8"),
  DimPlot(epi_1.5, reduction = "UMAP_on_CCA", label = T)+ggtitle("louvain_1.5"),
  DimPlot(epi_2, reduction = "UMAP_on_CCA", label = T)+ggtitle("louvain_2"),
  DimPlot(epi_2.5, reduction = "UMAP_on_CCA", label = T)+ggtitle("louvain_2"),
  DimPlot(epi_3, reduction = "UMAP_on_CCA", label = T)+ggtitle("louvain_3"),
  DimPlot(epi_3.5, reduction = "UMAP_on_CCA", label = T)+ggtitle("louvain_3.5"),
  DimPlot(epi_4, reduction = "UMAP_on_CCA", label = T)+ggtitle("louvain_4")
)
dev.off()

markers_genes_0.5 <- FindAllMarkers(epi_0.5, logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, min.diff.pct = 0.2, only.pos = TRUE, max.cells.per.ident = 50, assay = "RNA")
markers_genes_0.8 <- FindAllMarkers(epi_0.8, logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, min.diff.pct = 0.2, only.pos = TRUE, max.cells.per.ident = 50, assay = "RNA")
markers_genes_1.5 <- FindAllMarkers(epi_1.5, logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, min.diff.pct = 0.2, only.pos = TRUE, max.cells.per.ident = 50, assay = "RNA")
markers_genes_2 <- FindAllMarkers(epi_2, logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, min.diff.pct = 0.2, only.pos = TRUE, max.cells.per.ident = 50, assay = "RNA")
markers_genes_2.5 <- FindAllMarkers(epi_2.5, logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, min.diff.pct = 0.2, only.pos = TRUE, max.cells.per.ident = 50, assay = "RNA")
markers_genes_3 <- FindAllMarkers(epi_3, logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, min.diff.pct = 0.2, only.pos = TRUE, max.cells.per.ident = 50, assay = "RNA")
markers_genes_3.5 <- FindAllMarkers(epi_3.5, logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, min.diff.pct = 0.2, only.pos = TRUE, max.cells.per.ident = 50, assay = "RNA")
markers_genes_4 <- FindAllMarkers(epi_4, logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, min.diff.pct = 0.2, only.pos = TRUE, max.cells.per.ident = 50, assay = "RNA")



saveRDS(epi_res_0.5,"Epithelial_Subcluster_res_0.5_060220.rds")
saveRDS(a.epi.n,"Epithelial_Subcluster_res_0.8_060220.rds")
saveRDS(epi_res_1.5,"Epithelial_Subcluster_res_1.5_060220.rds")
saveRDS(epi_res_2,"Epithelial_Subcluster_res_2_060220.rds")
saveRDS(epi_res_2.5,"Epithelial_Subcluster_res_2.5_060220.rds")
saveRDS(epi_res_3,"Epithelial_Subcluster_res_3_060220.rds")
saveRDS(epi_res_3.5,"Epithelial_Subcluster_res_3.5_060220.rds")
saveRDS(epi_res_4,"Epithelial_Subcluster_res_4_060220.rds")
Save.image("Epi_all_Res_0.5-4_060220.RData")



write.table(markers_genes_0.5, file="Epi_subClust_Res_0.5.txt", sep="\t",append=F, row.names = F)
write.table(markers_genes_0.8, file="Epi_subClust_Res_0.8.txt", sep="\t",append=F, row.names = F)
#write.table(markers_genes_1, file="Epi_subClust_Res_1.txt", sep="\t",append=F, row.names = F)
write.table(markers_genes_1.5, file="Epi_subClust_Res_1.5.txt", sep="\t",append=F, row.names = F)
write.table(markers_genes_2, file="Epi_subClust_Res_2.txt", sep="\t",append=F, row.names = F)
write.table(markers_genes_2.5, file="Epi_subClust_Res_2.5.txt", sep="\t",append=F, row.names = F)
write.table(markers_genes_3, file="Epi_subClust_Res_3.txt", sep="\t",append=F, row.names = F)
write.table(markers_genes_3.5, file="Epi_subClust_Res_3.5.txt", sep="\t",append=F, row.names = F)
write.table(markers_genes_4, file="Epi_subClust_Res_4.txt", sep="\t",append=F, row.names = F)


savehistory("epithelial_subcluster.RHistory")
save.image("epi_subcluster_res_0.5-4.RData")

saveRDS(conta.minus, "contamination_minus_Epi_cluster_res_0.8_060320.rds")
save.image("contamination_minus_Epi_cluster_res_0.8_060320.RData")



## Singler


epi_meta$nCount_RNA <- NULL
epi_meta$nFeature_RNA <- NULL
epi_meta$percent_mito <- NULL
epi_meta$percent_ribo <- NULL
epi_meta$S.Score <- NULL
epi_meta$G2M.Score <- NULL
epi_meta$CCA_snn_res.0.8 <- NULL
epi_meta$seurat_clusters <- NULL
epi_meta$CCA_snn_res.3 <- NULL
epi_meta$RNA_snn_res.0.8 <- NULL


counts<-as.matrix(tcells@assays$RNA@counts)

library(SingleR)
singler.nk.cells = CreateSinglerObject(counts_nk, project.name = "NK-Cells_subclust", annot = "nk_meta.txt", min.genes = 200,
  technology = "10X", species = "Human", citation = "",
  ref.list = list(), normalize.gene.length = F, variable.genes = "de",
  fine.tune = T, do.signatures = F, do.main.types = T, 
  reduce.file.size = FALSE, numCores = 2)

save.image("SingleR_NK.RData")


singler_tcells$seurat <- alldata
*****************************************************************************************************************************************

***************************************************************************************************************************************
##Diffexp
***************************************************************************************************************************************


markers_genes <- FindAllMarkers(singler_sub_epithelial$seurat, logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, min.diff.pct = 0.2, only.pos = TRUE, max.cells.per.ident = 50, assay = "RNA")

write.csv(markers_genes, "Diff_exp_between_clusts.csv")


top25 <- markers_genes %>% group_by(cluster) %>% top_n(-25, p_val_adj)


png("Diff_exp_bar_plots.png", units="in", width=30, height=30, res=300)

mypar(1, 26, mar = c(4, 6, 3, 1))
for (i in unique(top25$cluster)) {
    barplot(sort(setNames(top25$p_val, top25$gene)[top25$cluster == i], F), horiz = T, las = 1, main = paste0(i, " vs. rest"), border = "red", yaxs = "i")
    abline(v = c(0, 0.25), lty = c(1, 2))
}

dev.off()


top5 <- markers_genes %>% group_by(cluster) %>% top_n(-5, p_val_adj)


singler_sub_epithelial$seurat <- ScaleData(singler_sub_epithelial$seurat, features = as.character(unique(top5$gene)), assay = "RNA")


png("Diff_exp_heatmap_top5.png", units="in", width=15, height=20, res=300)

DoHeatmap(singler_sub_epithelial$seurat, features = as.character(unique(top5$gene)), assay = "RNA")
dev.off()

png("Diff_exp_dotplot_top5.png", units="in", width=40, height=40, res=300)

DotPlot(singler_sub_epithelial$seurat, features = as.character(unique(top5$gene)), assay = "RNA") + coord_flip()

dev.off()

png("Diff_exp_vinplot_top5.png", units="in", width=20, height=20, res=300)
VlnPlot(singler_sub_epithelial$seurat, features = as.character(unique(top5$gene)), ncol = 5, assay = "RNA")
dev.off()



saveRDS(singler_sub_epithelial$seurat,"Stroma_surat_obj_052320.rds")
save.image("singler_Stroma_052320.RData")
savehistory("singler_Stroma_052320.RHistory")


##### diff expression b/w Conditions

alldata <- SetIdent(alldata, value = "Condition") # here Normal Vs Tumor

markers_genes <- FindAllMarkers(alldata, logfc.threshold = 0.2, test.use = "wilcox", only.pos = FALSE, assay = "RNA")

# Compute differentiall expression
markers_genes <- FindAllMarkers(alldata, logfc.threshold = 0.2, test.use = "wilcox", min.pct = 0.1, min.diff.pct = 0.2, only.pos = TRUE, max.cells.per.ident = 50, assay = "RNA")



write.csv(markers_genes, "Diff_exp_between_Normal_Tumor.csv")

top25 <- markers_genes %>% group_by(cluster) %>% top_n(-25, p_val_adj)



png("Diff_exp_heatmap_top25.png", units="in", width=15, height=20, res=300)

DoHeatmap(alldata, features = as.character(unique(top25$gene)), assay = "RNA")
dev.off()


# for all genes in Heatmap

png("Diff_exp_Epithelial_heatmap_ALL.png", units="in", width=40, height=50, res=300)

DoHeatmap(alldata, features = as.character(unique(top25$gene)), assay = "RNA")
dev.off()



png("Diff_exp_dotplot_top25.png", units="in", width=40, height=40, res=300)

DotPlot(alldata, features = as.character(unique(top5$gene)), assay = "RNA") + coord_flip()

dev.off()

png("Diff_exp_vinplot_top5.png", units="in", width=20, height=20, res=300)
VlnPlot(singler_sub_epithelial$seurat, features = as.character(unique(top5$gene)), ncol = 5, assay = "RNA")
dev.off()

saveRDS(alldata,"Epithelial_Normal_Vs_Tumor_052720.rds")
save.image("Epithelial_Normal_Vs_Tumor_052720.RData")
savehistory("Epithelial_Normal_Vs_Tumor_052720.RHistory")


SingleR.DrawHeatmap(singler.T.cells$singler[[1]]$SingleR.single, top.n = 50, clusters = singler.T.cells$meta.data$orig.ident)

***************************************************************************************************************************************

###########

***************************************************************************************************************************************

# for getting average expression for individual patients/sampels
Idents(epi_tumor) <- "orig.ident"
table(epi_tumor@active.ident)


Average <- AverageExpression(alldata, assay = "RNA", return.seurat = T) #Get Mean Values
expr_matrix <- as.matrix(Average@assays$RNA@counts[,names(Average@active.ident)]) #Get Expression


singler.epi.cells = CreateSinglerObject(epi.counts, project.name = "Epi-cell_subclust", annot = NULL, min.genes = 200,
  technology = "10X", species = "Human", citation = "",
  ref.list = list(), normalize.gene.length = F, variable.genes = "de",
  fine.tune = T, do.signatures = F, do.main.types = T, 
  reduce.file.size = FALSE, numCores =6)

singler_tcells$seurat <- alldata

SingleR.DrawHeatmap(singler.B.cells$singler[[1]]$SingleR.single, top.n = 50,clusters = bcells.0.8_cr@meta.data$seurat_clusters,order.by.clusters=TRUE,cells_order=TRUE,fontsize_row = 20)

***************************************************************************************************************************************
***************************************************************************************************************************************
***************************************************************************************************************************************

# saving matx file from Seurat object

***************************************************************************************************************************************
***************************************************************************************************************************************
***************************************************************************************************************************************


library(Seurat)
library(DropletUtils)
tmpdir <- tempfile()
write10xCounts(x = alldata@assays$RNA@counts, path = "tmpdir")

#then in a new terminal move the file to desired location from tmpdir


# write the metadata / annotations into csv

write.csv(alldata@meta.data, 'uncleaned_metadata.csv')

#####

T_cell <- Seurat::Read10X("/Users/akhaliq/Desktop/samples_final/03302020_masood_colon/CAC4_B_EXOM_191017_NB551624_0033_AH5MG5BDXX/outs/filtered/")

Adding metadata to the max file and creating Seurat Object

tcell.metadata <- read.csv("/Users/akhaliq/Desktop/scrna_analysis/Analysis_nick_Bioturing/mtx_files/cleaned_counts/cleaned_metadata.csv", row.names=1)

T_Cell_sub <- CreateSeuratObject(T_cell, meta.data = tcell.metadata, project = "T-Cell_SubCluster")
#you get a warning message just ignore that

#####

#Plot PCA,UMAP,TSNE

png("UMAP_TSNE_PCA_Fibroblast.png", units="in", width=15, height=15, res=300)
plot_grid(ncol = 3,
  DimPlot(fibro_res_0.8, reduction = "pca", group.by = "orig.ident",label=T),
  DimPlot(fibro_res_0.8, reduction = "tsne", group.by = "orig.ident",label=T),
  DimPlot(fibro_res_0.8, reduction = "umap", group.by = "orig.ident",label=T)
)
dev.off()

png("UMAP_TSNE_PCA_Fibroblast_Right.png", units="in", width=10, height=10, res=300)
plot_grid(ncol = 3,
  DimPlot(alldata, reduction = "PCA_on_CCA", label=T),
  DimPlot(alldata, reduction = "TSNE_on_CCA", label=T),
  DimPlot(alldata, reduction = "UMAP_on_CCA", label=T)
)
dev.off()



png("integrated_RD_Tcells.png", units="in", width=15, height=15, res=300)
plot_grid(ncol = 3,
  
  DimPlot(alldata.int, reduction = "PCA_on_CCA", group.by = "orig.ident"),
  DimPlot(alldata.int, reduction = "TSNE_on_CCA", group.by = "orig.ident"),
  DimPlot(alldata.int, reduction = "UMAP_on_CCA", group.by = "orig.ident")
)
dev.off()

png("UMAP_default_res_seurat.png", units="in", width=10, height=10, res=300)
DimPlot(alldata, reduction = "UMAP_on_CCA",label=T)
dev.off()


png("Diff_exp_Tcells.png", units="in", width=20, height=20, res=300)

mypar(1, 15, mar = c(4, 6, 3, 1))
for (i in unique(markers_genes$cluster)) {
    barplot(sort(setNames(markers_genes$avg_logFC, markers_genes$gene)[markers_genes$cluster == i], F), horiz = T, las = 1, main = paste0(i, " vs. rest"), border = "white", yaxs = "i")
    abline(v = c(0, 0.25), lty = c(1, 2))
}

dev.off()



alldata <- ScaleData(alldata, features = as.character(unique(markers_genes$gene)), assay = "RNA")

png("Diff_exp_Tcells_HeatMap.png", units="in", width=40, height=60, res=300)
DoHeatmap(alldata, features = as.character(unique(markers_genes$gene)), group.by = "CCA_snn_res.0.8", assay = "RNA")
dev.off()



################################. EXTract raw counts and save iT INTO A CSV FILE


counts <- as.matrix(alldata@assays$RNA@counts) # try with counts if not present try with data
write.csv(alldata@assays$RNA@counts,"B-cell_raw_counts_matrix_051920.csv")



#################


#!/usr/bin/env Rscript

library(SingleR)
library (Seurat)

setwd("/data/Subclustring_crc_051520_Res_0.8/singler/epithelial/")
alldata<-readRDS("/data/Subclustring_crc_051520_Res_0.8/Epithelial_SubClust/Objects_Rdata_README/Epithelial_subCluster_Seurat_052020.rds")

counts<- as.matrix(alldata@assays$RNA@counts)


singler_sub_epithelial <- CreateSinglerSeuratObject(counts = counts, project.name = "CRC_SubClust_epithelial", species = "Human", npca = 30, technology = "X10", citation = "", ref.list = list(), normalize.gene.length = F, variable.genes = "de", fine.tune = T, do.signatures = T, do.main.types = T, reduce.file.size = F, numCores = 6, annot="/data/Subclustring_crc_051520_Res_0.8/Metadata.txt")


singler = CreateSinglerSeuratObject(counts = as.matrix(alldata@assays$RNA@counts), annot = cell_type$Cell_Type, project.name="PDAC",
                                    min.genes = 200, min.cells = 3, npca = 10, technology, species = "Human", citation,
                                    normalize.gene.length = F, ref.list = list(), fine.tune = F,
                                    do.signatures = T, do.main.type = T, reduce.file.size = F, numCores = 4)






singler_T-cell_NEW <- CreateSinglerSeuratObject(counts = counts, project.name = "CRC_SubClust_T-cells", species = "Human", npca = 30, technology = "X10", citation = "", ref.list = list(), normalize.gene.length = F, variable.genes = "de", fine.tune = T, do.signatures = T, do.main.types = T, reduce.file.size = F, numCores = 4, annot="tcell_metdata.txt")

saveRDS(singler_T-cell_NEW,"singler_T-cell_NEW.rds")


save.image("epithelial_singler.RData")
savehistory("epithelial_singler.RHistory")




library(SingleR)
library(Seurat)

setwd("/data/Subclustring_crc_051520_Res_0.8/singler/tcell/")

alldata<-readRDS("/data/Subclustring_crc_051520_Res_0.8/T-Cell_SubClust/our_pipeline/Objects_Rdata_README/tcells_seurat_obj.rds")

counts<- as.matrix(alldata@assays$RNA@counts)


singler_sub_tcell <- CreateSinglerSeuratObject(counts = counts, project.name = "CRC_SubClust_Tcells", species = "Human", npca = 15, technology = "X10", citation = "", ref.list = list(), normalize.gene.length = F, variable.genes = "de", fine.tune = T, do.signatures = T, do.main.types = T, reduce.file.size = FALSE, numCores = 6, annot="Metadata.txt")

save.image("Tcell_singler.RData")
savehistory("Tcell_singler.RHistory")



***************************************************************************************************************************************

#BioTuring Suggested Pipeline

***************************************************************************************************************************************


library(Seurat)
library(dbplyr)

data <- Read10X(data.dir = "C:/Users/haan1/Downloads/vdj_v1_hs_nsclc_5gex_filtered_gene_bc_matrices/")
alldata <- CreateSeuratObject(counts = data, project = "Ha An", min.cells = 3, min.features = 200)
alldata[["percent.mt"]] <- PercentageFeatureSet(alldata, pattern = "^MT-")
alldata <- subset(alldata, subset = nFeature_RNA>200 & nFeature_RNA<2500 & percent.mt <5)
alldata <- NormalizeData(alldata, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
alldata <- FindVariableFeatures(alldata, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
alldata <- ScaleData(alldata, verbose = FALSE)
alldata <- RunPCA(alldata, verbose = FALSE, features = VariableFeatures(object = alldata))
alldata <- RunTSNE(object = alldata, dims.use = 1:30)
alldata <- RunUMAP(alldata, dims = 1:30)

alldata <- FindNeighbors(alldata, dims = 1:30, force.recalc = T)
alldata <- FindClusters(alldata, force.recalc = T)

rownames(myeloid@meta.data)
names(myeloid$RNA_snn_res.0.8)

markers_genes_0.8 <- FindAllMarkers(epi_0.8, logfc.threshold = 0.2, test.use = "wilcox", only.pos = TRUE,assay = "RNA")
write.table(markers_genes_0.8, file="marker_genes_0.8.txt", sep="\t",append=F, row.names = F)
cbind(table(markers_genes_0.8$cluster))




png("RidgePlot_CD8_Tumor_reactive.png", units="in", width=20, height=15, res=300)
RidgePlot(cd8.0.8, features = c("HAVCR2","LAG3","ENTPD1","CXCL13"), ncol = 5, assay = "RNA")
dev.off()

png("FeaturePlot_CD8_Tumor_reactive.png", units="in", width=10, height=10, res=300)
FeaturePlot(cd8.0.8, reduction = "umap",features = c("HAVCR2","LAG3","ENTPD1","CXCL13") ,min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()



png("RidgePlot_CD8_Effector_cells.png", units="in", width=20, height=15, res=300)
RidgePlot(cd8.0.8, features = c("FGFBP2","GNLY"), ncol = 2, assay = "RNA")
dev.off()

png("FeaturePlot_CD8_Effector_cells.png", units="in", width=10, height=10, res=300)
FeaturePlot(cd8.0.8, reduction = "umap",features = c("FGFBP2","GNLY") ,min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()




png("RidgePlot_CD8_Naive_cells.png", units="in", width=10, height=10, res=300)
RidgePlot(cd8.0.8, features = c("CCR7","GZMK"), ncol = 2, assay = "RNA")
dev.off()

png("FeaturePlot_CD8_Naive_cells.png", units="in", width=10, height=10, res=300)
FeaturePlot(cd8.0.8, reduction = "umap",features = c("CCR7","GZMK") ,min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()


png("RidgePlot_CD8_Cent_memory_cells.png", units="in", width=10, height=10, res=300)
RidgePlot(cd8.0.8, features = c("CCR7","SELL"), ncol = 2, assay = "RNA")
dev.off()

png("FeaturePlot_CD8_Cent_memory_cells.png", units="in", width=10, height=10, res=300)
FeaturePlot(cd8.0.8, reduction = "umap",features = c("CCR7","SELL") ,min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()


png("RidgePlot_CD8_Mucosal_assosiated_cells.png", units="in", width=10, height=10, res=300)
RidgePlot(cd8.0.8, features = c("KLRB1","TRAV1","TRAV2","TRAJ33"), ncol = 4, assay = "RNA")
dev.off()

png("FeaturePlot_CD8_Mucosal_assosiated_cells.png", units="in", width=10, height=10, res=300)
FeaturePlot(cd8.0.8, reduction = "umap",features = c("KLRB1","TRAV1","TRAV2","TRAJ33") ,min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()



png("RidgePlot_CD8_Prolifarative_cells.png", units="in", width=10, height=10, res=300)
RidgePlot(cd8.0.8, features = "MKI67", ncol = 1, assay = "RNA")
dev.off()

png("FeaturePlot_CD8_Proliferative_cells.png", units="in", width=10, height=10, res=300)
FeaturePlot(cd8.0.8, reduction = "umap",features = "MKI67" ,min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()



png("RidgePlot_CD8_Chemokine_cells.png", units="in", width=10, height=10, res=300)
RidgePlot(cd8.0.8, features = c("XCL1","XCL2"), ncol = 2, assay = "RNA")
dev.off()

png("FeaturePlot_CD8_Chemokine_cells.png", units="in", width=10, height=10, res=300)
FeaturePlot(cd8.0.8, reduction = "umap",features = c("XCL1","XCL2") ,min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()

# All Plots

png("RidgePlot_CD4_All_cells.png", units="in", width=20, height=20, res=300)
RidgePlot(cd4_0.2, features = c("CXCL13","IFNG","TOX","FOXP3","IL2RA","TIGIT","TNFRSF4","TNFRSF9","TNFRSF18","CD27","CTLA4","CD69","GZMB","PRF1","NKG7","IL17A","STMN1","TUBB","PCNA","HMGB1","HMGB2"), ncol = 4, assay = "RNA")
dev.off()

png("FeaturePlot_CD4_ALL_cells.png", units="in", width=20, height=20, res=300)
FeaturePlot(cd4_0.2, reduction = "umap",features = c("CXCL13","IFNG","TOX","FOXP3","IL2RA","TIGIT","TNFRSF4","TNFRSF9","TNFRSF18","CD27","CTLA4","CD69","GZMB","PRF1","NKG7","IL17A","STMN1","TUBB","PCNA","HMGB1","HMGB2") ,min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()

##########
Renaming Clusters
##########

# Save old identity classes (the cluster labels) for reference.
stem.combined[["old.ident"]] <- Idents(object = stem.combined)

# Rename classes.
stem.combined <- RenameIdents(object = stem.combined, `0` = "your cell type", `1` = "your other cell type", `2` = "your last cell type")

##########


png("UMAP.png", units="in", width=15, height=10, res=300)

DimPlot(object = alldata_str_res_0.8, reduction = "pca")
DimPlot(object = alldata_str_res_0.8, reduction = "tsne")
DimPlot(object = alldata_str_res_0.8, reduction = "umap", label=T)

dev.off()


table(alldata@meta.data$seurat_clusters)

##Diff_exp

suppressPackageStartupMessages({
    library(Seurat)
    library(venn)
    library(dplyr)
    library(cowplot)
    library(ggplot2)
    library(pheatmap)
    library(rafalib)
})


# Compute differentiall expression
markers_genes <- FindAllMarkers(alldata, logfc.threshold = 0.2, test.use = "wilcox", min.pct = 0.1, min.diff.pct = 0.2, only.pos = TRUE, max.cells.per.ident = 50, assay = "RNA")


 markers_genes_1.5 <- FindAllMarkers(epi_res_1.5, logfc.threshold = 0.2, test.use = "wilcox",only.pos = TRUE, assay = "RNA")
 markers_genes_2 <- FindAllMarkers(epi_res_2, logfc.threshold = 0.2, test.use = "wilcox",only.pos = TRUE, assay = "RNA")
 markers_genes_2.5 <- FindAllMarkers(epi_res_2.5, logfc.threshold = 0.2, test.use = "wilcox",only.pos = TRUE, assay = "RNA")
 markers_genes_3 <- FindAllMarkers(epi_res_3, logfc.threshold = 0.2, test.use = "wilcox",only.pos = TRUE, assay = "RNA")
 markers_genes_3.5 <- FindAllMarkers(epi_res_3.5, logfc.threshold = 0.2, test.use = "wilcox",only.pos = TRUE, assay = "RNA")
 markers_genes_4 <- FindAllMarkers(epi_res_4, logfc.threshold = 0.2, test.use = "wilcox",only.pos = TRUE, assay = "RNA")

png("cell_label_Epithelial_marker.png", units="in", width=50, height=50, res=300)

FeaturePlot(alldata, reduction = "UMAP_on_CCA",features = c(),min.cutoff="q9",cols.use = c("lightgrey", "blue"), pt.size = 0.5)





png("Epi_Marker_lgr5.png", units="in", width=10, height=10, res=300)
FeaturePlot(alldata, reduction = "umap",features = "LGR5",min.cutoff="q9", cols.use = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()


***************************************************************************************************************************************

###### Clusttree Analysis


alldata_allres <- FindClusters(alldata, resolution = c(0,0.2,0.5,0.8,1.5,2,2.5,3,3.5,4) , algorithm = 1)
library(clustree)
png("korean_Subclust_Clusttree.png", units="in", width=10, height=15, res=300)
clustree(alldata_str.res.all, prefix = "RNA_snn_res.")
dev.off()

#########





# Epi Marker_cell

png("cell_label_Epithelial_marker_res_0.8.png", units="in", width=20, height=20, res=300)
FeaturePlot(alldata, reduction = "UMAP_on_CCA",features = c("EPCAM","KRT8","KRT18"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()



png("cell_label_Epithelial_marker_res_0.5.png", units="in", width=20, height=20, res=300)
FeaturePlot(alldata_0.5, reduction = "UMAP_on_CCA",features = c("EPCAM","KRT8","KRT18"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()


png("cell_label_Epithelial_marker_res_1.5.png", units="in", width=20, height=20, res=300)
FeaturePlot(alldata_1.5, reduction = "UMAP_on_CCA",features = c("EPCAM","KRT8","KRT18"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()



png("cell_label_Epithelial_marker_res_2.png", units="in", width=20, height=20, res=300)
FeaturePlot(alldata_2, reduction = "UMAP_on_CCA",features = c("EPCAM","KRT8","KRT18"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()



png("cell_label_Epithelial_marker_res_2.5.png", units="in", width=20, height=20, res=300)
FeaturePlot(alldata_2.5, reduction = "UMAP_on_CCA",features = c("EPCAM","KRT8","KRT18"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()


png("cell_label_Epithelial_marker_res_3.png", units="in", width=20, height=20, res=300)
FeaturePlot(alldata_3, reduction = "UMAP_on_CCA",features = c("EPCAM","KRT8","KRT18"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()


png("cell_label_Epithelial_marker_res_3.5.png", units="in", width=20, height=20, res=300)
FeaturePlot(alldata_3.5, reduction = "UMAP_on_CCA",features = c("EPCAM","KRT8","KRT18"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()



png("cell_label_Epithelial_marker_res_4.png", units="in", width=20, height=20, res=300)
FeaturePlot(alldata_4, reduction = "UMAP_on_CCA",features = c("EPCAM","KRT8","KRT18"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()

png("cell_label_Epithelial_marker_res_3_tsne.png", units="in", width=20, height=20, res=300)
FeaturePlot(alldata_3, reduction = "TSNE_on_CCA",features = c("EPCAM","KRT8","KRT18"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()




####stromal_marker


png("cell_label_Stromal_marker_res_0.8.png", units="in", width=20, height=20, res=300)
FeaturePlot(alldata, reduction = "UMAP_on_CCA",features = c("COL1A1", "COL1A2", "COL6A1","COL6A2", "VWF", "PLVAP", "CDH5", "S100B"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()



png("cell_label_Stromal_marker_res_0.5.png", units="in", width=20, height=20, res=300)
FeaturePlot(alldata_0.5, reduction = "UMAP_on_CCA",features = c("COL1A1", "COL1A2", "COL6A1","COL6A2", "VWF", "PLVAP", "CDH5", "S100B"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()


png("cell_label_Stromal_marker_res_1.5.png", units="in", width=20, height=20, res=300)
FeaturePlot(alldata_1.5, reduction = "UMAP_on_CCA",features = c("COL1A1", "COL1A2", "COL6A1","COL6A2", "VWF", "PLVAP", "CDH5", "S100B"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()



png("cell_label_Stromal_marker_res_2.png", units="in", width=20, height=20, res=300)
FeaturePlot(alldata_2, reduction = "UMAP_on_CCA",features = c("COL1A1", "COL1A2", "COL6A1","COL6A2", "VWF", "PLVAP", "CDH5", "S100B"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()



png("cell_label_Stromal_marker_res_2.5.png", units="in", width=20, height=20, res=300)
FeaturePlot(alldata_2.5, reduction = "UMAP_on_CCA",features = c("COL1A1", "COL1A2", "COL6A1","COL6A2", "VWF", "PLVAP", "CDH5", "S100B"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()


png("cell_label_Stromal_marker_res_3.png", units="in", width=20, height=20, res=300)
FeaturePlot(alldata_3, reduction = "UMAP_on_CCA",features = c("COL1A1", "COL1A2", "COL6A1","COL6A2", "VWF", "PLVAP", "CDH5", "S100B"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()


png("cell_label_Stromal_marker_res_3.5.png", units="in", width=20, height=20, res=300)
FeaturePlot(alldata_3.5, reduction = "UMAP_on_CCA",features = c("COL1A1", "COL1A2", "COL6A1","COL6A2", "VWF", "PLVAP", "CDH5", "S100B"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()



png("cell_label_Stromal_marker_res_4.png", units="in", width=20, height=20, res=300)
FeaturePlot(alldata_4, reduction = "UMAP_on_CCA",features = c("COL1A1", "COL1A2", "COL6A1","COL6A2", "VWF", "PLVAP", "CDH5", "S100B"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()

######### Immune Marker
  

png("cell_label_immune_marker_res_0.8.png", units="in", width=20, height=20, res=300)
FeaturePlot(alldata, reduction = "UMAP_on_CCA",features = c("CD52", "CD2", "CD3D", "CD3G", "CD3E", "CD79A", "CD79B", "CD14", "CD16", "CD68", "CD83", "CSF1R", "FCER1G"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()



png("cell_label_immune_marker_res_0.5.png", units="in", width=20, height=20, res=300)
FeaturePlot(alldata_0.5, reduction = "UMAP_on_CCA",features = c("CD52", "CD2", "CD3D", "CD3G", "CD3E", "CD79A", "CD79B", "CD14", "CD16", "CD68", "CD83", "CSF1R", "FCER1G"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()


png("cell_label_immune_marker_res_1.5.png", units="in", width=20, height=20, res=300)
FeaturePlot(alldata_1.5, reduction = "UMAP_on_CCA",features = c("CD52", "CD2", "CD3D", "CD3G", "CD3E", "CD79A", "CD79B", "CD14", "CD16", "CD68", "CD83", "CSF1R", "FCER1G"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()



png("cell_label_immune_marker_res_2.png", units="in", width=20, height=20, res=300)
FeaturePlot(alldata_2, reduction = "UMAP_on_CCA",features = c("CD52", "CD2", "CD3D", "CD3G", "CD3E", "CD79A", "CD79B", "CD14", "CD16", "CD68", "CD83", "CSF1R", "FCER1G"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()

png("cell_label_immune_marker_res_2.5.png", units="in", width=20, height=20, res=300)
FeaturePlot(alldata_2.5, reduction = "UMAP_on_CCA",features = c("CD52", "CD2", "CD3D", "CD3G", "CD3E", "CD79A", "CD79B", "CD14", "CD16", "CD68", "CD83", "CSF1R", "FCER1G"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()
#### B cells

png("NATURE_MED_B-CELL_marker_res_2.5.png", units="in", width=20, height=20, res=300)

FeaturePlot(alldata_2.5, reduction = "UMAP_on_CCA",features = c("CD79A", "IGKC", "IGLC3", "IGHG3"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)

dev.off()

#### T cells
 "CD3D", "TRBC1", "TRBC2", "TRAC" "CD3E", "CD3G"

png("NATURE_MED_T-CELL_marker_res_3_tsne.png", units="in", width=20, height=20, res=300)

FeaturePlot(alldata_3, reduction = "TSNE_on_CCA",features = c("CD3D", "TRBC1", "TRBC2", "TRAC", "CD3E", "CD3G"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)

dev.off()


"LYZ", "MARCO", "CD68","FCGR3A"
png("NATURE_MED_MYLOID_marker_res_2.5.png", units="in", width=20, height=20, res=300)

FeaturePlot(alldata_2.5, reduction = "UMAP_on_CCA",features = c("LYZ", "MARCO", "CD68","FCGR3A"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)

dev.off()
####

png("cell_label_immune_marker_res_3.png", units="in", width=20, height=20, res=300)
FeaturePlot(alldata_3, reduction = "UMAP_on_CCA",features = c("CD52", "CD2", "CD3D", "CD3G", "CD3E", "CD79A", "CD79B", "CD14", "CD16", "CD68", "CD83", "CSF1R", "FCER1G"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()

#### B cells

png("NATURE_MED_B-CELL_marker_res_3.png", units="in", width=20, height=20, res=300)

FeaturePlot(alldata_3, reduction = "UMAP_on_CCA",features = c("CD79A", "IGKC", "IGLC3", "IGHG3", "SDC1", "CD19", "MS4A1"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)

dev.off()

#### T cells
 "CD3D", "TRBC1", "TRBC2", "TRAC" "CD3E", "CD3G"

png("NATURE_MED_T-CELL_marker_res_3.png", units="in", width=20, height=20, res=300)

FeaturePlot(alldata_3, reduction = "UMAP_on_CCA",features = c("CD3D", "TRBC1", "TRBC2", "TRAC", "CD3E", "CD3G"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)

dev.off()

png("Myeloid_marker_res_3.png", units="in", width=20, height=20, res=300)

FeaturePlot(alldata_3, reduction = "UMAP_on_CCA",features = c("LYZ", "MARCO", "CD68","FCGR3A", "KIT"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)

dev.off()


######

png("cell_label_immune_marker_res_3.5.png", units="in", width=20, height=20, res=300)
FeaturePlot(alldata_3.5, reduction = "UMAP_on_CCA",features = c("CD52", "CD2", "CD3D", "CD3G", "CD3E", "CD79A", "CD79B", "CD14", "CD16", "CD68", "CD83", "CSF1R", "FCER1G"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()



png("cell_label_immune_marker_res_4.png", units="in", width=20, height=20, res=300)
FeaturePlot(alldata_4, reduction = "UMAP_on_CCA",features = c("CD52", "CD2", "CD3D", "CD3G", "CD3E", "CD79A", "CD79B", "CD14", "CD16", "CD68", "CD83", "CSF1R", "FCER1G"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()

###### Smelie_epithelial

"ASCL2","LGR5","LEFTY1","CDC25C","SLC26A2","FABP1","CA1","EPCAM","RBP2","SLC26A3","CCL20","TNFRSF11A","BEST4","TFF1","MUC2","HCK","TRPM5","SCGN","CHGA","CCL7","CCL13"



png("cell_SMILLIE_Epithelial_marker_res_0.8.png", units="in", width=40, height=40, res=300)
FeaturePlot(alldata, reduction = "UMAP_on_CCA",features = c("ASCL2","LGR5","LEFTY1","CDC25C","SLC26A2","FABP1","CA1","EPCAM","RBP2","SLC26A3","CCL20","TNFRSF11A","BEST4","TFF1","MUC2","HCK","TRPM5","SCGN","CHGA","CCL7","CCL13"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()



png("cell_SMILLIE_Epithelial_marker_res_0.5.png", units="in", width=40, height=40, res=300)
FeaturePlot(alldata_0.5, reduction = "UMAP_on_CCA",features = c("ASCL2","LGR5","LEFTY1","CDC25C","SLC26A2","FABP1","CA1","EPCAM","RBP2","SLC26A3","CCL20","TNFRSF11A","BEST4","TFF1","MUC2","HCK","TRPM5","SCGN","CHGA","CCL7","CCL13"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()


png("cell_SMILLIE_Epithelial_marker_res_1.5.png", units="in", width=40, height=40, res=300)
FeaturePlot(alldata_1.5, reduction = "UMAP_on_CCA",features = c("ASCL2","LGR5","LEFTY1","CDC25C","SLC26A2","FABP1","CA1","EPCAM","RBP2","SLC26A3","CCL20","TNFRSF11A","BEST4","TFF1","MUC2","HCK","TRPM5","SCGN","CHGA","CCL7","CCL13"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()



png("cell_SMILLIE_Epithelial_marker_res_2.png", units="in", width=40, height=40, res=300)
FeaturePlot(alldata_2, reduction = "UMAP_on_CCA",features = c("ASCL2","LGR5","LEFTY1","CDC25C","SLC26A2","FABP1","CA1","EPCAM","RBP2","SLC26A3","CCL20","TNFRSF11A","BEST4","TFF1","MUC2","HCK","TRPM5","SCGN","CHGA","CCL7","CCL13"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()



png("cell_SMILLIE_Epithelial_marker_res_2.5.png", units="in", width=40, height=40, res=300)
FeaturePlot(lt_nr_subset, reduction = "umap",features = c("ASCL2","LGR5","LEFTY1","CDC25C","SLC26A2","FABP1","CA1","EPCAM","RBP2","SLC26A3","CCL20","TNFRSF11A","BEST4","TFF1","MUC2","HCK","TRPM5","SCGN","CHGA","CCL7","CCL13"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()


png("cell_SMILLIE_Epithelial_marker_res_3.png", units="in", width=40, height=40, res=300)
FeaturePlot(alldata_3, reduction = "UMAP_on_CCA",features = c("ASCL2","LGR5","LEFTY1","CDC25C","SLC26A2","FABP1","CA1","EPCAM","RBP2","SLC26A3","CCL20","TNFRSF11A","BEST4","TFF1","MUC2","HCK","TRPM5","SCGN","CHGA","CCL7","CCL13"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()


png("cell_SMILLIE_Epithelial_marker_res_3.5.png", units="in", width=40, height=40, res=300)
FeaturePlot(alldata_3.5, reduction = "UMAP_on_CCA",features = c("ASCL2","LGR5","LEFTY1","CDC25C","SLC26A2","FABP1","CA1","EPCAM","RBP2","SLC26A3","CCL20","TNFRSF11A","BEST4","TFF1","MUC2","HCK","TRPM5","SCGN","CHGA","CCL7","CCL13"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()



png("cell_SMILLIE_Epithelial_marker_res_4.png", units="in", width=40, height=40, res=300)
FeaturePlot(alldata_4, reduction = "UMAP_on_CCA",features = c("ASCL2","LGR5","LEFTY1","CDC25C","SLC26A2","FABP1","CA1","EPCAM","RBP2","SLC26A3","CCL20","TNFRSF11A","BEST4","TFF1","MUC2","HCK","TRPM5","SCGN","CHGA","CCL7","CCL13"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()



###### SMILLIE Stromal

"RSPO3",	"CCL11", "WNT5B", "BMP4", "CHI3L1", "TAGLN", "ACTA2", "WNT2B", "SOX17", "RBP7", "VWA1", "PLVAP", "CDH5", "CD36", "MADCAM1", "DARC", "RGS5", "SOX10", "S100B", "CD68", "XCR1", "CLEC9A"


png("cell_SMILLIE_STROMAL_marker_res_0.8.png", units="in", width=40, height=40, res=300)
FeaturePlot(alldata, reduction = "UMAP_on_CCA",features = c("RSPO3",	"CCL11", "WNT5B", "BMP4", "CHI3L1", "TAGLN", "ACTA2", "WNT2B", "SOX17", "RBP7", "VWA1", "PLVAP", "CDH5", "CD36", "MADCAM1", "DARC", "RGS5", "SOX10", "S100B", "CD68", "XCR1", "CLEC9A"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()



png("cell_SMILLIE_STROMAL_marker_res_0.5.png", units="in", width=40, height=40, res=300)
FeaturePlot(alldata_0.5, reduction = "UMAP_on_CCA",features = c("RSPO3",	"CCL11", "WNT5B", "BMP4", "CHI3L1", "TAGLN", "ACTA2", "WNT2B", "SOX17", "RBP7", "VWA1", "PLVAP", "CDH5", "CD36", "MADCAM1", "DARC", "RGS5", "SOX10", "S100B", "CD68", "XCR1", "CLEC9A"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()


png("cell_SMILLIE_STROMAL_marker_res_1.5.png", units="in", width=40, height=40, res=300)
FeaturePlot(alldata_1.5, reduction = "UMAP_on_CCA",features = c("RSPO3",	"CCL11", "WNT5B", "BMP4", "CHI3L1", "TAGLN", "ACTA2", "WNT2B", "SOX17", "RBP7", "VWA1", "PLVAP", "CDH5", "CD36", "MADCAM1", "DARC", "RGS5", "SOX10", "S100B", "CD68", "XCR1", "CLEC9A"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()



png("cell_SMILLIE_STROMAL_marker_res_2.png", units="in", width=40, height=40, res=300)
FeaturePlot(alldata_2, reduction = "UMAP_on_CCA",features = c("RSPO3",	"CCL11", "WNT5B", "BMP4", "CHI3L1", "TAGLN", "ACTA2", "WNT2B", "SOX17", "RBP7", "VWA1", "PLVAP", "CDH5", "CD36", "MADCAM1", "DARC", "RGS5", "SOX10", "S100B", "CD68", "XCR1", "CLEC9A"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()



png("cell_SMILLIE_STROMAL_marker_res_2.5.png", units="in", width=40, height=40, res=300)
FeaturePlot(alldata_2.5, reduction = "UMAP_on_CCA",features = c("RSPO3",	"CCL11", "WNT5B", "BMP4", "CHI3L1", "TAGLN", "ACTA2", "WNT2B", "SOX17", "RBP7", "VWA1", "PLVAP", "CDH5", "CD36", "MADCAM1", "DARC", "RGS5", "SOX10", "S100B", "CD68", "XCR1", "CLEC9A"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()


png("cell_SMILLIE_STROMAL_marker_res_3.png", units="in", width=40, height=40, res=300)
FeaturePlot(alldata_3, reduction = "UMAP_on_CCA",features = c("RSPO3",	"CCL11", "WNT5B", "BMP4", "CHI3L1", "TAGLN", "ACTA2", "WNT2B", "SOX17", "RBP7", "VWA1", "PLVAP", "CDH5", "CD36", "MADCAM1", "DARC", "RGS5", "SOX10", "S100B", "CD68", "XCR1", "CLEC9A"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()


png("cell_SMILLIE_STROMAL_marker_res_3.5.png", units="in", width=40, height=40, res=300)
FeaturePlot(alldata_3.5, reduction = "UMAP_on_CCA",features = c("RSPO3",	"CCL11", "WNT5B", "BMP4", "CHI3L1", "TAGLN", "ACTA2", "WNT2B", "SOX17", "RBP7", "VWA1", "PLVAP", "CDH5", "CD36", "MADCAM1", "DARC", "RGS5", "SOX10", "S100B", "CD68", "XCR1", "CLEC9A"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()



png("cell_SMILLIE_STROMAL_marker_res_4.png", units="in", width=40, height=40, res=300)
FeaturePlot(alldata_4, reduction = "UMAP_on_CCA",features = c("RSPO3",	"CCL11", "WNT5B", "BMP4", "CHI3L1", "TAGLN", "ACTA2", "WNT2B", "SOX17", "RBP7", "VWA1", "PLVAP", "CDH5", "CD36", "MADCAM1", "DARC", "RGS5", "SOX10", "S100B", "CD68", "XCR1", "CLEC9A"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()

###### Smylie Myeloid

"CLEC10A", "CD1C", "S100A8", "S100A9", "TPSAB1", "OSM", "NKG7", "KLRC1"

png("cell_SMILLIE_MYELOID_marker_res_0.8.png", units="in", width=40, height=40, res=300)
FeaturePlot(alldata, reduction = "UMAP_on_CCA",features = c("CLEC10A", "CD1C", "S100A8", "S100A9", "TPSAB1", "OSM", "NKG7", "KLRC1"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()



png("cell_SMILLIE_MYELOID_marker_res_0.5.png", units="in", width=40, height=40, res=300)
FeaturePlot(alldata_0.5, reduction = "UMAP_on_CCA",features = c("CLEC10A", "CD1C", "S100A8", "S100A9", "TPSAB1", "OSM", "NKG7", "KLRC1"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()


png("cell_SMILLIE_MYELOID_marker_res_1.5.png", units="in", width=40, height=40, res=300)
FeaturePlot(alldata_1.5, reduction = "UMAP_on_CCA",features = c("CLEC10A", "CD1C", "S100A8", "S100A9", "TPSAB1", "OSM", "NKG7", "KLRC1"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()



png("cell_SMILLIE_MYELOID_marker_res_2.png", units="in", width=40, height=40, res=300)
FeaturePlot(alldata_2, reduction = "UMAP_on_CCA",features = c("CLEC10A", "CD1C", "S100A8", "S100A9", "TPSAB1", "OSM", "NKG7", "KLRC1"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()



png("cell_SMILLIE_MYELOID_marker_res_2.5.png", units="in", width=40, height=40, res=300)
FeaturePlot(alldata_2.5, reduction = "UMAP_on_CCA",features = c("CLEC10A", "CD1C", "S100A8", "S100A9", "TPSAB1", "OSM", "NKG7", "KLRC1"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()


png("cell_SMILLIE_MYELOID_marker_res_3.png", units="in", width=40, height=40, res=300)
FeaturePlot(alldata_3, reduction = "UMAP_on_CCA",features = c("CLEC10A", "CD1C", "S100A8", "S100A9", "TPSAB1", "OSM", "NKG7", "KLRC1"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()


png("cell_SMILLIE_MYELOID_marker_res_3.5.png", units="in", width=40, height=40, res=300)
FeaturePlot(alldata_3.5, reduction = "UMAP_on_CCA",features = c("CLEC10A", "CD1C", "S100A8", "S100A9", "TPSAB1", "OSM", "NKG7", "KLRC1"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()



png("cell_SMILLIE_MYELOID_marker_res_4.png", units="in", width=40, height=40, res=300)
FeaturePlot(alldata_4, reduction = "UMAP_on_CCA",features = c("CLEC10A", "CD1C", "S100A8", "S100A9", "TPSAB1", "OSM", "NKG7", "KLRC1"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()


###### T- cellS

"FOXP3", "CTLA4", "CD8B", "CXCR6", "CD3D"

png("cell_SMILLIE_T-CELLS_marker_res_0.8.png", units="in", width=40, height=40, res=300)
FeaturePlot(alldata, reduction = "UMAP_on_CCA",features = c("FOXP3", "CTLA4", "CD8B", "CXCR6", "CD3D"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()



png("cell_SMILLIE_T-CELLS_marker_res_0.5.png", units="in", width=40, height=40, res=300)
FeaturePlot(alldata_0.5, reduction = "UMAP_on_CCA",features = c("FOXP3", "CTLA4", "CD8B", "CXCR6", "CD3D"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()


png("cell_SMILLIE_T-CELLS_marker_res_1.5.png", units="in", width=40, height=40, res=300)
FeaturePlot(alldata_1.5, reduction = "UMAP_on_CCA",features = c("FOXP3", "CTLA4", "CD8B", "CXCR6", "CD3D"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()



png("cell_SMILLIE_T-CELLS_marker_res_2.png", units="in", width=40, height=40, res=300)
FeaturePlot(alldata_2, reduction = "UMAP_on_CCA",features = c("FOXP3", "CTLA4", "CD8B", "CXCR6", "CD3D"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()



png("cell_SMILLIE_T-CELLS_marker_res_2.5.png", units="in", width=40, height=40, res=300)
FeaturePlot(alldata_2.5, reduction = "UMAP_on_CCA",features = c("FOXP3", "CTLA4", "CD8B", "CXCR6", "CD3D"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()


png("cell_SMILLIE_T-CELLS_marker_res_3.png", units="in", width=40, height=40, res=300)
FeaturePlot(alldata_3, reduction = "UMAP_on_CCA",features = c("FOXP3", "CTLA4", "CD8B", "CXCR6", "CD3D"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()


png("cell_SMILLIE_T-CELLS_marker_res_3.5.png", units="in", width=40, height=40, res=300)
FeaturePlot(alldata_3.5, reduction = "UMAP_on_CCA",features = c("FOXP3", "CTLA4", "CD8B", "CXCR6", "CD3D"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()



png("cell_SMILLIE_T-CELLS_marker_res_4.png", units="in", width=40, height=40, res=300)
FeaturePlot(alldata_4, reduction = "UMAP_on_CCA",features = c("FOXP3", "CTLA4", "CD8B", "CXCR6", "CD3D"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()



######## B-Cells

"MZB1", "IGHA1", "SELL", "CD19", "AICDA"


png("cell_SMILLIE_B-CELLS_marker_res_0.8.png", units="in", width=40, height=40, res=300)
FeaturePlot(alldata, reduction = "UMAP_on_CCA",features = c("MZB1", "IGHA1", "SELL", "CD19", "AICDA"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()



png("cell_SMILLIE_B-CELLS_marker_res_0.5.png", units="in", width=40, height=40, res=300)
FeaturePlot(alldata_0.5, reduction = "UMAP_on_CCA",features = c("MZB1", "IGHA1", "SELL", "CD19", "AICDA"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()


png("cell_SMILLIE_B-CELLS_marker_res_1.5.png", units="in", width=40, height=40, res=300)
FeaturePlot(alldata_1.5, reduction = "UMAP_on_CCA",features = c("MZB1", "IGHA1", "SELL", "CD19", "AICDA"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()



png("cell_SMILLIE_B-CELLS_marker_res_2.png", units="in", width=40, height=40, res=300)
FeaturePlot(alldata_2, reduction = "UMAP_on_CCA",features = c("MZB1", "IGHA1", "SELL", "CD19", "AICDA"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()



png("cell_SMILLIE_B-CELLS_marker_res_2.5.png", units="in", width=40, height=40, res=300)
FeaturePlot(alldata_2.5, reduction = "UMAP_on_CCA",features = c("MZB1", "IGHA1", "SELL", "CD19", "AICDA"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()


png("cell_SMILLIE_B-CELLS_marker_res_3.png", units="in", width=40, height=40, res=300)
FeaturePlot(alldata_3, reduction = "UMAP_on_CCA",features = c("MZB1", "IGHA1", "SELL", "CD19", "AICDA"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()


png("cell_SMILLIE_B-CELLS_marker_res_3.5.png", units="in", width=40, height=40, res=300)
FeaturePlot(alldata_3.5, reduction = "UMAP_on_CCA",features = c("MZB1", "IGHA1", "SELL", "CD19", "AICDA"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()



png("cell_SMILLIE_B-CELLS_marker_res_4.png", units="in", width=40, height=40, res=300)
FeaturePlot(alldata_4, reduction = "UMAP_on_CCA",features = c("MZB1", "IGHA1", "SELL", "CD19", "AICDA"),min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()

##################

Converting Seurat object to LOG 2 - Log normalised count data


#################


new.data = log2(exp(as.matrix(GetAssayData(object = tum.sub, slot = "data"))))



###################
Merging two data frames 
###################

library (Seurat)


alldata <- readRDS("/data/All_samples_053120/2_Subclustring/1_Epithelial_subclust_results/1_Contamination_removed/1_Markers_Clusters/RDS/contamination_minus_Epi_cluster_res_0.8_060320.rds")


#divide the alldata into normal and tumor samples
# later we can convert those counts into log2value and use it to plot the heatmap

# Step1 --Dividing the alldata obj into Normal 




# we have 2738 Normal Cells and 6227 Tumor Epithelial cells

# lets take out only normal cells


n.epi<- subset(x = alldata, subset = Condition == "Normal")

# Step 2 --now lets convert the counts into log2



# Step 2 --now lets convert the counts into log2 (n.epi)


n.log2 = log2(exp(as.matrix(GetAssayData(object = n.epi, slot = "data"))))

marker_genes <- read.csv("/data/All_samples_053120/2_Subclustring/1_Epithelial_subclust_results/1_Contamination_removed/5_Normal_log2/Marker_gene_list.csv", sep=",", header = TRUE)
# make the First col of Marker_genes as rownames

# convert the first column as row name
# https://stackoverflow.com/questions/5555408/convert-the-values-in-a-column-into-row-names-in-an-existing-data-frame-in-r


 rownames(marker_genes) <- marker_genes[,1]
#head(marker_genes)

marker_genes[,1] <- NULL
#head(marker_genes)


# now merge the markergenes with n.log2

n.merge <- merge(marker_genes, n.log2, by=0)
dim(n.merge)


write.csv(n.merge,"Normal_Log2Values_Marker_genes.csv",row.names=TRUE)

#######

rownames(counts) = sapply(rownames(counts),function(v) paste('XXXX',v,'YYYY',sep='_')) 
data_obj = dataConstruct(counts);


######


######################
Running GSVA analysis

######################


alldata <- readRDS("/data/All_samples_053120/2_Subclustring/1_Epithelial_subclust_results/1_Contamination_removed/1_Markers_Clusters/RDS/contamination_minus_Epi_cluster_res_0.8_060320.rds")
meta <- read.csv("Tumor_Nor_conditions_metadata.csv", row.names=1, header=TRUE, sep = ",")

Counts<-as.matrix(alldata@assays$RNA@counts)

org.Hs.eg.db

library(GSEABase)

H.gsea <- getGmt("h.all.v7.1.symbols.gmt", collectionType=BroadCollection(category="h"), geneIdType=SymbolIdentifier())


# Load package
library(Biobase)

# Create ExpressionSet object

eset <- ExpressionSet(assayData = x, phenoData = AnnotatedDataFrame(p), featureData = AnnotatedDataFrame(f))


filtered_eset <- nsFilter(eset, require.entrez=TRUE, remove.dupEntrez=TRUE, var.func=IQR, var.filter=TRUE, var.cutoff=0.5, filterByQuantile=TRUE, feature.exclude="^AFFX")


cache(leukemia_es <- gsva(filtered_eset, H.gsea, min.sz=10, max.sz=500, verbose=TRUE)$es.obs,dir=cacheDir, prefix=cachePrefix)

cache(es <- gsva(filtered.eset, H.gsea, min.sz=10, max.sz=500, verbose=TRUE)$es.obs, dir=cacheDir, prefix=cachePrefix)



####

maits_sce <- createSCE(assayFile = t(maits$expressionmat),
                       annotFile = maits$cdat,
                       featureFile = maits$fdat,
                       assayName = "logtpm",
                       inputDataFrames = TRUE,
                       createLogCounts = FALSE)


gsvaRes <- gsvaSCE(maits_entrez, useAssay = "logtpm", H.gsea, parallel.sz=1)


all_subset <- convertGeneIDs(all, inSymbol = "ENTREZID", outSymbol = "SYMBOL", database = "org.Hs.eg.db")

	
all_subset <- convertGeneIDs(all, inSymbol = "SYMBOL", outSymbol = "ENTREZID", database = "org.Hs.eg.db")


subsetDES <- DESeqDataSetFromMatrix(countData = subsetdat@assays$RNA@counts, colData = meta,  design = ~ 1)





# http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

subsetDES <- DESeq(subsetDES, sfType="poscounts", useT=TRUE, minmu=1e-6, minReplicatesForReplace=Inf)



suppressPackageStartupMessages(library(GSEABase))
suppressPackageStartupMessages(library(GSVAdata))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(genefilter))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(GSVA))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(DESeq2))


leukemia_es <- gsva(assay(rld), H.gsea,min.sz=10, max.sz=500, verbose=TRUE)


DEgeneSets <- topTable(fit, coef="NormalVsTumor", number=Inf, p.value=adjPvalueCutoff, adjust="BH")


###
#https://www.biostars.org/p/369221/

markers_genes_0.8 <- FindAllMarkers(epi_0.8, logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, min.diff.pct = 0.2, only.pos = TRUE, max.cells.per.ident = 50, assay = "RNA")
write.table(markers_genes_0.8, file="marker_genes_0.8.txt", sep="\t",append=F, row.names = F)
cbind(table(markers_genes_0.8$cluster))

# Subset cluster based on cluster numbers

subset(markers_genes_0.8,markers_genes_0.8$cluster == 0)
c1 <-subset(markers_genes_0.8,markers_genes_0.8$cluster == 0)
c2 <-subset(markers_genes_0.8,markers_genes_0.8$cluster == 1)
c3 <-subset(markers_genes_0.8,markers_genes_0.8$cluster == 2)
c4 <-subset(markers_genes_0.8,markers_genes_0.8$cluster == 3)
c5 <-subset(markers_genes_0.8,markers_genes_0.8$cluster == 4)
c6 <-subset(markers_genes_0.8,markers_genes_0.8$cluster == 5)
c7 <-subset(markers_genes_0.8,markers_genes_0.8$cluster == 6)
c8 <-subset(markers_genes_0.8,markers_genes_0.8$cluster == 7)
c9 <-subset(markers_genes_0.8,markers_genes_0.8$cluster == 8)
c10 <-subset(markers_genes_0.8,markers_genes_0.8$cluster == 9)
c11 <-subset(markers_genes_0.8,markers_genes_0.8$cluster == 10)
c12 <-subset(markers_genes_0.8,markers_genes_0.8$cluster == 11)
c13 <-subset(markers_genes_0.8,markers_genes_0.8$cluster == 12)
c14 <-subset(markers_genes_0.8,markers_genes_0.8$cluster == 13)
c15 <-subset(markers_genes_0.8,markers_genes_0.8$cluster == 14)
c16 <-subset(markers_genes_0.8,markers_genes_0.8$cluster == 15)
c17 <-subset(markers_genes_0.8,markers_genes_0.8$cluster == 16)
c18 <-subset(markers_genes_0.8,markers_genes_0.8$cluster == 17)
c19 <-subset(markers_genes_0.8,markers_genes_0.8$cluster == 18)
c20 <-subset(markers_genes_0.8,markers_genes_0.8$cluster == 19)

c1g<-c1$gene
c2g<-c2$gene
c3g<-c3$gene
c4g<-c4$gene
c5g<-c5$gene
c6g<-c6$gene
c7g<-c7$gene
c8g<-c8$gene
c9g<-c9$gene
c10g<-c10$gene
c11g<-c11$gene
c12g<-c12$gene
c13g<-c13$gene
c14g<-c14$gene
c15g<-c15$gene
c16g<-c16$gene
c17g<-c17$gene
c18g<-c18$gene
c19g<-c19$gene
c20g<-c20$gene



#Get Cluster data

a1 <- SubsetData(alldata, subset.name = "seurat_clusters", accept.value = 0)

count1 <-as.matrix(a1@assays$RNA@counts)
count2 <-as.matrix(a2@assays$RNA@counts)
count3 <-as.matrix(a3@assays$RNA@counts)
count4 <-as.matrix(a4@assays$RNA@counts)
count5 <-as.matrix(a5@assays$RNA@counts)


#Subset Marker genes from the diff genes

mg1 <-subset(markers_genes_0.8,markers_genes_0.8$cluster == 2)
mg2 <-subset(markers_genes_0.8,markers_genes_0.8$cluster == 6)
mg3 <-subset(markers_genes_0.8,markers_genes_0.8$cluster == 14)
mg4 <-subset(markers_genes_0.8,markers_genes_0.8$cluster == 15)
mg5 <-subset(markers_genes_0.8,markers_genes_0.8$cluster == 18)

# Get the gene column and convert that into data frame

mg1 <- as.data.frame(mg1$gene)
mg2 <- as.data.frame(mg2$gene)
mg3 <- as.data.frame(mg3$gene)
mg4 <- as.data.frame(mg4$gene)
mg5 <- as.data.frame(mg5$gene)

# make rownames for mg1

rownames(mg1) <- mg1[,1]
rownames(mg2) <- mg2[,1]
rownames(mg3) <- mg3[,1]
rownames(mg4) <- mg4[,1]
rownames(mg5) <- mg5[,1]

# merge to get final counts

count1 <- merge(mg1,count1, by=0)
count2 <- merge(mg2,count2, by=0)
count3 <- merge(mg3,count3, by=0)
count4 <- merge(mg4,count4, by=0)
count5 <- merge(mg5,count5, by=0)


# allrownames to the counts

rownames(count1) <- count1[,2]
rownames(count2) <- count2[,2]
rownames(count3) <- count3[,2]
rownames(count4) <- count4[,2]
rownames(count5) <- count5[,2]


# delete Row.names and mg4$gene columns in all count1,2,3,4,5

count1[,1] <- NULL
count1[,1] <- NULL
count2[,1] <- NULL
count2[,1] <- NULL
count3[,1] <- NULL
count3[,1] <- NULL
count4[,1] <- NULL
count4[,1] <- NULL
count5[,1] <- NULL
count5[,1] <- NULL

# check the row names for the safer side





dif.g.1 <- merge(c1g,count1, by=0)

rownames(dif.g.1) <- dif.g.1$c1g

dif.g.1$c1g <- NULL
dif.g.1$Row.names <- NULL

# Metadata

meta.1 <- as.matrix(a1@meta.data)



write.csv(dif.g.1,"diff_m_c1.csv")




# 

exampleSet <- ExpressionSet(assayData=exprs,
+ phenoData=phenoData,
+ experimentData=experimentData,
+ annotation="hgu95av2")


#

filtered_eset <- nsFilter(eset, require.entrez=FALSE, remove.dupEntrez=TRUE, var.func=IQR, var.filter=TRUE, var.cutoff=0.5, filterByQuantile=TRUE, feature.exclude="^AFFX")



design <- model.matrix(~ factor(meta_clusters$seurat_clusters))



##########



colVars <- list(Clusters=c("0"="forestgreen", 
                           "1"="darkorange", 
                           "2"="magenta4", 
                           "3"="hotpink", 
                           "4"="red3", 
                           "5"="skyblue", 
                           "6"="darkblue",
			   "7"="brown",
			   "8"="darkgreen",
			   "9"="red",
			   "10"="yellow",	
			   "11"="lightblue",
			   "12"="pink",
			   "13"="lightgreen",
			   "14"="cyan",
			   "15"="violet",
			   "16"="#ff6361"))



genesKept <- geneFiltering(expr, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(expr),
                           minSamples=ncol(expr)*.01)


Maximum value in the expression matrix: 6503
Ratio of detected vs non-detected: 0.047
Number of counts (in the dataset units) per gene:
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
      0       5      89    1213     546  887702 
Number of cells in which each gene is detected:
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    0.0     5.0    82.5   404.6   482.0  8511.0 

Number of genes left after applying the following filters (sequential):
	8831	genes with counts per gene 268.95
	8819	genes detected in more than 89.65 cells
	7849	genes available in RcisTarget database
Gene list saved in int/1.1_genesKept.Rds

############################################################################################################################################################

# USE Screen

screen
nohup Rscript SCENIC_right_left.R &
Jobs -l
first press control a and then d to detach the screen
screen -list # to check the screen active
screen -r # to resume the terminal

screen -X -S SCREENID kill# to kill the screen 
pkill screen # to kill all screen




################################################################################################################################################################################

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=system.file("extdata", "oligodendroglioma_expression_downsampled.counts.matrix.gz", package = "infercnv"),
                                    annotations_file=system.file("extdata", "oligodendroglioma_annotations_downsampled.txt", package = "infercnv"),
                                    delim="\t",
                                    gene_order_file=system.file("extdata", "gencode_downsampled.EXAMPLE_ONLY_DONT_REUSE.txt", package = "infercnv"),
                                    ref_group_names=c("Microglia/Macrophage","Oligodendrocytes (non-malignant)")) 



infercnv_obj = CreateInfercnvObject( raw_counts_matrix = "exp_raw.txt",
  annotations_file = "meta.txt",
  gene_order_file="genes_code_mod.txt",
  delim="\t",
  ref_group_names=c("Normal-Left","Normal-Right"))


"/Users/akhaliq/Desktop/conics/normal_tumor"
# for group names

infercnv_obj.all = CreateInfercnvObject(raw_counts_matrix=epi@assays$RNA@counts, annotations_file="all_epi_meta.txt", delim="\t", gene_order_file="genes_code_mod.txt", ref_group_names="Normal")
infercnv_obj.all = infercnv::run(infercnv_obj.all,cutoff=0.1, out_dir="/Users/akhaliq/Desktop/conics/normal_tumor/",plot_steps=F, HMM=TRUE, denoise=TRUE, cluster_by_groups=TRUE, num_threads=3, no_prelim_plot=TRUE,png_res=300)



#
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=epi.subset@assays$RNA@counts, annotations_file="epi_tumor_meta.txt", delim="\t", gene_order_file="genes_code_mod.txt", ref_group_names=NULL)


# for my system ()

infercnv_obj = infercnv::run(infercnv_obj,cutoff=0.1, out_dir="/Users/akhaliq/Desktop/conics/",plot_steps=F, HMM=T, denoise=T, cluster_by_groups=T, num_threads=3, no_prelim_plot=TRUE,png_res=300)

 
# for server
infercnv_obj = infercnv::run(infercnv_obj,cutoff=0.1, out_dir="/data/infercnv/",plot_steps=T,HMM=F, denoise=T,cluster_by_groups=T, analysis_mode="subclusters", num_threads=2,png_res=300)



########

library(inferCNV)
infercnv_obj = CreateInfercnvObject(raw_counts_matrix= epi.subset@assays$RNA@counts,
                                    annotations_file= "epi_tumor_meta.txt",
                                    delim="\t",
                                    gene_order_file= "inferCNV_gene_order_file.txt",
                                    ref_group_names=NULL)

# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq/Fluidigm C1, 0.1 for 10x-genomics/SNRS/MARS-Seq/Microwell/Drop-seq
                             out_dir="Epi_Tumor1",
			     min_cells_per_gene = 50,  # dir is auto-created for storing outputs
                             cluster_by_groups=FALSE,   # cluster
                             HMM=TRUE,
			     k_obs_groups = 2,
                             denoise=TRUE
		             )


seurat_obj = infercnv::add_to_seurat(infercnv_output_path="/Users/akhaliq/Desktop/infercnv/Epi_Tumor1/",
                                     seurat_obj=epi.subset, # optional
                                     top_n=10
                                     )

####

seurat_obj = infercnv::add_to_seurat(infercnv_output_path=/Users/akhaliq/Desktop/conics/, seurat_obj=epi.subset, top_n=10)

#

out_dir = tempfile()
infercnv_obj_default = infercnv::run(
    infercnv_obj,
    cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
    out_dir="/data/All_samples_053120/2_Subclustring/1_Epithelial_subclust_results/1_Contamination_removed/2_Normal_vs_Tumor/infercnv_normal_tumor/output/",
    cluster_by_groups=TRUE, 
    plot_steps=FALSE,
    denoise=TRUE,
    HMM=FALSE,
    no_prelim_plot=TRUE,
    png_res=300
)

###

infercnv_obj = CreateInfercnvObject( raw_counts_matrix = kor.epi@assays$RNA@counts,
  annotations_file = "kor.epi.txt",
  gene_order_file="genes_code_mod.txt",
  delim="\t",
  ref_group_names=c("Normal-left-MSS","Normal-right-MSI-H","Normal-right-MSS"))

out_dir = tempfile()
infercnv_obj_default = infercnv::run(
    infercnv_obj,
    cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
    out_dir="/Users/akhaliq/Desktop/scrna_analysis/Korean_data/output/",
    cluster_by_groups=TRUE, 
    plot_steps=FALSE,
    denoise=TRUE,
    HMM=TRUE,
    no_prelim_plot=TRUE,
    png_res=300
)


/Users/akhaliq/Desktop/scrna_analysis/Korean_data/output
#####################################################################################################################################################################

Downloding File through R -WGET

####################################################################################################################################################################################

download.file("https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/7.1/h.all.v7.1.symbols.gmt", destfile = "h.all.v7.1.symbols.gmt", method = "wget")

########################################################################################################################################################################################################

###
Selecting only one or listed columns

meta <- subset(meta, select = Condition) # for only one Column
meta <- subset(meta, select = c(orig.ident, nCount_RNA, nFeature_RNA, Condition, Location) # for set of coloumns



###########

GSVA


###########



library(GSEABase)
library(GSVAdata)
library(Biobase)
library(genefilter)
library(limma)
library(RColorBrewer)
library(GSVA)
library(DESeq2)



setwd("/data/All_samples_053120/2_Subclustring/2_Stroma_subclust/2.1_Fibroblast_subclust/2.1.1_Fibro_removal_MSiH/GSVA")
H.gsea <- getGmt("/Users/akhaliq/Desktop/scrna_analysis/misc/ppt_ku/h.all.v7.1.symbols.gmt", collectionType=BroadCollection(category="h"), geneIdType=SymbolIdentifier())



meta.left <- mss.left_0.8@meta.data
meta.left <- subset(meta.left, select = Condition)
head(meta.left)

subsetDES.left <- DESeqDataSetFromMatrix(countData = mss.left_0.8@assays$RNA@counts, colData = meta.left,  design = ~ 1)
subsetDES.left <- DESeq(subsetDES.left, sfType="poscounts", useT=TRUE, minmu=1e-6, minReplicatesForReplace=Inf)

rld.left <- varianceStabilizingTransformation(subsetDES.left)
#ls()
topMatrixGSVA.left <- gsva(assay(rld.left), H.gsea, min.sz=10, max.sz=999999, abs.ranking=FALSE, verbose=TRUE)
##ls()
#head(topMatrixGSVA.right)
#ls()
merge.left <-  merge(as.data.frame(mss.left_0.8@meta.data), t(topMatrixGSVA.left), by = 0)

install.packages('rjags',configure.args='--with-jags-prefix=/home/masoodlab/JAGS-4.3.0 --with-jags-includedir=/home/emontene/JAGS-4.3.0/include/ --with-jags-libdir=/home/emontene/JAGS-4.3.0/lib/ --enable-rpath')



#head(merge.right)
dim(merge.left)
dim(topMatrixGSVA.left)

write.csv(merge.left, "GSVA_MSS_left_Fibro.csv")
#ls()
save.image("GSVA_MSS_LEFT.RData")
savehistory("GSVA_MSS_LEFT.RHistory")


######### 
Plotting Marker Genes
#########


png("Marker_genes_right.png", units="in", width=40, height=40, res=300)

mypar(1, 5, mar = c(4, 6, 3, 1))
for (i in unique(markers_genes_mss_right_0.8$cluster)) {
    barplot(sort(setNames(markers_genes_mss_right_0.8$avg_logFC, markers_genes_mss_right_0.8$gene)[markers_genes_mss_right_0.8$cluster == i], F), horiz = T, las = 1, main = paste0(i, " vs. rest"), border = "white", yaxs = "i")
    abline(v = c(0, 0.25), lty = c(1, 2))
}

dev.off()



top25_left <- markers_genes_mss_left_0.8 %>% group_by(cluster) %>% top_n(-25, p_val_adj)

mss.left_0.8.heatmap <- ScaleData(mss.left_0.8, features = as.character(unique(top25_left$gene)), assay = "RNA")

png("Marker_genes_Left_heatmap.png", units="in", width=20, height=20, res=300)
DoHeatmap(mss.left_0.8.heatmap, features = as.character(unique(top25_left$gene)), assay = "RNA")
dev.off()

############
Seurat clusters from 1 
############


seurat_object$seurat_clusters <- as.factor(as.numeric(as.character(seurat_object$seurat_clusters)) + 1)
Idents(alldata) <- "seurat_clusters"

png("UMAP_TSNE_PCA_Fibroblast_Right.png", units="in", width=10, height=10, res=300)
DimPlot(alldata, reduction = "UMAP_on_CCA", label=T, group.by="seurat_clusters")
dev.off()

########

########################################################################################################################################################################

#####################. STACKED VIOLIN PLOT

Ref : https://divingintogeneticsandgenomics.rbind.io/post/stacked-violin-plot-for-visualizing-single-cell-data-in-seurat/

library(Seurat)
library(patchwork)
library(ggplot2)

modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "mm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(2), angle = 0, face= "bold"), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "mm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(size = rel(2),angle = 90, face= "bold"), axis.ticks.x = element_line()) # for 90 degree naming use this code instead theme(axis.text.x=element_text(angle = 90, face= "bold"), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))

  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}



png("Stacked_violin_plot_Fibro_SPATIAL_new_cafs4_NEW.png", units="in", width=10, height=20, res=300)
StackedVlnPlot(obj = fibro, features = features_spatial_vcaf)+ theme(plot.title = element_text(size=5))
dev.off()




# CD4 #features<- c("IFNG","TOX", "FOXP3","IL2RA","TIGIT","TNFRSF4","TNFRSF9","TNFRSF18", "CD27", "CTLA4", "CD69", "GZMB", "PRF1","NKG7", "IL17A", "STMN1", "TUBB", "PCNA", "HMGB1", "HMGB2")
#CD8


#features = c("HAVCR2","LAG3","ENTPD1","FGFBP2","SELL","KLRB1","TRAV2","XCL1","XCL2")

features <- c("TNFRSF4", "FOXP3", "TNFRSF18", "CTLA4", "IL2RA", "TIGIT", "ANXA1", "IL7R", "CCR7", "SELL", "GPR18", "GZMK", "GZMH", "CCL4L2", "CST7", "XCL1", "XCL2", "RPS26", "MALAT1", "HSPA6", "ZFAND2A", "HSPA1A", "HSPA1B", "DNAJB4", "BAG3", "KLF2","IL17A", "TRAV1-2", "GZMB", "IL22", "CXCL13", "TOX2", "ICA1", "GZMA", "FGFBP2", "NKG7", "GNLY","FCGR3A", "PRF1", "CCL3", "CCL4", "CXCL13")

png("6_Stacked_violin_plot_Tcells.png", units="in", width=15, height=20, res=300)
StackedVlnPlot(obj = tcell, features = features)
dev.off()



#####################


Renaming Seurat Clusters

#####################
current.cluster.ids <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23)

new.cluster.ids <- c("CD4-TREG", "CD4-ANXA1-IL7R", "CD4-CM", "CD8-GZMK", "CD8-CHEMOKINES", "CD4-CLUSTER-6", "CD4-CM", "CD4-CM","CD4-HSP", "CD4-TH17", "CD4-TH17", "CD4-CXCL13", "CD8-GZMB", "CD8-GZMA", "CD8-GZMK", "NK-CELLS", "CD4-CM","CD8-CM","NK-CELLS","CD8-CXCL13","CD8-CHEMOKINES","CD8-CM")


new.cluster.ids <- c("CD4-CM", "CD4-TH17", "CD4-TREG", "CD4-CM", "CD8-Naive", "CD8-CM", "CD8-Chemokines", "CD8-Effector","CD8-TumorReactive", "CD4-TREG")

tcell$seurat_clusters <- plyr::mapvalues(x = tcell$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)

png("UMAP_TCELLS_Newcelltypes.png", units="in", width=10, height=10, res=300)
DimPlot(tcells, reduction = "umap", label=T, group.by="seurat_clusters",label.size =6)
dev.off()


new.cluster.ids <- c("CD4-TREG", "CD4-ANXA1-IL7R", "CD4-CM", "CD8-GZMK", "CD8-CHEMOKINES", "CD4-CLUSTER 6", "CD4-CM", "CD4-CM","CD4-HSP", "CD4-TH17", "CD4-TH17", "CD4-CXCL13", "CD8-GZMB", "CD8-GZMA", "CD8-GZMK", "NK-CELLS", "CD4-CM","CD8-CM","NK-CELLS","CD8-CXCL13","CD8-CHEMOKINES","CD8-S100B")

png("2_UMAP_TCELLS_celltypes.png", units="in", width=10, height=10, res=300)
dittoDimPlot(tcell, "ident",do.label = TRUE, labels.repel = FALSE)
dev.off()


png("3_UMAP_TCELLS_celltypes.png", units="in", width=10, height=10, res=300)
dittoDimPlot(tcell, "seurat_clusters",do.label = TRUE, labels.repel = FALSE)
dev.off()


png("4_TCELLS_BarPlot.png", units="in", width=10, height=10, res=300)
dittoBarPlot(tcells, "Condition", group.by = "seurat_clusters")
dev.off()

png("4.1_TCELLS_BarPlot.png", units="in", width=10, height=10, res=300)
dittoBarPlot(tcell, "Location", group.by = "seurat_clusters")
dev.off()

png("4.2_TCELLS_BarPlot.png", units="in", width=10, height=10, res=300)
dittoBarPlot(tcell, "MSI_Status", group.by = "seurat_clusters")
dev.off()

png("5_TCELLS_dimPlot.png", units="in", width=10, height=10, res=300)
dittoDimPlot(tcell, "Condition", split.by = "seurat_clusters")
dev.off()

png("5.1_TCELLS_dimPlot.png", units="in", width=10, height=10, res=300)
dittoDimPlot(tcell, "Location", split.by = "seurat_clusters")
dev.off()

png("5.2_TCELLS_dimPlot.png", units="in", width=10, height=10, res=300)
dittoDimPlot(tcell, "MSI_Status", split.by = "seurat_clusters")
dev.off()

#########

#####################################################################################################################


 Seurat subsetting ,Clustering and marker genes identification

#####################################################################################################################



library(Seurat)
library(dbplyr)

tcell <- readRDS("/data/All_samples_053120/2_Subclustring/3_T_cell/T-cells_res-0.8.rds")

ct26 <- SubsetData(tcell, subset.name = "seurat_clusters", accept.value = c(8,15,18,22,19,17))



ct26 <- NormalizeData(ct26, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
ct26 <- FindVariableFeatures(ct26, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
ct26 <- ScaleData(ct26, verbose = FALSE)

ct26 <- RunPCA(ct26, verbose = FALSE, features = VariableFeatures(object = ct26),npcs=35)
ct26 <- RunTSNE(object = ct26, dims.use = 1:35)
ct26 <- RunUMAP(ct26, dims = 1:35)
ct26 <- FindNeighbors(ct26, dims = 1:35, force.recalc = T)

ct26_allres <- FindClusters(ct26, resolution = c(0,0.2,0.5,0.8,1.5,2,2.5,3,3.5,4) , algorithm = 1, force.recalc=T)

library(clustree)
png("ct26_Clusttree.png", units="in", width=10, height=15, res=300)
clustree(ct26_allres, prefix = "SCT_snn_res.")
dev.off()

ct26 <- FindClusters(ct26, resolution = 0.8, algorithm = 1, force.recalc=T)
ct26$seurat_clusters <- as.factor(as.numeric(as.character(ct26$seurat_clusters)) + 1)
Idents(ct26) <- "seurat_clusters"

table(ct26$seurat_clusters)
png("UMAP_Orig_ident_ct26.png", units="in", width=10, height=10, res=300)
DimPlot(ct26, reduction = "umap", label=F, group.by="seurat_clusters")
dev.off()

markers_genes_ct26 <- FindAllMarkers(ct26, logfc.threshold = 0.2, test.use = "wilcox", only.pos = TRUE,assay = "RNA")
write.table(markers_genes_ct26, file="marker_genes_ct26.txt", sep="\t",append=F, row.names = F)
cbind(table(markers_genes_ct26$cluster))

top10 <- markers_genes_ct26 %>% group_by(cluster) %>% top_n(10, avg_log2FC)

png("heatmap_top10_ct26.png", units="in", width=10, height=10, res=300)
DoHeatmap(object = ct26, top10$gene)
dev.off()
 

 use.pcs = 1:35

experiment.aggregate <- FindClusters(
    object = ct26, 
    reduction.type = "pca", 
    dims.use = use.pcs, 
    resolution = seq(0.5,4,0.5), 
    print.output = FALSE
)
PrintFindClustersParams(object = experiment.aggregate)


FeaturePlot( ct26, "Lgals7", cols.use = c("lightgrey", "blue") )




png("Alldata_markergenes_compartment.png", units="in", width=40, height=10, res=300)
cowplot::plot_grid(ncol = 6,
FeaturePlot(alldata,features = "Krt18" ,min.cutoff="q9", cols=c("lightgrey", "#02a7c0"), label=FALSE)+ labs(title = "Epithelial Cells"),
FeaturePlot(alldata,features = "Lgals7" ,min.cutoff="q9", cols=c("lightgrey", "#02a7c0"), label=FALSE)+ labs(title = "Tumor Cells"),
FeaturePlot(alldata,features = "Cd3d" ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE)+ labs(title = "T Cells"),
FeaturePlot(alldata,features = "Cd79a" ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE)+ labs(title = "B Cells"),
FeaturePlot(alldata,features = "Col1a1" ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE)+ labs(title = "Fibroblast"),
FeaturePlot(alldata,features = "Cldn5" ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE)+ labs(title = "Endothelial"))
dev.off()


genes = c("Krt8","Krt18","Epcam","Lgals7","Cd3d","Cd79a","Col1a1","Fap","Cldn5")
png("heatmap_ct26.png", units="in", width=10, height=10, res=300)
dittoHeatmap(ct26, genes, annot.by = c("seurat_clusters"), heatmap.colors = colorRampPalette(c("blue", "red"))(50))
dev.off()


## SingleR
# https://github.com/mousepixels/sanbomics_scripts/blob/main/single_r.Rmd

library(celldex)
library(SingleR)

ref <- ImmGenData()
results <- SingleR(test = as.SingleCellExperiment(ct26), ref = ref, labels = ref$label.main)
ct26$singlr_labels_Immgen <- results$labels

png("umap_immigen_singler_ct26.png", units="in", width=10, height=10, res=300)
DimPlot(ct26, reduction = 'umap', group.by = 'singlr_labels_Immgen', label = TRUE)
dev.off()

#mouseRNAseq
ref2 <- celldex::MouseRNAseqData()

results <- SingleR(test = as.SingleCellExperiment(ct26), ref = ref2, labels = ref2$label.main)

ct26$singlr_mouseRNA <- results$labels

png("umap_immigen_singler_ct26.png", units="in", width=20, height=10, res=300)
DimPlot(ct26, reduction = 'umap', group.by = 'singlr_mouseRNA', label = TRUE)
dev.off()

results <- SingleR(test = as.SingleCellExperiment(ct26), ref = ref2, labels = ref2$label.main, de.method="wilcox")
ct26$singlr_mouseRNA <- results$labels

png("heatmap_immigen_singler_ct26.png", units="in", width=10, height=15, res=300)
plotScoreHeatmap(results)
dev.off()

###
png("dimplot_immegen_ct263.png", units="in", width=20, height=15, res=300)
dittoDimPlot(ct26, "seurat_clusters", split.by = "singlr_labels_Immgen")
dev.off()

genes3<- c("Igfbp3", "Egfl7", "Igfbp7", "Fabp4", "Col1a2", "Dcn", "Sparc", "Col1a1", "Col3a1", "Cited1", "Lcn2", "Krt19", "Krt18", "Ptn", "Cd63", "Rbp1", "Cald1", "Lgals1", "Hspb1", "Gzmb", "AW112010", "Naaa", "H2-Eb1", "Irf8", "Cst3", "Ctla2a", "Xcl1", "Gzma", "Nkg7", "Ccl5", "Mafb", "Ms4a4c", "Ms4a6d", "Plac8", "Trac", "Trbc1", "Cd3d", "Cd3g", "Trbc2", "C1qc", "Lyz2", "Apoe", "C1qb", "C1qa", "Il1b", "Cxcl2", "G0s2", "S100a8", "S100a9", "H2-Aa", "Cd83", "Cd74", "Cd79a", "Igkc")

genes3 <- rev(genes3)

png("heatmap_immegen_ct263.png", units="in", width=20, height=15, res=300)
DoHeatmap(ct26, group.by="singlr_labels_Immgen",features = genes3) + NoLegend()
dev.off()

 Idents(ct26)<- "singlr_labels_Immgen"

fibro <- subset(ct26, ident=c("Stromal cells","Fibroblasts"))
exp.rawdata <- as.matrix(fibro@assays$RNA@counts)

library(copykat)
copykat.mouse <- copykat(rawmat=exp.rawdata, id.type="S", ngene.chr=5, win.size=25, KS.cut=0.1, sam.name="mouse_fibro", distance="euclidean", norm.cell.names="",output.seg="FLASE", plot.genes="TRUE", genome="mm10",n.cores=1)

##################

suppressPackageStartupMessages({
  library(Seurat)
  library(cowplot)
  library(ggplot2)
})

png("CD4_Central_memory_cells_0.8.png", units="in", width=10, height=10, res=300)
FeaturePlot(ct26, reduction = "umap",features = c("IFNG","TOX"), order = T, ncol = 3, min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()

png("CD4_Central_memory_cells_0.8.png", units="in", width=10, height=10, res=300)
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()


png("CD4_TREGS_cells_0.8.png", units="in", width=10, height=10, res=300)
FeaturePlot(ct26, reduction = "umap",features = c("FOXP3","IL2RA","TIGIT","TNFRSF4","TNFRSF9","TNFRSF18","CD27","CTLA4"), order = T, ncol = 3, min.cutoff="q9", cols=c("lightgrey", "green"), label=TRUE)
dev.off()

png("CD4_Activated_cells_0.8.png", units="in", width=10, height=10, res=300)
FeaturePlot(ct26, reduction = "umap",features = "CD69", order = T, ncol = 3, min.cutoff="q9", cols=c("lightgrey", "blue"), label=TRUE)
dev.off()


png("CD4_Cytotoxic_cells_0.8.png", units="in", width=10, height=10, res=300)
FeaturePlot(ct26, reduction = "umap",features = c("GZMB","PRF1","NKG7"), order = T, ncol = 3, min.cutoff="q9", cols=c("lightgrey", "magenta"), label=TRUE)
dev.off()

png("CD4_AntiTumor_cells_0.8.png", units="in", width=10, height=10, res=300)
FeaturePlot(ct26, reduction = "umap",features = "IL17A", order = T, ncol = 3, min.cutoff="q9", cols=c("lightgrey", "darkred"), label=TRUE)
dev.off()

png("CD4_Prolif_cells_0.8.png", units="in", width=10, height=10, res=300)
FeaturePlot(ct26, reduction = "umap",features = c("SRMN1","TUBB","PCNA","HMGB1","HMGB2"), order = T, ncol = 3, min.cutoff="q9", cols=c("lightgrey", "darkgreen"), label=TRUE)
dev.off()


png("CD4_Central_memory_cells_split_by condition_0.8.png", units="in", width=10, height=10, res=300)
FeaturePlot(ct26, reduction = "umap",features = c("IFNG","TOX"), order = T, ncol = 3, min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE,split.by="Condition")
dev.off()

png("CD4_CD8_0.8.png", units="in", width=10, height=10, res=300)
FeaturePlot(ct26, reduction = "umap",features = c("CD4","CD8A"), order = T, ncol = 3, min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE)
dev.off()

########

png("SingleR-T-cells_0.8_seurat_clusters.png", units="in", width=10, height=10, res=300)
SingleR.DrawHeatmap(singler.T.cells$singler[[1]]$SingleR.single, top.n = 50, clusters = ct26@meta.data$seurat_clusters, order.by.clusters=TRUE )
dev.off()



png("SingleR-T-cells_0.8__enocde_seurat_clusters.png", units="in", width=10, height=10, res=300)
SingleR.DrawHeatmap(singler.T.cells$singler[[2]]$SingleR.single, top.n = 50, clusters = ct26@meta.data$seurat_clusters, order.by.clusters=TRUE )
dev.off()


##############

Ditto Plot

##############

library(dittoSeq)
library(scRNAseq)
library(SingleCellExperiment)
library(Seurat)


genes <- c("TNFRSF4",
"FOXP3",
"TNFRSF18")

genes <- c("TNFRSF4", "FOXP3", "TNFRSF18", "CTLA4", "IL2RA", "TIGIT", "ANXA1", "IL7R", "CCR7", "SELL", "GPR18", "GZMK", "GZMH", "CCL4L2", "CST7", "XCL1", "XCL2", "RPS26", "MALAT1", "HSPA6", "ZFAND2A", "HSPA1A", "HSPA1B", "DNAJB4", "BAG3", "KLF2","IL17A", "TRAV1-2", "GZMB", "IL22", "CXCL13", "TOX2", "ICA1", "GZMA", "FGFBP2", "NKG7", "GNLY","FCGR3A", "PRF1", "CCL3", "CCL4")

png("8_T-cells_HEATMAP_seurat_clusters_marker.png", units="in", width=15, height=15, res=300)
dittoHeatmap(alldata, genes, annot.by = c("seurat_clusters", "Condition", "Location", "MSI_Status"), order.by = "seurat_clusters", cluster_cols = FALSE)
dev.off()

png("7_Dimplot_celltypes_marker.png", units="in", width=10, height=10, res=300)
dittoDimPlot(tells, "seurat_clusters",do.label = TRUE, labels.repel = FALSE)
dev.off()

png("10_T-cells_multiditto_cd4.png", units="in", width=15, height=15, res=300)
multi_dittoDimPlotVaryCells(tcell, "CD4", vary.cells.meta = "seurat_clusters")
dev.off()

#######
T test . 2 Groups

#######


library(tidyverse)
library(rstatix)
library(ggpubr)
   

cell.count <- read.csv("tcell_count.csv", header=TRUE, sep=',',row.names=1)

mydata.long <- cell.count %>% pivot_longer(-Group, names_to = "variables", values_to = "value")
mydata.long %>% sample_n(6)


stat.test <- mydata.long %>% group_by(variables) %>% t_test(value ~ Group) %>% adjust_pvalue(method = "BH") %>% add_significance()
write.csv(stat.test, "Tumor-normal-T-test.csv")

png("ttest_TN.png", units="in", width=10, height=10, res=300)

myplot <- ggboxplot(
  mydata.long, x = "Group", y = "value",
  fill = "Group", palette = "npg", legend = "none",
  ggtheme = theme_pubr(border = TRUE)
  ) +
  facet_wrap(~variables)
# Add statistical test p-values
stat.test <- stat.test %>% add_xy_position(x = "Group")
myplot + stat_pvalue_manual(stat.test, label = "p.adj.signif")

dev.off()

##multi group##

library(tidyverse)
library(rstatix)
library(ggpubr)
   

cell.count <- read.csv("tcell_count_multi.csv", header=TRUE, sep=',',row.names=1)

pwc <- cell.count %>% pairwise_t_test(CD4.TREG ~ Group, p.adjust.method = "bonferroni")

pwc <- cell.count %>% pairwise_t_test(Group ~ counts, p.adjust.method = "bonferroni")

mydata.long <- cell.count %>% pivot_longer(-Group, names_to = "variables", values_to = "value")
mydata.long %>% sample_n(6)


stat.test <- mydata.long %>% group_by(variables) %>% t_test(value ~ Group) %>% adjust_pvalue(method = "BH") %>% add_significance()
pwc <- mydata.long %>% group_by(variables) %>% pairwise_t_test(value ~ Group, p.adjust.method = "bonferroni")

png("ttest_TN.png", units="in", width=10, height=10, res=300)

myplot <- ggboxplot(
  mydata.long, x = "Group", y = "value",
  fill = "Group", palette = "npg", legend = "none",
  ggtheme = theme_pubr(border = TRUE)
  ) +
  facet_wrap(~variables)
# Add statistical test p-values
stat.test <- stat.test %>% add_xy_position(x = "Group")
myplot + stat_pvalue_manual(stat.test, label = "p.adj.signif")

dev.off()

##multi group##



###########


png("simple_perm_test.png", units="in", width=10, height=10, res=300)

(ggplot(cell.count,aes(MSS.Right,MSS.Left))
    + geom_boxplot()
    + stat_sum(colour="darkgray",alpha=0.5)
    + scale_size(breaks=1:2, range=c(3,6))
)

dev.off()



boxplot(aframe, horizontal=FALSE, varwidth=TRUE, notch=FALSE, range=1.5, outline=TRUE, names=c("A","B","C"), boxwex=0.3, border=c("blue"), col=c("red"))


#########. PERMUTATION T TEst 




library("ggplot2"); theme_set(theme_bw())
library("lmPerm")
library("coin")
library("gtools")
require(reshape2)

cell.count <-read.csv("mss-rt-lt.csv", header=TRUE, sep=',')

df.m <- melt(cell.count, id.var = "cell.type")


set.seed(101) ## for reproducibility
nsim <- 100000
res <- numeric(nsim)

#1-CD4-TREG

for (i in 1:nsim) {
    ## standard approach: scramble response value
    perm <- sample(nrow(df.m))
    bdat <- transform(df.m,value=value[perm])
    ## compute & store difference in means; store the value
    res[i] <- mean(bdat[bdat$cell.type =="CD4-TREG","value"])
}

obs <- mean(df.m[df.m$cell.type =="CD4-TREG","value"])

res <- c(res,obs)


png("ptest_mss_RT_LT_CD4-TREG.png", units="in", width=20, height=20, res=300)
hist(res,col="gray",las=1,main="")
abline(v=obs,col="red")
dev.off()

Pv.cd4treg <- mean(abs(res)>=abs(obs))  ## count both tails: matches lmPerm


#2-CD4-ANXA1-IL7R

for (i in 1:nsim) {
    ## standard approach: scramble response value
    perm <- sample(nrow(df.m))
    bdat <- transform(df.m,value=value[perm])
    ## compute & store difference in means; store the value
    res[i] <- mean(bdat[bdat$cell.type =="CD4-ANXA1-IL7R","value"])
}

obs <- mean(df.m[df.m$cell.type =="CD4-ANXA1-IL7R","value"])

res <- c(res,obs)


png("ptest_mss_RT_LT_CD4-ANXA1-IL7R.png", units="in", width=20, height=20, res=300)
hist(res,col="gray",las=1,main="")
abline(v=obs,col="red")
dev.off()

Pv.cd4anxa1 <- mean(abs(res)>=abs(obs))  ## count both tails: matches lmPerm



#3-CD4-CM  
for (i in 1:nsim) {
    ## standard approach: scramble response value
    perm <- sample(nrow(df.m))
    bdat <- transform(df.m,value=value[perm])
    ## compute & store difference in means; store the value
    res[i] <- mean(bdat[bdat$cell.type =="CD4-CM ","value"])
}

obs <- mean(df.m[df.m$cell.type =="CD4-CM ","value"])

res <- c(res,obs)


png("ptest_mss_RT_LT_CD4-CM.png", units="in", width=20, height=20, res=300)
hist(res,col="gray",las=1,main="")
abline(v=obs,col="red")
dev.off()

Pv.cd4cm <- mean(abs(res)>=abs(obs))  ## count both tails: matches lmPerm

#4-CD8-GZMK


for (i in 1:nsim) {
    ## standard approach: scramble response value
    perm <- sample(nrow(df.m))
    bdat <- transform(df.m,value=value[perm])
    ## compute & store difference in means; store the value
    res[i] <- mean(bdat[bdat$cell.type =="CD8-GZMK","value"])
}

obs <- mean(df.m[df.m$cell.type =="CD8-GZMK","value"])

res <- c(res,obs)


png("ptest_mss_RT_LT_CD8-GZMK.png", units="in", width=20, height=20, res=300)
hist(res,col="gray",las=1,main="")
abline(v=obs,col="red")
dev.off()

Pv.CD8GZMK <- mean(abs(res)>=abs(obs))


#5CD8-CHEMOKINES
for (i in 1:nsim) {
    ## standard approach: scramble response value
    perm <- sample(nrow(df.m))
    bdat <- transform(df.m,value=value[perm])
    ## compute & store difference in means; store the value
    res[i] <- mean(bdat[bdat$cell.type =="CD8-CHEMOKINES","value"])
}

obs <- mean(df.m[df.m$cell.type =="CD8-CHEMOKINES","value"])

res <- c(res,obs)


png("ptest_mss_RT_LT_CD8-CHEMOKINES.png", units="in", width=10, height=10, res=300)
hist(res,col="gray",las=1,main="")
abline(v=obs,col="red")
dev.off()

Pv.CD8CHEMOKINES <- mean(abs(res)>=abs(obs))


#6CD4-CLUSTER 6 

for (i in 1:nsim) {
    ## standard approach: scramble response value
    perm <- sample(nrow(df.m))
    bdat <- transform(df.m,value=value[perm])
    ## compute & store difference in means; store the value
    res[i] <- mean(bdat[bdat$cell.type =="CD4-CLUSTER 6","value"])
}

obs <- mean(df.m[df.m$cell.type =="CD4-CLUSTER 6","value"])

res <- c(res,obs)


png("ptest_mss_RT_LT_CD4-CLUSTER 6.png", units="in", width=10, height=10, res=300)
hist(res,col="gray",las=1,main="")
abline(v=obs,col="red")
dev.off()

Pv.cluster6 <- mean(abs(res)>=abs(obs))

#7CD4-HSP

for (i in 1:nsim) {
    ## standard approach: scramble response value
    perm <- sample(nrow(df.m))
    bdat <- transform(df.m,value=value[perm])
    ## compute & store difference in means; store the value
    res[i] <- mean(bdat[bdat$cell.type =="CD4-HSP","value"])
}

obs <- mean(df.m[df.m$cell.type =="CD4-HSP","value"])

res <- c(res,obs)


png("ptest_mss_RT_LT_CD4-HSP.png", units="in", width=20, height=20, res=300)
hist(res,col="gray",las=1,main="")
abline(v=obs,col="red")
dev.off()

Pv.CD4HSP <- mean(abs(res)>=abs(obs))


#8CD4-TH17
for (i in 1:nsim) {
    ## standard approach: scramble response value
    perm <- sample(nrow(df.m))
    bdat <- transform(df.m,value=value[perm])
    ## compute & store difference in means; store the value
    res[i] <- mean(bdat[bdat$cell.type =="CD4-TH17","value"])
}

obs <- mean(df.m[df.m$cell.type =="CD4-TH17","value"])

res <- c(res,obs)


png("ptest_mss_RT_LT_CD4-TH17.png", units="in", width=20, height=20, res=300)
hist(res,col="gray",las=1,main="")
abline(v=obs,col="red")
dev.off()

Pv.CD4TH17 <- mean(abs(res)>=abs(obs))


#9CD4-CXCL13
for (i in 1:nsim) {
    ## standard approach: scramble response value
    perm <- sample(nrow(df.m))
    bdat <- transform(df.m,value=value[perm])
    ## compute & store difference in means; store the value
    res[i] <- mean(bdat[bdat$cell.type =="CD4-CXCL13","value"])
}

obs <- mean(df.m[df.m$cell.type =="CD4-CXCL13","value"])

res <- c(res,obs)


png("ptest_mss_RT_LT_CD4-CXCL13.png", units="in", width=20, height=20, res=300)
hist(res,col="gray",las=1,main="")
abline(v=obs,col="red")
dev.off()

Pv.CD4CXCL13 <- mean(abs(res)>=abs(obs))



#10CD8-GZMB
for (i in 1:nsim) {
    ## standard approach: scramble response value
    perm <- sample(nrow(df.m))
    bdat <- transform(df.m,value=value[perm])
    ## compute & store difference in means; store the value
    res[i] <- mean(bdat[bdat$cell.type =="CD8-GZMB","value"])
}

obs <- mean(df.m[df.m$cell.type =="CD8-GZMB","value"])

res <- c(res,obs)


png("ptest_mss_RT_LT_CD8-GZMB.png", units="in", width=20, height=20, res=300)
hist(res,col="gray",las=1,main="")
abline(v=obs,col="red")
dev.off()

Pv.CD8GZMB <- mean(abs(res)>=abs(obs))



#11CD8-GZMA

for (i in 1:nsim) {
    ## standard approach: scramble response value
    perm <- sample(nrow(df.m))
    bdat <- transform(df.m,value=value[perm])
    ## compute & store difference in means; store the value
    res[i] <- mean(bdat[bdat$cell.type =="CD8-GZMA","value"])
}

obs <- mean(df.m[df.m$cell.type =="CD8-GZMA","value"])

res <- c(res,obs)


png("ptest_mss_RT_LT_CD8-GZMA.png", units="in", width=20, height=20, res=300)
hist(res,col="gray",las=1,main="")
abline(v=obs,col="red")
dev.off()

Pv.CD8GZMA <- mean(abs(res)>=abs(obs))


#12NK-CELLS

for (i in 1:nsim) {
    ## standard approach: scramble response value
    perm <- sample(nrow(df.m))
    bdat <- transform(df.m,value=value[perm])
    ## compute & store difference in means; store the value
    res[i] <- mean(bdat[bdat$cell.type =="NK-CELLS","value"])
}

obs <- mean(df.m[df.m$cell.type =="NK-CELLS","value"])

res <- c(res,obs)


png("ptest_mss_RT_LT_NK-CELLS.png", units="in", width=20, height=20, res=300)
hist(res,col="gray",las=1,main="")
abline(v=obs,col="red")
dev.off()

Pv.NKCELLS <- mean(abs(res)>=abs(obs))


#13CD8-CM


for (i in 1:nsim) {
    ## standard approach: scramble response value
    perm <- sample(nrow(df.m))
    bdat <- transform(df.m,value=value[perm])
    ## compute & store difference in means; store the value
    res[i] <- mean(bdat[bdat$cell.type =="CD8-CM","value"])
}

obs <- mean(df.m[df.m$cell.type =="CD8-CM","value"])

res <- c(res,obs)


png("ptest_mss_RT_LT_CD8-CM.png", units="in", width=20, height=20, res=300)
hist(res,col="gray",las=1,main="")
abline(v=obs,col="red")
dev.off()

Pv.CD8CM <- mean(abs(res)>=abs(obs))


#14CD8-CXCL13


for (i in 1:nsim) {
    ## standard approach: scramble response value
    perm <- sample(nrow(df.m))
    bdat <- transform(df.m,value=value[perm])
    ## compute & store difference in means; store the value
    res[i] <- mean(bdat[bdat$cell.type =="CD8-CXCL13","value"])
}

obs <- mean(df.m[df.m$cell.type =="CD8-CXCL13","value"])

res <- c(res,obs)


png("ptest_mss_RT_LT_CD8-CXCL13.png", units="in", width=20, height=20, res=300)
hist(res,col="gray",las=1,main="")
abline(v=obs,col="red")
dev.off()

Pv.CD8CXCL13 <- mean(abs(res)>=abs(obs))


#15CD8-S100B


for (i in 1:nsim) {
    ## standard approach: scramble response value
    perm <- sample(nrow(df.m))
    bdat <- transform(df.m,value=value[perm])
    ## compute & store difference in means; store the value
    res[i] <- mean(bdat[bdat$cell.type =="CD8-S100B","value"])
}

obs <- mean(df.m[df.m$cell.type =="CD8-S100B","value"])

res <- c(res,obs)


png("ptest_mss_RT_LT_CD8-S100B.png", units="in", width=20, height=20, res=300)
hist(res,col="gray",las=1,main="")
abline(v=obs,col="red")
dev.off()

Pv.CD8S100B <- mean(abs(res)>=abs(obs))

####
# Gini Index

https://www.r-bloggers.com/gini-index-and-lorenz-curve-with-r/
ineq(c(284,288),type="Gini")


#####
png("tsne_all_Korean_Data.png", units="in", width=30, height=30, res=300)
cowplot::plot_grid(ncol = 4,dittoDimPlot(alldata, "seurat_clusters",do.label = TRUE, labels.repel = FALSE,reduction="tsne"),
dittoDimPlot(alldata, "Cell_type",do.label = TRUE, labels.repel = FALSE,reduction="tsne"),
dittoDimPlot(alldata, "Condition",do.label = TRUE, labels.repel = FALSE,reduction="tsne"),
dittoDimPlot(alldata, "Location",do.label = TRUE, labels.repel = FALSE,reduction="tsne"),
dittoDimPlot(alldata, "MSI_Status",do.label = TRUE, labels.repel = FALSE,reduction="tsne"),
dittoDimPlot(alldata, "Cell_subtype",do.label = TRUE, labels.repel = FALSE,reduction="tsne")
)



######################################################################################################################################
Querying Multiple Columns and to get counts
######################################################################################################################################



library(dplyr)

meta %>% filter(Condition == "Tumor" & MSI_Status == "MSS" & Location == "Left") %>% count("orig.ident")







library(mAPKL)

###################

Convert symbols into Entrez genes

##################

library(CMScaller)
log2_entrez <- CMScaller::replaceGeneId(log2_all, id.in='symbol', id.out="entrez")
log2_entrez <- log2_entrez[grepl("^NA", rownames(log2_entrez))==F,]# remove those NA if required



#####################

ADDing METADATA - New COLUMN to the existing Seurat object

######################
#ref https://stackoverflow.com/questions/42279766/add-metadata-to-seurat-object

library(Seurat)
library(dplyr)
library(knitr)
library(janitor)

# Read the CMS prediction file
cms.meta <- read.csv("/Users/akhaliq/Desktop/scrna_analysis/Subclustring_crc_051520_Final/CMScaller/CMS_classification_V2.csv",header=TRUE,sep=',', row.names=1)
cms.bulk <- read.csv("/Users/akhaliq/Desktop/All_Clusters/bulk_cms_meta_all.csv",header=TRUE, sep=',', row.names=1)
# read the Meta data into separate object

fibro.meta <- fibro@meta.data

# merge

merge.meta <- merge(fibro.meta,cms.meta,by=0)

#Make the first column as row name

rownames(merge.meta) <- merge.meta[,1]

merge.meta[,1] <- NULL


#keep only the predictions column 
merge1 <- merge.meta[,"prediction", drop=FALSE]

# add Meta data to the object
fibro2 <- AddMetaData(fibro, merge1)

# check
fibro2@meta.data%>% count(orig.ident,prediction,sort=TRUE)

# to check the total at the end

kable(epi2@meta.data%>% count(prediction,sort=TRUE)%>% adorn_totals("row"))
kable(epi2@meta.data%>% count(MSI_Status,Location,Condition,prediction,sort=TRUE)%>% adorn_totals("row"))


###### CMS-Bulk MRNA

cms.bulk <- read.csv("/Users/akhaliq/Desktop/All_Clusters/bulk_cms_meta_all.csv",header=TRUE, sep=',', row.names=1)

merge.meta <- merge(tcells@meta.data,cms.bulk,by=0)
rownames(merge.meta) <- merge.meta[,1]
merge.meta[,1] <- NULL
merge1 <- merge.meta[,"bulk_prediction", drop=FALSE]
tcells <- AddMetaData(tcells, merge1)
tcells@meta.data%>% count(orig.ident,bulk_prediction,sort=TRUE)

kable(tcells@meta.data%>% count(MSI_Status,Location,Condition,prediction, bulk_prediction,sort=TRUE)%>% adorn_totals("row"))


########################################################
CMAP  Drug study
########################################################

library("DrInsight")

#file is in
# /Users/akhaliq/Desktop/trajectory/epithelial/rankMatrix.txt

cmap.ref.profiles = get.cmap.ref(cmap.data.path = 'rankMatrix.txt', probe.to.genes = probe.to.genes, drug.info = drug.info)


counts <- as.matrix(epi2@assays$RNA@counts)
head(counts)
meta <- epi2@meta.data
head(meta)
meta <- meta[,"Condition", drop=FALSE]
head(meta)
library("DESeq2")
subsetDES <- DESeqDataSetFromMatrix(countData = counts, colData = meta,  design = ~ Condition)


res1 <- res[,"pvalue", drop=FALSE]

res2<-res1[complete.cases(res1), ]
dim(res2)
res2<-na.omit(res1)
dim(res2)
head(res2)
library(data.table)


res3 <-setDT(res2, keep.rownames = "geneSymbol")

colnames(res3)[2] <- "score"


drug.ident.res = drug.ident(query.data = res3, cmap.ref.profiles = cmap.ref.profiles, repurposing.unit = "treatment", connectivity = "negative")

drug.pvals%>%filter(FDR  <= 0.1)

data("pathway.PID")
pathway.PID
path.analysis.res = pathway.analysis(drug.ident.res = drug.ident.res, pathway.list =pathway.PID, pathway.list.path = NULL, drug.FDR.cutoff = 0.1)


#NOTE: we had to abort the analysis because of FDR >0.1
##################################



################
Writing multiple objects in separate sheets of Xls
################


library(xlsx)
write.xlsx(dataframe1, file="filename.xlsx", sheetName="sheet1", row.names=FALSE)
write.xlsx(dataframe2, file="filename.xlsx", sheetName="sheet2", append=TRUE, row.names=FALSE)



########
Graph - GGplots
#######
Ref : https://www.r-graph-gallery.com/48-grouped-barplot-with-ggplot2.html


state_counts <- read.csv("states_counts1.csv", sep=',', header=TRUE)

png("trajectory_fibro_counts.png", units="in", width=10, height=10, res=300)
cowplot::plot_grid(ncol = 2,
ggplot(state_counts, aes(fill=Condition, y=Counts, x=States)) + geom_bar(position="fill", stat="identity"),
ggplot(state_cluster_counts, aes(fill=Clusters, y=Counts, x=States)) + geom_bar(position="fill", stat="identity"),
ggplot(state_countsCMS, aes(fill=Prediction, y=Counts, x=States)) + geom_bar(position="fill", stat="identity")
)



cowplot::plot_grid(ncol = 2,
dittoDimPlot(epi,"CD8A",do.label = TRUE, labels.repel = FALSE),
dittoDimPlot(epi,"CD8B",do.label = TRUE, labels.repel = FALSE),
dittoDimPlot(epi,"CD3D",do.label = TRUE, labels.repel = FALSE),
dittoDimPlot(epi,"CD4",do.label = TRUE, labels.repel = FALSE),
dittoDimPlot(epi,"seurat_clusters",do.label = TRUE, labels.repel = FALSE)
)

cowplot::plot_grid(ncol = 2,
dittoDimPlot(tcells_new,"seurat_clusters",do.label = TRUE, labels.repel = FALSE),
dittoDimPlot(tcells_new,"MSI_Status",do.label = TRUE, labels.repel = FALSE)
)

* For multiple graphs



png("trajectory_fibro_all_Clusters.png", units="in", width=20, height=20, res=300)
ggplot(state_counts_all, aes(fill=Condition, y=Counts, x=States)) + geom_bar(position="fill", stat="identity") + scale_fill_viridis(discrete = T, option = "viridis") +
    ggtitle("Cell Counts in Each Cluster") +
    facet_wrap(~Cluster) +
    theme_ipsum() +
    theme(legend.position="bottom") +
    xlab("")
dev.off()


########
markers <- findMarkers(endo.sce, clusters = colData(endo.sce)$seurat_clusters,
lfc = 1.5, direction = 'up', log.p = TRUE, 
                       BPPARAM = BiocParallel::MulticoreParam(1))

genes <- lapply(markers, function(x) {
    rownames(x)[x$Top <= 20]
})



plotHeatmap(endo.sce, genes, colour_columns_by = "seurat_clusters", show_colnames = FALSE, clustering_method = 'ward.D2',fontsize_row = 6)


#####


if (all(names(x = mean.exp) == rownames(x = endo1.2NT@assays$RNA))) {
  cat("Cell names order match in 'mean.exp' and 'endo1.2NT@assays$RNA':\n", 
      "adding gene set mean expression values in 'endo1.2NT@meta.data$gene.set.score'")
  endo1.2NT@meta.data$gene.set.score <- mean.exp
}




#######


png("Relative_abundance_tcells_TumorNormal.png", units="in", width=20, height=15, res=300)

ggplot(ggdf) +
  geom_boxplot(aes(x = condition, y = proportion, color = condition, 
    fill = condition),  position = position_dodge(), alpha = 0.5, 
    outlier.color = NA) +
  geom_point(aes(x = condition, y = proportion, color = condition, 
    shape = patient_id), alpha = 0.8, position = position_jitterdodge()) +
  facet_wrap(~ cluster, scales = "free", nrow = 2) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
    strip.text = element_text(size = 6)) +
  scale_color_manual(values = color_conditions) +
  scale_fill_manual(values = color_conditions) +
  scale_shape_manual(values = c(16, 17, 8, 3, 12, 0, 1, 2,4,5,6,7,9,10,45,22,33,44,55,66,77,23,34))


dev.off()


Relative abundance of the 12 Clusters in each sample, in the TCells Compartment, represented with boxplots. Different colors are used for the two conditions: Normal and Tumor. Values for each patient are indicated with different shape.


png("Relative_abundance_tcells_MSI.png", units="in", width=20, height=15, res=300)

ggplot(ggdf) +
  geom_boxplot(aes(x = msi, y = proportion, color = msi, 
    fill = msi),  position = position_dodge(), alpha = 0.5, 
    outlier.color = NA) +
  geom_point(aes(x = msi, y = proportion, color = msi, 
    shape = patient_id), alpha = 0.8, position = position_jitterdodge()) +
  facet_wrap(~ cluster, scales = "free", nrow = 2) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
    strip.text = element_text(size = 6)) +
  scale_color_manual(values = color_conditions) +
  scale_fill_manual(values = color_conditions) +
  scale_shape_manual(values = c(16, 17, 8, 3, 12, 0, 1, 2,4,5,6,7,9,10,45,22,33,44,55,66,77,23,34))

dev.off()

png("Relative_abundance_tcells_MSI-normlTumor.png", units="in", width=20, height=15, res=300)
ggplot(ggdf) +
  geom_boxplot(aes(x = msi_condition, y = proportion, color = msi_condition, 
    fill = msi_condition),  position = position_dodge(), alpha = 0.5, 
    outlier.color = NA) +
  geom_point(aes(x = msi_condition, y = proportion, color = msi_condition, 
    shape = patient_id), alpha = 0.8, position = position_jitterdodge()) +
  facet_wrap(~ cluster, scales = "free", nrow = 2) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
    strip.text = element_text(size = 6)) +
  scale_color_manual(values = color_conditions) +
  scale_fill_manual(values = color_conditions) +
  scale_shape_manual(values = c(16, 17, 8, 3, 12, 0, 1, 2,4,5,6,7,9,10,45,22,33,44,55,66,77,23,34))

dev.off()

color_conditions <- c("#6A3D9A", "#FF7F00","#33a02c","#e31a1c")



######
#dotPlot
DotPlot(object = tcells, features = genes, assay="RNA")+guides(color = guide_colorbar(title = 'Scaled Average Expression')) + theme(axis.text.x = element_text(angle=90))+ scale_colour_gradient2(low = "#FF00FF", mid = "#000000", high = "#FFFF00")




png("Spearman_correlation_Tcells_Our_Korean_1.png", units="in", width=20, height=20, res=300)
cowplot::plot_grid(ncol = 2,
P1, P2
)
dev.off()

#######

Merge unequal dataframes and replace missing rows with 0

#######


bulk <- read.csv("blk_cms.csv", header=T, sep = ',')
require(plyr)
zz<-join(all.meta, bulk, type="left")
write.csv(zz,"bulk_cms_meta_all.csv")

###

Monocle 2

###
Ref : https://nbisweden.github.io/workshop-scRNAseq/oldlabs/monocle_analysis


suppressMessages(library(monocle))
suppressMessages(library(stringr))
suppressMessages(library(plyr))
suppressMessages(library(netbiov))
library(Seurat)
library(cowplot)

load("/Users/akhaliq/Desktop/trajectory/epithelial/Location/Trajectory_analysis_epi_Location.RData")
ls()
R <- seurat@assays[["RNA"]]@counts
M <- seurat@meta.data
num_cells <- apply(R,1,function(x) sum(x>1))
dim(num_cells)

#num_cells

genes <- data.frame(gene_short_name = rownames(R),num_cells_expressed=num_cells)
rownames(genes)<-rownames(R)
pd <- new("AnnotatedDataFrame", data = M)
fd <- new("AnnotatedDataFrame", data = genes)



cds1 <- newCellDataSet(as.matrix(R),
                       phenoData = pd,
                       featureData = fd,
                       lowerDetectionLimit=0.5,
                       expressionFamily=negbinomial.size())
expressed_genes <- row.names(subset(fData(cds1), num_cells_expressed >= 10))
expressed_genes
fData(cds1)
cds1 <- estimateSizeFactors(cds1)
cds1
savefile <- "monocle_de_genes_NEW12.Rdata"
if (file.exists(savefile)){
  load(savefile)
}else{
  diff_test_res <- differentialGeneTest(cds1[expressed_genes,],
                                      fullModelFormulaStr="~Condition")
  save(diff_test_res,file=savefile)
}
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
length(ordering_genes)

cds1 <- setOrderingFilter(cds1, ordering_genes)
cds1 <- reduceDimension(cds1, max_components=2)
#Now that the space is reduced, it???s time to order the cells using the orderCells function as shown below.
cds1 <- orderCells(cds1)

plot_cell_trajectory(cds1, color_by="Condition")
plot_cell_trajectory(cds1, color_by="seurat_clusters")
plot_cell_trajectory(cds1, color_by="Location")
png("trajectory_epi_normal_tumor.png", units="in", width=10, height=10, res=300)
plot_cell_trajectory(cds1, color_by="Condition")
dev.off()

png("Epithelial_Trajectory.png", units="in", width=20, height=20, res=300)
cowplot::plot_grid(ncol = 2,
plot_cell_trajectory(cds1, color_by="Condition"),
plot_cell_trajectory(cds1, color_by="State"),
plot_cell_trajectory(cds1, color_by="MSI_Status"),
plot_cell_trajectory(cds1, color_by="Location"),
plot_cell_trajectory(cds1, color_by="prediction"),
plot_cell_trajectory(cds1, color_by="bulk_prediction"))
dev.off()


png("Epithelial_Trajectory_perdictions.png", units="in", width=20, height=20, res=300)
cowplot::plot_grid(ncol = 2,
plot_cell_trajectory(cds1, color_by="prediction",cell_name_size=6),
plot_cell_trajectory(cds1, color_by="bulk_prediction",cell_name_size=6))
dev.off()

# for user specific colours

svg("trajectory_epi_bulk_cms.svg", width=5, height=5)
plot_cell_trajectory(cds1, color_by="bulk_prediction")+ scale_color_manual(values=c("#dcc134", "#008000", "#37c8ab","#ff5555"))
dev.off()

svg("trajectory_belgian_tumor.svg", width=5, height=5)
plot_cell_trajectory(cds1, color_by="Cell_subtype")+ scale_color_manual(values=c("#dcc134", "#008000", "blue","#ff5555"))
dev.off()



##########

HEatmap of expression pattern of genes in Seurat

##########

Ref: https://github.com/satijalab/seurat/issues/1589

genes <- c("CD37","CD40","CD28","CD274","TNFRSF4","CTLA4","ICOS","TNFRSF18","PDCD1LG2","CEACAM1","CD80","IL2","CD200R1","BTLA","PDCD1","LAG3","TIGIT","IL10","TNFRSF9","HAVCR2")
png("HeatmapofExpressionPatternsTargetedImmunotherapies.png", units="in", width=20, height=10, res=300)
DoHeatmap(tcells,feature= genes,group.by = "seurat_clusters",raster=F)+scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "YIOrRd")) ) + guides(color=FALSE)
dev.off()



#########

Merging Multiple Dataframes by now names

#########
Ref: https://stackoverflow.com/questions/22617593/merge-multiple-data-frames-by-row-names

merge.all <- function(x, ..., by = "row.names") {
  L <- list(...)
  for (i in seq_along(L)) {
    x <- merge(x, L[[i]], by = by)
    rownames(x) <- x$Row.names
    x$Row.names <- NULL
  }
  return(x)
}

merge.all(x,y,z)


#################

Plotting 

################
png("Condition_bcells_BarPlot.png", units="in", width=10, height=10, res=300)
dittoBarPlot(bcells, "Condition", group.by = "seurat_clusters",data.out = TRUE)
dev.off()

png("Location_bcells_BarPlot.png", units="in", width=10, height=10, res=300)
dittoBarPlot(bcells, "Location", group.by = "seurat_clusters")
dev.off()

png("MSI_Status_bcells_BarPlot.png", units="in", width=10, height=10, res=300)
dittoBarPlot(bcells, "MSI_Status", group.by = "seurat_clusters")
dev.off()


png("samples_bcells_BarPlot.png", units="in", width=10, height=10, res=300)
dittoBarPlot(bcells, "orig.ident", group.by = "seurat_clusters")
dev.off()


png("CMS_bcells_BarPlot.png", units="in", width=10, height=10, res=300)
dittoBarPlot(bcells, "prediction", group.by = "seurat_clusters")
dev.off()


png("Bulk_cms_Tcells_BarPlot.png", units="in", width=10, height=10, res=300)
dittoBarPlot(tcells_new, "bulk_prediction", group.by = "seurat_clusters")
dev.off()




###


png("Epithelial_No_transcripts.png", units="in", width=10, height=10, res=300)
dittoBoxPlot(epi, "nFeature_RNA",group.by = "seurat_clusters",jitter.size=0.0000005)
dev.off()


png("Endo_No_transcripts.png", units="in", width=10, height=10, res=300)
dittoBoxPlot(endo, "nFeature_RNA",group.by = "seurat_clusters",jitter.size=0.0000005)
dev.off()


png("fibro_No_transcripts.png", units="in", width=10, height=10, res=300)
dittoBoxPlot(fibro, "nFeature_RNA",group.by = "seurat_clusters",jitter.size=0.0000005)
dev.off()


png("Tcells_No_transcripts.png", units="in", width=10, height=10, res=300)
dittoBoxPlot(tcells, "nFeature_RNA",group.by = "seurat_clusters",jitter.size=0.0000005)
dev.off()


png("Bcells_No_transcripts.png", units="in", width=10, height=10, res=300)
dittoBoxPlot(bcells, "nFeature_RNA",group.by = "seurat_clusters",jitter.size=0.0000005)
dev.off()


png("myeloid_No_transcripts.png", units="in", width=10, height=10, res=300)
dittoBoxPlot(myeloid, "nFeature_RNA",group.by = "seurat_clusters",jitter.size=0.0000005)
dev.off()

###




png("Alldata_graphs.png", units="in", width=40, height=10, res=300)
cowplot::plot_grid(ncol = 7,
dittoDimPlot(alldata, "Condition",reduction.use="UMAP_on_CCA"),
dittoDimPlot(alldata, "Location",reduction.use="UMAP_on_CCA"),
dittoDimPlot(alldata, "MSI_Status",reduction.use="UMAP_on_CCA"),
dittoDimPlot(alldata, "prediction",reduction.use="UMAP_on_CCA"),
dittoDimPlot(alldata, "bulk_prediction",reduction.use="UMAP_on_CCA"),
dittoDimPlot(alldata, "orig.ident",reduction.use="UMAP_on_CCA"),
dittoDimPlot(alldata, "seurat_clusters",reduction.use="UMAP_on_CCA"))

dev.off()

png("Alldata_markergenes_compartment.png", units="in", width=40, height=10, res=300)
cowplot::plot_grid(ncol = 6,
FeaturePlot(alldata,features = "KRT18" ,min.cutoff="q9", cols=c("lightgrey", "#02a7c0"), label=FALSE)+ labs(title = "Epithelial Cells"),
FeaturePlot(alldata,features = "CD3D" ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE)+ labs(title = "T Cells"),
FeaturePlot(alldata,features = "CD79A" ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE)+ labs(title = "B Cells"),
FeaturePlot(alldata,features = "LYZ" ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE)+ labs(title = "Myeloid"),
FeaturePlot(alldata,features = "COL1A1" ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE)+ labs(title = "Fibroblast"),
FeaturePlot(alldata,features = "CLDN5" ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE)+ labs(title = "Endothelial"))
dev.off()



png("Alldata_markergenes_compartment.png", units="in", width=40, height=10, res=300)
cowplot::plot_grid(ncol = 6,
FeaturePlot(alldata,features = "Krt18" ,min.cutoff="q9", cols=c("lightgrey", "#02a7c0"), label=FALSE)+ labs(title = "Epithelial Cells"),
FeaturePlot(alldata,features = "Lgals7" ,min.cutoff="q9", cols=c("lightgrey", "#02a7c0"), label=FALSE)+ labs(title = "Tumor Cells"),
FeaturePlot(alldata,features = "Cd3d" ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE)+ labs(title = "T Cells"),
FeaturePlot(alldata,features = "Cd79a" ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE)+ labs(title = "B Cells"),
FeaturePlot(alldata,features = "Col1a1" ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE)+ labs(title = "Fibroblast"),
FeaturePlot(alldata,features = "Cldn5" ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE)+ labs(title = "Endothelial"))
dev.off()


##########

Filtering using dplyr 

##########

#Filter p value < 0.05
#Filter  logFC between -2/+2
# Remove any NA
# mark UP or Down expression

library(dplyr)

xlsx::write.xlsx(c1c2 %>% filter(P.Value < 0.05) %>% filter(logFC < -2 | logFC 2) %>% filter(Symbol != "#N/A") %>% mutate(FOLD_UP_DOWN = ifelse((logFC) 2, "up", "down")) ,file="Diff_exp_MicroArray_All.xlsx", sheetName="c1c2", row.names=FALSE,append=TRUE)


##########

Deleting rows with Duplicated column based on condition of the other



##########

Ref: https://stackoverflow.com/questions/24011246/deleting-rows-that-are-duplicated-in-one-column-based-on-the-conditions-of-anoth

c5 = c5[order(c5[,'Symbol'],-c5[,'logFC']),]
c5= c5[!duplicated(c5$Symbol),]

, B-cells, Myeloid cells, Epithelial cells, Endothelial cells and Fibroblasts were identified.

/Users/akhaliq/Desktop/Bulk_rna/geo_data/GSE33532_RAW/GSM835269_02_B_TK.CEL.gz

#######
#https://bioinformatics.stackexchange.com/questions/5281/how-to-deal-with-duplicate-genes-having-different-expression-values
How to deal with duplicate genes having different expression values?


df2 <- aggregate(. ~ gene_name, data = df, max)

zeynebkurt@gmail.com, cihaterdogan@gmail.com, sultansevgiturgut@gmail.com, elcinguveyi@gmail.com
#####

|Cell Type   |Counts|
|:-----------|-----:|
|B-cells     |  9219|
|CAF         |   819|
|Endothelial |   352|
|Epithelial  |  8965|
|Myeloid     |   816|
|T-cells     | 22525|

|Cluster|Counts|
|:------|----: |
|CAF-S1 |  573 |
|CAF-S4 |    0 |
|Normal |   37 |

png("barchart_all_clusters.png", units="in", width=25, height=6, res=300)
cowplot::plot_grid(ncol = 3,
dittoBarPlot(all_merged, "seurat_clusters", group.by = "orig.ident"),
dittoBarPlot(all_merged, "seurat_clusters", group.by = "prediction"),
dittoBarPlot(all_merged, "seurat_clusters", group.by = "bulk_prediction"))
dev.off()



png("barchart_all_clusters2.png", units="in", width=25, height=6, res=300)
cowplot::plot_grid(ncol = 3,
dittoBarPlot(all_merged, "seurat_clusters", group.by = "orig.ident"),
dittoBarPlot(all_merged, "prediction", group.by = "orig.ident"),
dittoBarPlot(all_merged, "bulk_prediction", group.by = "orig.ident"))
dev.off()


#######

DefaultAssay(kefeer1) <- "RNA"
DefaultAssay(fibro) <- "RNA"

ifnb.list <- lapply(X = c(fibro_cafs1,kefeer1), FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features1 <- SelectIntegrationFeatures(object.list = ifnb.list)
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
# this command creates an 'integrated' data assay
immune.combined <- IntegrateData(anchorset = immune.anchors)
# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunTSNE(object = alldata, dims.use = 1:30)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined_0.5 <- FindClusters(immune.combined, resolution = 0.5)

#modify the metadata include CRC and BC
kefeer_meta1 <- read.csv("Meta_data_combined.csv",row.names=1)
immune.combined_0.5<- AddMetaData(immune.combined_0.5, kefeer_meta1)

DimPlot(immune.combined_0.5, reduction = "umap", split.by = "disease_type")

# For performing differential expression after integration, we switch back to the original data
DefaultAssay(immune.combined) <- "RNA"
nk.markers <- FindConservedMarkers(immune.combined_0.5, ident.1 = 6, grouping.var = "disease_type", verbose = FALSE)
head(nk.markers)
dim(nk.markers)

png("combined_BC_CRC_5.png", units="in", width=10, height=10, res=400)
FeaturePlot(immune.combined_0.5, features = c("GJB2", "SCARA5", "ADH1B", "SEMA3C", "CST1","TGFB1"), split.by = "disease_type", max.cutoff = 3, cols = c("grey", "red"),label.size = 4)
dev.off()



cell_typeA_marker_gene_list <- list(c("ASPN", "COL3A1", "THY1", "SFRP2", "COL10A1", "COL6A3", "LRRC17", "CILP", "GRP", "ITGBL1", "COL8A1", "COL14A1", "ADAM12", "OLFML2B", "ELN", "PLPP4", "CREB3L1", "FBN1", "LOXL1", "MATN3", "LRRC15", "COMP", "ISLR", "P3H1", "COL11A1", "SEPT11", "NBL1", "SPON1", "SULF1", "FNDC1", "CNN1", "MIAT", "MMP23B", "CPXM1", "FIBIN", "P4HA3", "GXYLT2", "CILP2", "P3H4", "CCDC80"))
object <- AddModuleScore(object = kefeer, features = cell_typeA_marker_gene_list, name = "ecm_myCAF")
FeaturePlot(object = object, features = "ecm_myCAF1", max.cutoff = 6, 
    cols = c("grey", "red"),label.size = 4)

library(patchwork)

plots <- VlnPlot(immune.combined_0.5, features = c("GJB2", "SCARA5", "ADH1B", "SEMA3C", "CST1","TGFB1"), split.by = "disease_type", pt.size = 0, combine = FALSE)

png("combined_BC_yuan_CRC_violin_plot.png", units="in", width=15, height=20, res=400)
wrap_plots(plots = plots, ncol = 1)
dev.off()



Idents(immune.combined_0.5) <- factor(Idents(immune.combined_0.5))
markers.to.plot <- c("GJB2", "SCARA5", "ADH1B", "SEMA3C", "CST1","TGFB1")
DotPlot(immune.combined_0.5, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "disease_type") + RotatedAxis()

dittoDimPlot(immune.combined_0.5, markers.to.plot, split.by = "disease_type")



png("combined_BC_CRC_4.png", units="in", width=20, height=20, res=400)
cowplot::plot_grid(ncol = 1,
FeaturePlot(immune.combined_0.5, features = "GJB2", split.by = "disease_type", max.cutoff = 3, 
    cols = c("grey", "#cc0000"),pt.size = 2,label.size = 4),
FeaturePlot(immune.combined_0.5, features = "SCARA5", split.by = "disease_type", max.cutoff = 3, 
    cols = c("grey", "#009933")),
FeaturePlot(immune.combined_0.5, features = "ADH1B", split.by = "disease_type", max.cutoff = 3, 
    cols = c("grey", "#ff9966")),
FeaturePlot(immune.combined_0.5, features = "SEMA3C", split.by = "disease_type", max.cutoff = 3, 
    cols = c("grey", "#0000ff")),
FeaturePlot(immune.combined_0.5, features = "CST1", split.by = "disease_type", max.cutoff = 3, 
    cols = c("grey", "#cc9900")),
FeaturePlot(immune.combined_0.5, features ="TGFB1", split.by = "disease_type", max.cutoff = 3, 
    cols = c("grey", "#055755")))
dev.off()


png("integrated_analysis1.png", units="in", width=6, height=4, res=400)
cowplot::plot_grid(ncol = 2,
DimPlot(immune.combined_0.5, reduction = "umap", group.by = "disease_type"),
DimPlot(immune.combined_0.5, reduction = "umap", label = FALSE, repel = TRUE))
dev.off()


png("integrated_analysis1.png", units="in", width=6, height=6, res=400)
p1+p2
dev.off()

######

############# CONICSmat Copy number variation analysis

Ref : https://github.com/diazlab/CONICS/wiki/Tutorial---CONICSmat;---Dataset:-SmartSeq2-scRNA-seq-of-Oligodendroglioma

library(CONICSmat)

# The expression is provided as log2(CPM/10+1)

epi_both<- NormalizeData(epi,scale.factor = 1e6)
epi_both_log_data <- GetAssayData(epi_both)
epi_log2_data <- log(expm1(epi_both_log_data) + 1, 2)

suva_expr = as.matrix(epi_log2_data)

suva_expr [which(is.na(suva_expr ))]=0

suva_expr[1:5,1:5]

gene_pos=getGenePositions(rownames(suva_expr))
suva_expr=filterMatrix(suva_expr,gene_pos[,"hgnc_symbol"],minCells=5)
normFactor=calcNormFactors(suva_expr)
l=plotAll(suva_expr,normFactor,regions,gene_pos,"SUVA_CNVs_ALL")


png("patients_heatmap.png", units="in", width=30, height=20, res=400)
hi=plotHistogram(l,suva_expr,clusters=2,zscoreThreshold=4,patients)
dev.off()

png("patients_heatmap2.png", units="in", width=30, height=20, res=400)
hi=plotHistogram(l[,candRegions],suva_expr,clusters=4,zscoreThreshold=4,patients)
dev.off()


### CASPER


epi.subset <- subset(epi.subset, subset = nFeature_RNA 200 & nFeature_RNA < 2500 & percent.mt < 5)
epi.subset <- NormalizeData(epi.subset , scale.factor = 1e6, normalization.method = "RC")
epi.subset <- FindVariableFeatures(epi.subset, do.plot = T, nfeatures = 1000)
epi.subset <- ScaleData(epi.subset)

epi.subset <- RunPCA(epi.subset, features = VariableFeatures(object = epi.subset),npcs = 100)
epi.subset <- RunTSNE(epi.subset, dims.use = 1:10)
DimPlot(epi.subset, reduction = "tsne")
FeaturePlot(epi.subset, features = c("SDC1", "CD38"))

epi.subset <- FindNeighbors(epi.subset, dims = 1:10)
epi.subset <- FindClusters(epi.subset, resolution = 0.5)
DimPlot(epi.subset, reduction = "tsne", label=T)

log.ge <- as.matrix(epi.subset@assays$RNA@data)
control <- names(Idents(epi.subset) )[Idents(epi.subset) %in% c(2,7)]


genes <- rownames(log.ge)
annotation <- generateAnnotation(id_type="hgnc_symbol", genes=genes, centromere=centromere, ishg19 = T)
log.ge <- log.ge[match( annotation$Gene,rownames(log.ge)) , ]
rownames(log.ge) <- annotation$Gene
log.ge <- log2(log.ge +1)


load("maf.rda") ## from https://github.com/akdess/CaSpER/blob/master/data/maf.rda
loh<- list()
loh[[1]] <- maf
names(loh) <- "epi.subset"

loh.name.mapping <- data.frame (loh.name= "epi.subset" , sample.name=colnames(log.ge))

object <- CreateCasperObject(raw.data=log.ge,loh.name.mapping=loh.name.mapping, sequencing.type="single-cell", 
cnv.scale=3, loh.scale=3, 
expr.cutoff=0.1, filter="median", matrix.type="normalized",
annotation=annotation, method="iterative", loh=loh, 
control.sample.ids=control, cytoband=cytoband)



pdf("epi.subset_tumor.Distrubution.pdf")
plot(density(as.vector(object@control.normalized[[3]])))
plot(density(log2(object@control.normalized.noiseRemoved[[3]]+1)))
dev.off()


## runCaSpER
final.objects <- runCaSpER(object, removeCentromere=T, cytoband=cytoband, method="iterative")

## summarize large scale events 
finalChrMat <- extractLargeScaleEvents (final.objects, thr=0.75) 

obj <- final.objects[[9]]
plotHeatmap10x(object=obj, fileName="heatmap.png",cnv.scale= 3, cluster_cols = F, cluster_rows = T, show_rownames = T, only_soi = T)

#### VISUALIZATION 
chrMat <- finalChrMat
plot.data <- melt(chrMat)
plot.data$value2 <- "neutral"
plot.data$value2[plot.data$value 0] <- "amplification"
plot.data$value2[plot.data$value < 0] <- "deletion"
plot.data$value2 <- factor(plot.data$value2, levels = c("amplification", 
    "deletion", "neutral"))
plot.data$X2 <- factor(plot.data$X2, levels = colnames(chrMat))
p <- ggplot(aes(x = X2, y = X1, fill = value2), data = plot.data) + 
    geom_tile(colour = "white", size = 0.01) + 
    labs(x = "", 
    y = "") + scale_fill_manual(values = c(amplification = muted("red"), 
    deletion = muted("blue"), neutral = "white")) + theme_grey(base_size = 6) + 
    theme(legend.position = "right", legend.direction = "vertical", 
        legend.title = element_blank(), strip.text.x = element_blank(), 
        legend.text = element_text(colour = "black", size = 7, 
            face = "bold"), legend.key.height = grid::unit(0.8, 
            "cm"), legend.key.width = grid::unit(0.5, "cm"), 
        axis.text.x = element_text(size = 5, colour = "black", 
            angle = -45, hjust = 0), axis.text.y = element_text(size = 6, 
            vjust = 0.2, colour = "black"), axis.ticks = element_line(size = 0.4), 
        plot.title = element_text(colour = "black", hjust = 0, 
            size = 6, face = "bold"))



######### running Docker images

#### To give permission for a new docker image
sudo chmod 666 /var/run/docker.sock
Ref : https://www.digitalocean.com/community/questions/how-to-fix-docker-got-permission-denied-while-trying-to-connect-to-the-docker-daemon-socket
#####

docker images

It will show
REPOSITORY             TAG       IMAGE ID       CREATED        SIZE
satijalab/seurat       latest    6051045e47e1   3 weeks ago    3.62GB
hello-world            latest    d1165f221234   5 weeks ago    13.3kB
hello-world            latest    d1165f221234   5 weeks ago    13.3kB
trinityctat/infercnv   latest    2913b95aecf9   2 months ago   4.57GB

docker run -i -t 6051045e47e1

############################################################################################################################################

SingleR : DVRN's Pipeline

#############################################################################################################################################

#Create singleR object
#Add the Existing seurat object into that
#plot the umap
#then add the singler annotations to the seurat metadata.


singler.epi_tumor.cells = CreateSinglerObject(counts=as.matrix(epi.subset@assays$RNA@counts), project.name = "Epithelial Tumor Annotations", annot = "Epithelial_metadata.txt" , min.genes = 200,
  technology = "10X", species = "Human", citation = "",
  ref.list = list(), normalize.gene.length = F, variable.genes = "de",
  fine.tune = F, do.signatures = T, do.main.types = T, 
  reduce.file.size = FALSE, numCores = 16)

table_singler_kefwfer <- table(singler.epi_tumor.cells$singler[[2]]$SingleR.single$labels,epi.subset@meta.data$orig.ident))

singler.epi_tumor.cells$seurat <- epi.subset
singler.epi_tumor.cells$meta.data$orig.ident = epi.subset@meta.data$orig.ident
singler.epi_tumor.cells$meta.data$xy = epi.subset@reductions$umap@cell.embeddings
singler.epi_tumor.cells$meta.data$clusters = epi.subset@active.ident



singler.korean.epi = CreateSinglerObject(counts=as.matrix(kor.epi@assays$RNA@counts), project.name = "korean data", annot = data.frame(kor.epi@meta.data) , min.genes = 200, technology = "10X", species = "Human", citation = "",ref.list = list(), normalize.gene.length = F, variable.genes = "de",fine.tune = F, do.signatures = T, do.main.types = T, reduce.file.size = FALSE, numCores = 16)


out_epi_all_ident = SingleR.PlotTsne(singler.epi_tumor.cells$singler[[1]]$SingleR.single,
singler.epi_tumor.cells$meta.data$xy, do.label = FALSE, do.letters = F,
labels=singler.epi_tumor.cells$meta.data$orig.ident, label.size = 15, dot.size = 3)

out_epi_all = SingleR.PlotTsne(singler.epi_tumor.cells$singler[[1]]$SingleR.single,
singler.epi_tumor.cells$meta.data$xy, do.label = FALSE, do.letters = F,
labels=singler.epi_tumor.cells$singler[[1]]$SingleR.single$labels, label.size = 20, dot.size = 4)



png("SingleR_epi_all_UMAP.png", units="in", width=20, height=20, res=300)
out_epi_all$p + xlab('UMAP1') + ylab('UMAP2')
dev.off()

png("SingleR_OUR_all_UMAP.png", units="in", width=20, height=20, res=300)
out_epi_all_ident$p + xlab('UMAP1') + ylab('UMAP2')
dev.off()

## Exporting singler annotations to the csv and then adding it to the seurat object metadata

# for database [1] HPCA
write.csv(singler.epi_tumor.cells$singler[[1]]$SingleR.single.main$labels, "singler_1_types.csv") 

# for database [2] Blueprint-Encode
write.csv(singler.epi_tumor.cells$singler[[2]]$SingleR.single.main$labels, "singler_2_types.csv") 


## adding Metadata

singler1 <- read.csv("/Users/akhaliq/Desktop/infercnv/singler_1.csv",header=TRUE,sep=',', row.names=1)
epi.subset1 <- AddMetaData(epi.subset, singler1)

singler2 <- read.csv("/Users/akhaliq/Desktop/infercnv/singler_2.csv",header=TRUE,sep=',', row.names=1)
epi.subset1 <- AddMetaData(epi.subset1, singler2)


epi_sub_encode <- subset(x = epi, subset = singler_encode == "Epithelial cells")

####################################################

Converting Raw Counts Gene Name to Ensembl IDs
# Run scCancer

####################################################

biocLite("biomaRt")
library(biomaRt)

library("org.Hs.eg.db")

#BiocManager::install("EnsDb.Hsapiens.v79")

library(EnsDb.Hsapiens.v79)


gene_ids = rownames(epi.counts)
head(gene_ids)

epi.counts <-as.matrix(epi.subset.only@assays$RNA@counts)

geneIDs2 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= gene_ids, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
head(geneIDs2)
dim(geneIDs2)
dim(epi.counts)

library(tidyverse)

# geneIDs2 %>% remove_rownames %>% column_to_rownames(var="SYMBOL")
#epi.counts1<- tibble::rownames_to_column(epi.counts, "SYMBOL")
ls()

epi.counts1<- tibble::rownames_to_column(as.data.frame(epi.counts), "SYMBOL")
head(epi.counts1[,1:3])

merged <- merge(epi.counts1,geneIDs2,by="SYMBOL")
head(merged[,1:3])

# remove duplicate

merged1 <- merged[!duplicated(merged$SYMBOL),] 
dim(merged1)

head(merged1[,1:3])


# 1.Select two coloumns to make the gene.info file
 
library(dplyr)
gene.info <- merged1 %>% select (SYMBOL,GENEID)

# 2. flip two coloms
gene.info <- gene.info[c("GENEID","SYMBOL")]

#3. remove header
names(gene.info) <- NULL

head(gene.info)

# Remove the last Ensemble ID coloumn in the merged1 Data Matrix

merged1$GENEID <- NULL


# make the first coloumn as rowname in merged1

library(tidyverse)

merged1 <- merged1 %>% remove_rownames %>% column_to_rownames(var="SYMBOL")
head(merged1[,1:5])

# save rds 

saveRDS(merged1,"counts.rds")
saveRDS(gene.info,"gene.info.rds")

# then read it in the server where u want to run scCancer

counts <- readRDS("counts.rds")
gene.info <- readRDS("gene.info.rds")

# Generate a 10X-like data folder based on the data matrix and gene information, which can be used directly to perform scCancer analysis.

library(scCancer)

generate10Xdata(counts, gene.info, outPath="/data/sccancer/input", overwrite = F)

	
library(scCancer)

https://github.com/velocyto-team/velocyto.R/issues/57#issuecomment-506671728

dataPath <- "/data/sctyper/raiser/input"     # The path of cell ranger processed data
savePath <- "//data/sctyper/raiser/results"  # A path to save the results
sampleName <- "Reiser Sc Dataset"          # The sample name
authorName <- "Ateeq Khaliq"           # The author name to mark the report
statPath <- "/data/sctyper/raiser/results" 


stat.results <- runScStatistics(
    dataPath = dataPath,
    savePath = savePath,
    sampleName = sampleName,
    authorName = authorName
)

anno.results <- runScAnnotation(
    dataPath = dataPath,
    statPath = statPath,


    savePath = savePath,
    authorName = authorName,
    sampleName = sampleName,
    geneSet.method = "average"   # or "GSVA"
)

###


################################################################################################################################

Running HCL - Human cell Landscape to annotate the cells in the scExpression Matrix

################################################################################################################################

library(scHCL)

epi_hcl_tumor.only <- scHCL(scdata = epi.subset.only@assays$RNA@counts, numbers_plot = 3)

scHCL_vis(epi_hcl)
scHCL_vis(epi_hcl_tumor.only)
#open the browser in the plotly section save the image and csv file for further analysis

head(as.matrix(epi_hcl_tumor.only$scHCL))
write.csv(epi_hcl_tumor.only$scHCL,"relevant.cell.type_HCL.csv") # name the column 
hcl_anno <- read.csv("/Users/akhaliq/Desktop/hcl/only_tumor/relevant.cell.type_HCL.csv", header=TRUE, sep = ',',row.names=1)

# add metadata to the seurat object

epi.subset.only <- AddMetaData(epi.subset.only, hcl_anno)

library(dittoSeq)

png("Epi_only_tumor_samples_HCL_anno.png", units="in", width=50, height=20, res=400)
dittoBarPlot(epi.subset.only, "orig.ident", group.by = "HCL_anno")
dev.off()

/Users/akhaliq/Desktop/merged_epi/gsva_graphs/GSVA_Signature.svg

png("Epi_only_tumor_samples_singleR_HPCA_anno.png", units="in", width=10, height=10, res=400)
dittoBarPlot(epi.subset.only, "orig.ident", group.by = "SingleR.HPCA")
dev.off()

png("Epi_only_tumor_samples_Singler_Blueprint_Encode_anno.png", units="in", width=10, height=10, res=400)
dittoBarPlot(epi.subset.only, "orig.ident", group.by = "Singler.Blueprint.Encode")
dev.off()

png("Epi_only_tumor_samples_HCL_anno_dim_plot.png", units="in", width=20, height=10, res=400)
dittoDimPlot(epi.subset.only, "HCL_anno")
dev.off()


png("Epi_only_tumor_samples_SingleR.HPCA_dim_plot.png", units="in", width=20, height=10, res=400)
dittoDimPlot(epi.subset.only, "SingleR.HPCA")
dev.off()


png("Epi_only_tumor_samples_Singler.Blueprint.Encode.png", units="in", width=10, height=10, res=400)
dittoDimPlot(epi.subset.only, "Singler.Blueprint.Encode")
dev.off()


########
Create directory in R
########

 dir.create("/Users/akhaliq/Desktop/epithelial_data")


#########
Adding user to server 
########

sudo su
adduser cihat
passwd cihat
usermod -aG wheel cihat # giving Root permission

# logging in 

ssh masoodlab2@157.55.253.201 --pssd "sunrise_1712" # root
ssh cihat@157.55.253.201
ssh zeyneb@157.55.253.201
ssh sevgi@157.55.253.201
ssh elcin@157.55.253.201

## to delete a usr
userdel zeynab
######


####

Cellphonedb

####

conda activate cpdb
cellphonedb method statistical_analysis alldata_meta.txt log_counts.txt --counts-data gene_name --threads 8 
#[ ][CORE][06/06/21-03:13:45][INFO] [Cluster Statistical Analysis] Threshold:0.1 Iterations:1000 Debug-seed:-1 Threads:8 Precision:3
cellphonedb plot dot_plot --rows rows.txt --columns columns.txt



#####
ccfindR cNMF
#####


library(ccfindR)

sc <- scNMFSet(count = epi_final@assays$RNA@counts)
sc <- filter_cells(sc, umi.min = 10^2.6, umi.max = 10^3.4)

	
######
Cancer Subtypes
######

library(Seurat)
library("CancerSubtypes")
epi_tumor<- subset(x = epi_final, subset = Condition == "Tumor")

epi_log2= log2(exp(as.matrix(GetAssayData(object = epi_tumor, slot = "data"))))
data.checkDistribution(epi_log2)
epi_log2_fsvar=FSbyVar(epi_log2, cut.type = "topk",4000)
epi_log2_fsvar_norm=data.normalization(epi_log2_fsvar)
epi_cc=ExecuteCC(clusterNum=3,d=epi_log2_fsvar_norm, maxK=5,clusterAlg="hc", distance="pearson",title="Epi Tumor")
save.image("epi_cluster.RData")
epi_cnmf =ExecuteCNMF(epi_log2_fsvar_norm,clusterNum=3,nrun=30)
save.image("epi_cluster.RData")

subtypes <- as.data.frame(epi_cnmf2$group)
colnames(subtypes) <- "type"
subtypes$type <- paste0('subtype_', subtypes$type)
table(subtypes)
subtypes
subtype_1 subtype_2 
     2806      2982 
write.csv(as.data.frame(subtypes,colnames(epi_log2_fsvar_norm)),"nmf_N2.csv")

# Changing colnames 

colnames(subtypes)[1] <- "NMF_N2"

# add metadata
epi_tumor<- AddMetaData(epi_tumor, subtypes)


epi_cnmf3 =ExecuteCNMF(epi_log2_fsvar_norm,clusterNum=3,nrun=30)
gc()
epi_cnmf4 =ExecuteCNMF(epi_log2_fsvar_norm,clusterNum=4,nrun=30)
gc()
save.image("epi_cluster.RData")

epi_cnmf5 =ExecuteCNMF(epi_log2_fsvar_norm,clusterNum=5,nrun=30)
gc()
save.image("epi_cluster.RData")

epi_cnmf6 =ExecuteCNMF(epi_log2_fsvar_norm,clusterNum=6,nrun=30)
gc()
save.image("epi_cluster.RData")

epi_cnmf7 =ExecuteCNMF(epi_log2_fsvar_norm,clusterNum=7,nrun=30)
gc()
save.image("epi_cluster.RData")

epi_cnmf8 =ExecuteCNMF(epi_log2_fsvar_norm,clusterNum=8,nrun=30)
gc()
save.image("epi_cluster.RData")

epi_cnmf9 =ExecuteCNMF(epi_log2_fsvar_norm,clusterNum=9,nrun=30)
gc()
save.image("epi_cluster.RData")

epi_cnmf10 =ExecuteCNMF(epi_log2_fsvar_norm,clusterNum=10,nrun=30)
gc()
save.image("epi_cluster.RData")


plot(sil(epi_cnmf2$originalResult), col=c("maroon","pink"), border=NA)

png("Sil_nmf.png", units="in", width=10, height=20, res=300)
cowplot::plot_grid(ncol = 3,
plot(sil2),
plot(sil3),
plot(sil4),
plot(sil5),
plot(sil6)
)

png("epi_tumor_nmf.png", units="in", width=30, height=20, res=300)
cowplot::plot_grid(ncol = 4,
dittoDimPlot(epi_tumor, "seurat_clusters",do.label = TRUE, labels.repel = FALSE,reduction="umap"),
dittoDimPlot(epi_tumor, "Condition",do.label = TRUE, labels.repel = FALSE,reduction="umap"),
dittoDimPlot(epi_tumor, "Location",do.label = TRUE, labels.repel = FALSE,reduction="umap"),
dittoDimPlot(epi_tumor, "MSI_Status",do.label = TRUE, labels.repel = FALSE,reduction="umap"),
dittoDimPlot(epi_tumor, "prediction",do.label = TRUE, labels.repel = FALSE,reduction="umap"),
dittoDimPlot(epi_tumor, "bulk_prediction",do.label = TRUE, labels.repel = FALSE,reduction="umap"),

dittoDimPlot(epi_tumor, "NMF_N2",do.label = TRUE, labels.repel = FALSE,reduction="umap")

)
dev.off()



png("epi_tumor_nmf_tsne.png", units="in", width=20, height=10, res=300)
cowplot::plot_grid(ncol = 4,
dittoDimPlot(epi_tumor, "seurat_clusters",do.label = TRUE, labels.repel = FALSE,reduction="tsne"),
dittoDimPlot(epi_tumor, "Condition",do.label = TRUE, labels.repel = FALSE,reduction="tsne"),
dittoDimPlot(epi_tumor, "Location",do.label = TRUE, labels.repel = FALSE,reduction="tsne"),
dittoDimPlot(epi_tumor, "MSI_Status",do.label = TRUE, labels.repel = FALSE,reduction="tsne"),
dittoDimPlot(epi_tumor, "prediction",do.label = TRUE, labels.repel = FALSE,reduction="tsne"),
dittoDimPlot(epi_tumor, "bulk_prediction",do.label = TRUE, labels.repel = FALSE,reduction="tsne"),

dittoDimPlot(epi_tumor, "NMF_N2",do.label = TRUE, labels.repel = FALSE,reduction="tsne")

)
dev.off()



png("epi_tumor_nmf_sil2.png", units="in", width=10, height=10, res=300)
plot(silhouette(epi_cnmf2$originalResult), col=c("maroon","pink"), border=NA)
dev.off()

png("epi_tumor_nmf_sil3.png", units="in", width=10, height=10, res=300)
plot(silhouette(epi_cnmf3$originalResult), col=c("maroon","pink"), border=NA)
dev.off()

png("epi_tumor_nmf_sil4.png", units="in", width=10, height=10, res=300)
plot(silhouette(epi_cnmf4$originalResult), col=c("maroon","pink"), border=NA)
dev.off()

png("epi_tumor_nmf_sil5.png", units="in", width=10, height=10, res=300)
plot(silhouette(epi_cnmf5$originalResult), col=c("maroon","pink"), border=NA)
dev.off()

png("epi_tumor_nmf_sil6.png", units="in", width=10, height=10, res=300)
plot(silhouette(epi_cnmf6$originalResult), col=c("maroon","pink"), border=NA)
dev.off()




The enriched ligand–receptor interactions between two cell states on the basis of expression of a receptor by one cell state and a ligand by another cell state. For
each gene in the cluster, the percentage of cells expressing the gene and the gene expression mean were calculated. We consider the expression levels of ligands and receptors within each cell state
and use empirical shuffling to calculate which ligand–receptor pairs display significant cell-state specificity. Specificity of the ligand–receptor interaction ligand–receptor pairs are ubiquitously expressed by the cells in a tissue and
therefore are not informative regarding specific communication between particular cell states.


png("epi_tumor_gsva.png", units="in", width=10, height=10, res=300)
plot_gsva_pathway(gsva_result, pathway_id = rownames(max_difference)[1])
dev.off()

png("epi_tumor_gsva_heatmap.png", units="in", width=10, height=10, res=300)
plot_gsva_heatmap(gsva_result, max_pathways = 20, margins = c(6,20))
dev.off()

####################################################################
# Load all the Celltypes

load("/Users/akhaliq/Desktop/yuan_data/new_data/Recent_alldata_bk.031821.RData")

####################################################################

#June 2021

# Single-cell RNA-seq analysis - Pseudobulk DE analysis with DESeq2  

####################################################################

Ref: https://hbctraining.github.io/scRNA-seq/lessons/pseudobulk_DESeq2_scrnaseq.html

extract the cells from a Seurat object, which we had created at the end of the single-cell analysis workflow, we could use code similar to that below

library(SingleCellExperiment)

# Bring in Seurat object
seurat <- readRDS("path/to/seurat.rds")

# Extract raw counts and metadata to create SingleCellExperiment object
counts <- epi_tumor@assays$RNA@counts 

metadata <- epi_tumor@meta.data

# Set up metadata as desired for aggregation and DE analysis
metadata$cluster_id <- factor(epi_tumor@meta.data$bulk_prediction)
metadata$sample_id <- factor(epi_tumor@meta.data$orig.ident)

# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts), 
                           colData = metadata)

# Identify groups for aggregation of counts
groups <- colData(sce)[, c("cluster_id", "sample_id")]



dir.create("/data/pseudobulk/data")
dir.create("/data/pseudobulk/results")
dir.create("/data/pseudobulk/figures")

suppressPackageStartupMessages({
library(scater)
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)
#library(scPipe)
})



*** For every cell, we have information about the associated condition (ctrl or stim), sample ID, and cell type. We will use this information to perform the differential expression analysis between conditions for any particular cell type of interest.

# updated_myData <- subset(metadata, cluster_id!="NA")

/Users/akhaliq/Desktop/pseudobulk

clusters <- levels(metadata$cluster_id)
clusters
metadata <- read.csv("without_na.csv", header=T,row.names=1,sep=',')
head(metadata)
sce <- SingleCellExperiment(assays = list(counts = counts), 
                           colData = metadata)
sce
groups <- colData(sce)[, c("cluster_id", "sample_id")]
groups
assays(sce)
counts(sce)[1:6, 1:6]
dim(colData(sce))
head(colData(sce))
# Named vector of cluster names
kids <- purrr::set_names(levels(sce$cluster_id))
kids
# Total number of clusters
nk <- length(kids)
nk
# Named vector of sample names
sids <- purrr::set_names(levels(sce$sample_id))
# Total number of samples 
ns <- length(sids)
ns
## Determine the number of cells per sample
table(sce$sample_id)
## Turn named vector into a numeric vector of number of cells per sample
n_cells <- as.numeric(table(sce$sample_id))
## Determine how to reoder the samples (rows) of the metadata to match the order of sample names in sids vector
m <- match(sids, sce$sample_id)
## Create the sample level metadata by combining the reordered metadata with the number of cells corresponding to each sample.
ei <- data.frame(colData(sce)[m, ], 
                  n_cells, row.names = NULL) %>% 
                select(-"cluster_id")
ei
library(scater)
qc <- perCellQCMetrics(sce)
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce <- sce[, !ol]
dim(sce)
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
dim(sce)
sce <- computeLibraryFactors(sce)
sce <- logNormCounts(sce)
sce
# Subset metadata to only include the cluster and sample IDs to aggregate across
groups <- colData(sce)[, c("cluster_id", "sample_id")]
# Aggregate across cluster-sample groups
pb <- aggregate.Matrix(t(counts(sce)), 
                       groupings = groups, fun = "sum") 
class(pb)
dim(pb)
pb[1:6, 1:6]
splitf <- sapply(stringr::str_split(rownames(pb), 
                                    pattern = "_",  
                                    n = 2), 
                 `[`, 1)
pb <- split.data.frame(pb, 
                       factor(splitf)) %>%
        lapply(function(u) 
                set_colnames(t(u), 
                             stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))
class(pb)
# Explore the different components of list
str(pb)
options(width = 100)
table(sce$cluster_id, sce$sample_id)
get_sample_ids <- function(x){
        pb[[x]] %>%
                colnames()
}
de_samples <- map(1:length(kids), get_sample_ids) %>%
        unlist()
samples_list <- map(1:length(kids), get_sample_ids)
get_cluster_ids <- function(x){
        rep(names(pb)[x], 
            each = length(samples_list[[x]]))
}
de_cluster_ids <- map(1:length(kids), get_cluster_ids) %>%
        unlist()
# Create a data frame with the sample IDs, cluster IDs and condition
gg_df <- data.frame(cluster_id = de_cluster_ids,
                    sample_id = de_samples)
gg_df <- left_join(gg_df, ei[, c("sample_id", "group_id")]) 
metadata <- gg_df %>%
        dplyr::select(cluster_id, sample_id, group_id) 
      
        
metadata        
clusters <- levels(metadata$cluster_id)
clusters
clusters[1]
cluster_metadata <- metadata[which(metadata$cluster_id == clusters[1]), ]
head(cluster_metadata)
cluster_metadata
rownames(cluster_metadata) <- cluster_metadata$sample_id
head(cluster_metadata)
# Subset the counts to only the B cells
counts <- pb[[clusters[1]]]
cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])
# Check that all of the row names of the metadata are the same and in the same order as the column names of the counts in order to use as input to DESeq2
all(rownames(cluster_metadata) == colnames(cluster_counts))   
dds <- DESeqDataSetFromMatrix(cluster_counts, 
                              colData = cluster_metadata, 
                              design = ~ group_id)
dds
rld <- rlog(dds, blind=TRUE)
DESeq2::plotPCA(rld, intgroup = "group_id")
png("pca_epi_cms.png", units="in", width=10, height=10, res=300)
DESeq2::plotPCA(rld, intgroup = "sample_id")
dev.off()
rld_mat <- assay(rld)
head(rld_mat)
write.csv(assay(rld),"logcoutnts_pseudobulk.csv")
pheatmap(rld_cor, annotation = cluster_metadata[, c("group_id"), drop=F])
rld_cor <- cor(rld_mat)
pheatmap(rld_cor, annotation = cluster_metadata[, c("group_id"), drop=F])
png("pairwaise_epi_cms.png", units="in", width=10, height=10, res=300)
pheatmap(rld_cor, annotation = cluster_metadata[, c("group_id"), drop=F])
dev.off()
dds <- DESeq(dds)
plotDispEsts(dds)
plotDispEsts(dds)
png("Dispersion_estimates_epi_cms.png", units="in", width=10, height=10, res=300)
plotDispEsts(dds)
dev.off()
levels(cluster_metadata$group_id)
levels(cluster_metadata$group_id)[2]
levels(cluster_metadata$group_id)[1]
levels(cluster_metadata$group_id)[3]
levels(cluster_metadata$group_id)[4]


group_id
ls()
contrast.1_2 <- c("group_id", levels(cluster_metadata$group_id)[1], levels(cluster_metadata$group_id)[2])
cms1-cms2
contrast.1_2 <- c("group_id", levels(cluster_metadata$group_id)[1], levels(cluster_metadata$group_id)[2])
res_1_2 <- results(dds,contrast = contrast.1_2, alpha = 0.05)
res_1_2 <- lfcShrink(dds, contrast =  contrast.1_2,res=res_1_2)
#cms1-cms3
contrast.1_2 <- c("group_id", levels(cluster_metadata$group_id)[1], levels(cluster_metadata$group_id)[3])
res_1_3 <- results(dds,contrast = contrast.1_3, alpha = 0.05)
res_1_3 <- lfcShrink(dds, contrast =  contrast.1_3,res=res_1_3)
#cms1-cms4
contrast.1_4 <- c("group_id", levels(cluster_metadata$group_id)[1], levels(cluster_metadata$group_id)[4])
res_1_4 <- results(dds,contrast = contrast.1_4, alpha = 0.05)
res_1_4 <- lfcShrink(dds, contrast =  contrast.1_4,res=res_1_4)
#cms2-cms3
contrast.2_3 <- c("group_id", levels(cluster_metadata$group_id)[2], levels(cluster_metadata$group_id)[3])
res_2_3 <- results(dds,contrast = contrast.2_3, alpha = 0.05)
res_2_3 <- lfcShrink(dds, contrast =  contrast.2_3,res=res_2_3)
#cms2-cms4
contrast.2_4 <- c("group_id", levels(cluster_metadata$group_id)[2], levels(cluster_metadata$group_id)[4])
res_2_4 <- results(dds,contrast = contrast.2_4, alpha = 0.05)
res_2_4 <- lfcShrink(dds, contrast =  contrast.2_4,res=res_2_4)
#cms3-cms4
contrast.3_4 <- c("group_id", levels(cluster_metadata$group_id)[3], levels(cluster_metadata$group_id)[4])
res_3_4 <- results(dds,contrast = contrast.3_4, alpha = 0.05)
res_3_4 <- lfcShrink(dds, contrast =  contrast.3_4,res=res_3_4)
contrast.1_2 <- c("group_id", levels(cluster_metadata$group_id)[1], levels(cluster_metadata$group_id)[2])
res_1_2 <- results(dds,contrast = contrast.1_2, alpha = 0.05)
res_1_2 <- lfcShrink(dds, contrast =  contrast.1_2,res=res_1_2)
contrast.1_3 <- c("group_id", levels(cluster_metadata$group_id)[1], levels(cluster_metadata$group_id)[3])
res_1_3 <- results(dds,contrast = contrast.1_3, alpha = 0.05)
res_1_3 <- lfcShrink(dds, contrast =  contrast.1_3,res=res_1_3)
contrast.1_3 <- c("group_id", levels(cluster_metadata$group_id)[1], levels(cluster_metadata$group_id)[3])
res_1_3 <- results(dds,contrast = contrast.1_3, alpha = 0.05)
res_1_3 <- lfcShrink(dds, contrast =  contrast.1_3,res=res_1_3)
contrast.1_3 <- c("group_id", levels(cluster_metadata$group_id)[1], levels(cluster_metadata$group_id)[3])
res_1_3 <- results(dds,contrast = contrast.1_3, alpha = 0.05)
res_1_3 <- lfcShrink(dds, contrast =  contrast.1_3,res=res_1_3)
#cms2-cms4
contrast.2_4 <- c("group_id", levels(cluster_metadata$group_id)[2], levels(cluster_metadata$group_id)[4])
res_2_4 <- results(dds,contrast = contrast.2_4, alpha = 0.05)
res_2_4 <- lfcShrink(dds, contrast =  contrast.2_4,res=res_2_4)
#cms3-cms4
contrast.3_4 <- c("group_id", levels(cluster_metadata$group_id)[3], levels(cluster_metadata$group_id)[4])
res_3_4 <- results(dds,contrast = contrast.3_4, alpha = 0.05)
res_3_4 <- lfcShrink(dds, contrast =  contrast.3_4,res=res_3_4)
# Turn the results object into a tibble for use with tidyverse functions
res_tbl_1_2 <- res_1_2 %>%
        data.frame() %>%
        rownames_to_column(var="gene") %>%
        as_tibble()
res_tbl_1_2
# Write all results to file
write.csv(res_tbl_1_2,
          paste0("results/", clusters[1], "_", levels(cluster_metadata$group_id)[1], "_vs_", levels(cluster_metadata$group_id)[2], "_all_genes.csv"),
          quote = FALSE, 
          row.names = FALSE)
# Write all results to file
write.csv(res_tbl_1_2,
          paste0("results/", clusters[1], "_", "_all_genes.csv"),
          quote = FALSE, 
          row.names = FALSE)
# Write all results to file
write.csv(res_tbl_1_2,
          paste0("results/", clusters[1], "_", levels(cluster_metadata$group_id)[1], "_vs_", levels(cluster_metadata$group_id)[2], "_all_genes.csv"),
          quote = FALSE, 
          row.names = FALSE)
#cms1-cms3
# Turn the results object into a tibble for use with tidyverse functions
res_tbl_1_3 <- res_1_3  %>%
        data.frame() %>%
        rownames_to_column(var="gene") %>%
        as_tibble()
# Check results output
res_tbl_1_2
# Write all results to file
write.csv(res_tbl_1_3,
          paste0("results/", clusters[1], "_", levels(cluster_metadata$group_id)[1], "_vs_", levels(cluster_metadata$group_id)[3], "_all_genes.csv"),
          quote = FALSE, 
          row.names = FALSE)
#cms1-cms4
 # Turn the results object into a tibble for use with tidyverse functions
res_tbl_1_4 <- res_1_4 %>%
        data.frame() %>%
        rownames_to_column(var="gene") %>%
        as_tibble()
# Check results output
res_tbl_1_4
# Write all results to file
write.csv(res_tbl_1_4,
          paste0("results/", clusters[1], "_", levels(cluster_metadata$group_id)[1], "_vs_", levels(cluster_metadata$group_id)[4], "_all_genes.csv"),
          quote = FALSE, 
          row.names = FALSE)
#cms2-cms3
# Turn the results object into a tibble for use with tidyverse functions
res_tbl_2_3 <- res_2_3 %>%
        data.frame() %>%
        rownames_to_column(var="gene") %>%
        as_tibble()
# Check results output
res_tbl_2_3
# Write all results to file
write.csv(res_tbl_2_3,
          paste0("results/", clusters[1], "_", levels(cluster_metadata$group_id)[2], "_vs_", levels(cluster_metadata$group_id)[3], "_all_genes.csv"),
          quote = FALSE, 
          row.names = FALSE)
#cms2-cms4
# Turn the results object into a tibble for use with tidyverse functions
res_tbl_2_4 <- res_2_4 %>%
        data.frame() %>%
        rownames_to_column(var="gene") %>%
        as_tibble()
# Check results output
res_tbl_2_4
# Write all results to file
write.csv(res_tbl_2_4,
          paste0("results/", clusters[1], "_", levels(cluster_metadata$group_id)[2], "_vs_", levels(cluster_metadata$group_id)[4], "_all_genes.csv"),
          quote = FALSE, 
          row.names = FALSE)
#cms3-cms4
# Turn the results object into a tibble for use with tidyverse functions
res_tbl_3_4 <- res_3_4 %>%
        data.frame() %>%
        rownames_to_column(var="gene") %>%
        as_tibble()
# Check results output
res_tbl_3_4
# Write all results to file
write.csv(res_tbl_3_4,
          paste0("results/", clusters[1], "_", levels(cluster_metadata$group_id)[3], "_vs_", levels(cluster_metadata$group_id)[4], "_all_genes.csv"),
          quote = FALSE, 
          row.names = FALSE)
         
savehistory("PseudoBulk_epi_tumor.RHistory")
##
padj_cutoff <- 0.05
sig_res_1_2 <- dplyr::filter(res_tbl_1_2, padj < padj_cutoff) %>%
        dplyr::arrange(padj)
write.csv(sig_res_1_2,
          paste0("results", clusters[1], "_", levels(cluster_metadata$group_id)[1], "_vs_", levels(cluster_metadata$group_id)[2], "_sig_genes.csv"),
          quote = FALSE, 
          row.names = FALSE)
sig_res_1_3 <- dplyr::filter(res_tbl_1_3, padj < padj_cutoff) %>%
        dplyr::arrange(padj)
write.csv(sig_res_1_3,
          paste0("results", clusters[1], "_", levels(cluster_metadata$group_id)[1], "_vs_", levels(cluster_metadata$group_id)[3], "_sig_genes.csv"),
          quote = FALSE, 
          row.names = FALSE)
sig_res_1_4 <- dplyr::filter(res_tbl_1_4, padj < padj_cutoff) %>%
        dplyr::arrange(padj)
write.csv(sig_res_1_4,
          paste0("results", clusters[1], "_", levels(cluster_metadata$group_id)[1], "_vs_", levels(cluster_metadata$group_id)[4], "_sig_genes.csv"),
          quote = FALSE, 
          row.names = FALSE)
sig_res_2_3 <- dplyr::filter(res_tbl_2_3, padj < padj_cutoff) %>%
        dplyr::arrange(padj)
write.csv(sig_res_2_3,
          paste0("results", clusters[1], "_", levels(cluster_metadata$group_id)[2], "_vs_", levels(cluster_metadata$group_id)[3], "_sig_genes.csv"),
          quote = FALSE, 
          row.names = FALSE)
sig_res_2_4 <- dplyr::filter(res_tbl_2_4, padj < padj_cutoff) %>%
        dplyr::arrange(padj)
write.csv(sig_res_2_4,
          paste0("results", clusters[1], "_", levels(cluster_metadata$group_id)[2], "_vs_", levels(cluster_metadata$group_id)[4], "_sig_genes.csv"),
          quote = FALSE, 
          row.names = FALSE)
 
sig_res_3_4 <- dplyr::filter(res_tbl_3_4, padj < padj_cutoff) %>%
        dplyr::arrange(padj)
write.csv(sig_res_3_4,
          paste0("results", clusters[1], "_", levels(cluster_metadata$group_id)[3], "_vs_", levels(cluster_metadata$group_id)[4], "_sig_genes.csv"),
          quote = FALSE, 
          row.names = FALSE)
sig_res_3_4
ls()
sig_norm <- data.frame(normalized_counts) %>%
        rownames_to_column(var = "gene") %>%
        dplyr::filter(gene %in% sig_res_1_2$gene)
dds
normalized_counts <- counts(dds, 
                            normalized = TRUE)
head(normalized_counts)
sig_norm_1_2 <- data.frame(normalized_counts) %>%
        rownames_to_column(var = "gene") %>%
        dplyr::filter(gene %in% sig_res_1_2$gene)
heat_colors <- brewer.pal(6, "YlOrRd")
pheatmap(sig_norm_1_2[ , 2:length(colnames(sig_norm_1_2))], 
    color = heat_colors, 
    cluster_rows = T, 
    show_rownames = F,
    annotation = cluster_metadata[, c("group_id", "cluster_id")], 
    border_color = NA, 
    fontsize = 10, 
    scale = "row", 
    fontsize_row = 10, 
    height = 20) 
dev.off()
dev.off()
dev.off()
pheatmap(sig_norm_1_2[ , 2:length(colnames(sig_norm_1_2))], 
    color = heat_colors, 
    cluster_rows = T, 
    show_rownames = F,
    annotation = cluster_metadata[, c("group_id", "cluster_id")], 
    border_color = NA, 
    fontsize = 10, 
    scale = "row", 
    fontsize_row = 10, 
    height = 20) 
top20_sig_genes_1_2  <- sig_res_1_2  %>%
        dplyr::arrange(padj) %>%
        dplyr::pull(gene) %>%
        head(n=20)
top20_sig_norm_1_2 <- data.frame(normalized_counts) %>%
        rownames_to_column(var = "gene") %>%
        dplyr::filter(gene %in% top20_sig_genes_1_2 )
dev.off()
gathered_top20_sig_1_2 <- top20_sig_norm_1_2 %>%
        gather(colnames(top20_sig_norm)[2:length(colnames(top20_sig_norm_1_2))], key = "samplename", value = "normalized_counts")
        
gathered_top20_sig_1_2 <- inner_join(ei[, c("sample_id", "group_id" )], gathered_top20_sig_1_2, by = c("sample_id" = "samplename"))
## plot using ggplot2
ggplot(gathered_top20_sig_1_2) +
        geom_point(aes(x = gene, 
                       y = normalized_counts, 
                       color = group_id), 
                   position=position_jitter(w=0.1,h=0)) +
        scale_y_log10() +
        xlab("Genes") +
        ylab("log10 Normalized Counts") +
        ggtitle("Top 20 Significant DE Genes for CMS1 Vs. CMS2") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
        theme(plot.title = element_text(hjust = 0.5))
gathered_top20_sig_1_2 <- top20_sig_norm_1_2 %>%
        gather(colnames(top20_sig_norm_1_2)[2:length(colnames(top20_sig_norm_1_2))], key = "samplename", value = "normalized_counts")
        
gathered_top20_sig_1_2 <- inner_join(ei[, c("sample_id", "group_id" )], gathered_top20_sig_1_2, by = c("sample_id" = "samplename"))
## plot using ggplot2
ggplot(gathered_top20_sig_1_2) +
        geom_point(aes(x = gene, 
                       y = normalized_counts, 
                       color = group_id), 
                   position=position_jitter(w=0.1,h=0)) +
        scale_y_log10() +
        xlab("Genes") +
        ylab("log10 Normalized Counts") +
        ggtitle("Top 20 Significant DE Genes for CMS1 Vs. CMS2") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
        theme(plot.title = element_text(hjust = 0.5))
####        
png("epi_tumor_cms1_csm2_scatterplot.png", units="in", width=10, height=10, res=300)
ggplot(gathered_top20_sig_1_2) +
        geom_point(aes(x = gene, 
                       y = normalized_counts, 
                       color = group_id), 
                   position=position_jitter(w=0.1,h=0)) +
        scale_y_log10() +
        xlab("Genes") +
        ylab("log10 Normalized Counts") +
        ggtitle("Top 20 Significant DE Genes for CMS1 Vs. CMS2") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
        theme(plot.title = element_text(hjust = 0.5))
dev.off()
# Run pheatmap using the metadata data frame for the annotation
pheatmap(sig_norm_1_2[ , 2:length(colnames(sig_norm_1_2))], 
    color = heat_colors, 
    cluster_rows = T, 
    show_rownames = F,
    annotation = cluster_metadata[, c("group_id", "cluster_id")], 
    border_color = NA, 
    fontsize = 10, 
    scale = "row", 
    fontsize_row = 10, 
    height = 20) 
save.image("pseudobulk_Epi_tumor_063021.RData")
savehistory("pseudobulk_Epi_tumor_063021.RHistory")


########################################################################################################################################
Using dplyr to Subset and selec padj value < 0.05
########################################################################################################################################

padj_cutoff = 0.05

clust_tgfb <- subset(markers_genes_fibro,cluster=="TGFB-myCAF")%>% dplyr::filter(p_val_adj < padj_cutoff)%>%dplyr::arrange(p_val_adj)


#########################################################################################################################################################
GO analysis
######
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library("cowplot")
library(grid)
library(ggpubr)

#TGFB-Mycaf
clust_tgfb <- subset(markers_genes_fibro,cluster=="TGFB-myCAF")%>% dplyr::filter(p_val_adj < padj_cutoff)%>%dplyr::arrange(p_val_adj) 
dim(clust_tgfb)
tgfb_gene_with_fc <- dplyr::select(clust_tgfb, gene, avg_logFC)

gene_with_fc_vector_tgfb <- tgfb_gene_with_fc[,2]
names(gene_with_fc_vector_tgfb) = as.character(tgfb_gene_with_fc[,1])
gene_with_fc_vector_tgfb = sort(gene_with_fc_vector_tgfb, decreasing = TRUE)
symbol_tgfb <- as.vector(clust_tgfb$gene)
ids_tgfb=mapIds(org.Hs.eg.db, symbol_tgfb, 'ENTREZID', 'SYMBOL')

kegg_tgfb <- enrichKEGG(ids_tgfb, organism="hsa", pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.2)
BP_tgfb <- enrichGO(ids_tgfb, OrgDb= org.Hs.eg.db, pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.2, ont="BP")
CC_tgfb <- enrichGO(ids_tgfb, OrgDb= org.Hs.eg.db, pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.2, ont="CC")
MF_tgfb <- enrichGO(ids_tgfb, OrgDb= org.Hs.eg.db, pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.2, ont="MF")


write.csv(kegg_tgfb,"kegg_tgfb.csv")
write.csv(BP_tgfb,"BP_tgfb.csv")
write.csv(CC_tgfb,"CC_tgfb.csv")
write.csv(MF_tgfb,"MF_tgfb.csv")

png("TGFB_mycaf_GO_Dotplot_with_FC.png", height = 25, width = 30, units = "in", res = 300)
plot_grid(dotplot(kegg_tgfb, showCategory=50), dotplot(BP_tgfb, showCategory=50), dotplot(CC_tgfb, showCategory=50), dotplot(MF_tgfb, showCategory=50) + rremove("x.text"), labels = c("A: KEGG", "B: Biological Process", "C: Cellular Component", "D: Molecular Functions"), ncol = 2, nrow = 2)
dev.off()


#CAFS4
clust_CAFS4 <- subset(markers_genes_fibro,cluster=="CAF-S4")%>% dplyr::filter(p_val_adj < padj_cutoff)%>%dplyr::arrange(p_val_adj)

CAFS4_gene_with_fc <- dplyr::select(clust_CAFS4, gene, avg_logFC)

gene_with_fc_vector_CAFS4 <- CAFS4_gene_with_fc[,2]
names(gene_with_fc_vector_CAFS4) = as.character(CAFS4_gene_with_fc[,1])
gene_with_fc_vector_CAFS4 = sort(gene_with_fc_vector_CAFS4, decreasing = TRUE)
symbol_CAFS4 <- as.vector(clust_CAFS4$gene)
ids_CAFS4=mapIds(org.Hs.eg.db, symbol_CAFS4, 'ENTREZID', 'SYMBOL')

kegg_CAFS4 <- enrichKEGG(ids_CAFS4, organism="hsa", pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.2)
BP_CAFS4 <- enrichGO(ids_CAFS4, OrgDb= org.Hs.eg.db, pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.2, ont="BP")
CC_CAFS4 <- enrichGO(ids_CAFS4, OrgDb= org.Hs.eg.db, pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.2, ont="CC")
MF_CAFS4 <- enrichGO(ids_CAFS4, OrgDb= org.Hs.eg.db, pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.2, ont="MF")

write.csv(kegg_CAFS4,"kegg_CAFS4.csv")
write.csv(BP_CAFS4,"BP_CAFS4.csv")
write.csv(CC_CAFS4,"CC_CAFS4.csv")
write.csv(MF_CAFS4,"MF_tgfb.csv")

png("CAFS4_GO_Dotplot_with_FC.png", height = 25, width = 30, units = "in", res = 300)
plot_grid(dotplot(kegg_CAFS4, showCategory=50), dotplot(BP_CAFS4, showCategory=50), dotplot(CC_CAFS4, showCategory=50),dotplot(MF_CAFS4, showCategory=50) + rremove("x.text"), labels = c("A: KEGG", "B: Biological Process", "C: Cellular Component", "D: Molecular Functions"), ncol = 2, nrow = 2)
dev.off()


#detox-IL-iCAF
clust_detoxILiCAF <- subset(markers_genes_fibro,cluster=="detox-IL-iCAF")%>% dplyr::filter(p_val_adj < padj_cutoff)%>%dplyr::arrange(p_val_adj)
dim(clust_detoxILiCAF)
detoxILiCAF_gene_with_fc <- dplyr::select(clust_detoxILiCAF, gene, avg_logFC)

gene_with_fc_vector_detoxILiCAF <- detoxILiCAF_gene_with_fc[,2]
names(gene_with_fc_vector_detoxILiCAF) = as.character(detoxILiCAF_gene_with_fc[,1])
gene_with_fc_vector_detoxILiCAF = sort(gene_with_fc_vector_detoxILiCAF, decreasing = TRUE)
symbol_detoxILiCAF <- as.vector(clust_detoxILiCAF$gene)
ids_detoxILiCAF=mapIds(org.Hs.eg.db, symbol_detoxILiCAF, 'ENTREZID', 'SYMBOL')

kegg_detoxILiCAF <- enrichKEGG(ids_detoxILiCAF, organism="hsa", pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.2)
BP_detoxILiCAF <- enrichGO(ids_detoxILiCAF, OrgDb= org.Hs.eg.db, pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.2, ont="BP")
CC_detoxILiCAF <- enrichGO(ids_detoxILiCAF, OrgDb= org.Hs.eg.db, pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.2, ont="CC")
MF_detoxILiCAF <- enrichGO(ids_detoxILiCAF, OrgDb= org.Hs.eg.db, pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.2, ont="MF")


write.csv(kegg_detoxILiCAF,"kegg_detoxILiCAF.csv")
write.csv(BP_detoxILiCAF,"BP_detoxILiCAF.csv")
write.csv(CC_detoxILiCAF,"CC_detoxILiCAF.csv")
write.csv(MF_detoxILiCAF,"MF_detoxILiCAF.csv")

png("detoxILiCAF_GO_Dotplot_with_FC.png", height = 25, width = 30, units = "in", res = 300)
plot_grid(dotplot(kegg_detoxILiCAF, showCategory=50), dotplot(BP_detoxILiCAF, showCategory=50), dotplot(CC_detoxILiCAF, showCategory=50), dotplot(MF_detoxILiCAF, showCategory=50) + rremove("x.text"), labels = c("A: KEGG", "B: Biological Process", "C: Cellular Component", "D: Molecular Functions"), ncol = 2, nrow = 2)
dev.off()


#ecm-myCAF
clust_ecmmyCAF <- subset(markers_genes_fibro,cluster=="ecm-myCAF")%>% dplyr::filter(p_val_adj < padj_cutoff)%>%dplyr::arrange(p_val_adj)
dim(clust_ecmmyCAF)
ecmmyCAF_gene_with_fc <- dplyr::select(clust_ecmmyCAF, gene, avg_logFC)

gene_with_fc_vector_ecmmyCAF <- ecmmyCAF_gene_with_fc[,2]
names(gene_with_fc_vector_ecmmyCAF) = as.character(ecmmyCAF_gene_with_fc[,1])
gene_with_fc_vector_ecmmyCAF = sort(gene_with_fc_vector_ecmmyCAF, decreasing = TRUE)
symbol_ecmmyCAF <- as.vector(clust_ecmmyCAF$gene)
ids_ecmmyCAF=mapIds(org.Hs.eg.db, symbol_ecmmyCAF, 'ENTREZID', 'SYMBOL')

kegg_ecmmyCAF <- enrichKEGG(ids_ecmmyCAF, organism="hsa", pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.2)
BP_ecmmyCAF <- enrichGO(ids_ecmmyCAF, OrgDb= org.Hs.eg.db, pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.2, ont="BP")
CC_ecmmyCAF <- enrichGO(ids, OrgDb= org.Hs.eg.db, pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.2, ont="CC")
MF_ecmmyCAF <- enrichGO(ids_ecmmyCAF, OrgDb= org.Hs.eg.db, pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.2, ont="MF")


write.csv(kegg_ecmmyCAF,"kegg_ecmmyCAF.csv")
write.csv(BP_ecmmyCAF,"BP_ecmmyCAF.csv")
write.csv(CC_ecmmyCAF,"CC_ecmmyCAF.csv")
write.csv(MF_ecmmyCAF,"MF_ecmmyCAF.csv")

png("ecmmyCAF_GO_Dotplot_with_FC.png", height = 25, width = 30, units = "in", res = 300)
plot_grid(dotplot(kegg_ecmmyCAF, showCategory=50), dotplot(BP_ecmmyCAF, showCategory=50), dotplot(CC_ecmmyCAF, showCategory=50) , dotplot(MF_ecmmyCAF, showCategory=50) + rremove("x.text"), labels = c("A: KEGG", "B: Biological Process", "C: Cellular Component", "D: Molecular Functions"), ncol = 2, nrow = 2)
dev.off()


#wound-myCAF
clust_woundmyCAF <- subset(markers_genes_fibro,cluster=="wound-myCAF")%>% dplyr::filter(p_val_adj < padj_cutoff)%>%dplyr::arrange(p_val_adj)
dim(clust_woundmyCAF)
woundmyCAF_gene_with_fc <- dplyr::select(clust_woundmyCAF, gene, avg_logFC)

gene_with_fc_vector_woundmyCAF <- woundmyCAF_gene_with_fc[,2]
names(gene_with_fc_vector_woundmyCAF) = as.character(woundmyCAF_gene_with_fc[,1])
gene_with_fc_vector_woundmyCAF = sort(gene_with_fc_vector_woundmyCAF, decreasing = TRUE)
symbol_woundmyCAF <- as.vector(clust_woundmyCAF$gene)
ids_woundmyCAF=mapIds(org.Hs.eg.db, symbol_woundmyCAF, 'ENTREZID', 'SYMBOL')

kegg_woundmyCAF <- enrichKEGG(ids_woundmyCAF, organism="hsa", pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.2)
BP_woundmyCAF <- enrichGO(ids_woundmyCAF, OrgDb= org.Hs.eg.db, pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.2, ont="BP")
CC_woundmyCAF <- enrichGO(ids_woundmyCAF, OrgDb= org.Hs.eg.db, pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.2, ont="CC")
MF_woundmyCAF <- enrichGO(ids_woundmyCAF, OrgDb= org.Hs.eg.db, pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.2, ont="MF")


write.csv(kegg_woundmyCAF,"kegg_woundmyCAF.csv")
write.csv(BP_woundmyCAF,"BP_woundmyCAF.csv")
write.csv(CC_woundmyCAF,"CC_woundmyCAF.csv")
write.csv(MF_woundmyCAF,"MF_woundmyCAF.csv")

png("woundmyCAF_GO_Dotplot_with_FC.png", height = 25, width = 20, units = "in", res = 300)
plot_grid(dotplot(kegg_woundmyCAF, showCategory=50), dotplot(BP_woundmyCAF, showCategory=50), dotplot(CC_woundmyCAF, showCategory=50), dotplot(MF_woundmyCAF, showCategory=50) + rremove("x.text"), labels = c("A: KEGG", "B: Biological Process", "C: Cellular Component", "D: Molecular Functions"), ncol = 2, nrow = 2)
dev.off()


# for all ontologies in one

clust_woundmyCAF_t100 <- subset(markers_genes_fibro,cluster=="wound-myCAF")%>% dplyr::filter(p_val_adj < padj_cutoff)%>%dplyr::arrange(p_val_adj)  %>% head(100)
woundmyCAF_gene_with_fc_t100  <- dplyr::select(clust_woundmyCAF_t100, gene, avg_logFC)
gene_with_fc_vector_woundmyCAF_t100  <- woundmyCAF_gene_with_fc_t100 [,2]
names(gene_with_fc_vector_woundmyCAF_t100 ) = as.character(woundmyCAF_gene_with_fc_t100 [,1])
gene_with_fc_vector_woundmyCAF_t100  = sort(gene_with_fc_vector_woundmyCAF_t100 , decreasing = TRUE)
symbol_woundmyCAF_t100  <- as.vector(clust_woundmyCAF_t100$gene)
ids_woundmyCAF_t100 =mapIds(org.Hs.eg.db, symbol_woundmyCAF_t100, 'ENTREZID', 'SYMBOL')

woundmyCAF_all <- enrichGO(ids_woundmyCAF_t100, OrgDb = "org.Hs.eg.db", pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.2,ont="all")
kegg_woundmyCAF_t100 <- enrichKEGG(ids_woundmyCAF_t100, organism="hsa", pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.2)
write.csv(woundmyCAF_all,"woundmyCAF_all.csv")
write.csv(kegg_woundmyCAF_t100,"woundmyCAF_kegg.csv")

png("woundmyCAF.png", height = 20, width = 30, units = "in", res = 300)
plot_grid(dotplot(kegg_woundmyCAF_t100, showCategory=50), dotplot(woundmyCAF_all, split="ONTOLOGY", showCategory=30) + facet_grid(ONTOLOGY~., scale="free") + rremove("x.text"), labels = c("A: KEGG - wound-myCAF ", "B: Gene Ontologies - wound-myCAF"), ncol = 2, nrow = 1)
dev.off()



#########
Fig Publication
#########



png("tsne_all_clusters.png", units="in", width=10, height=8, res=300)
dittoDimPlot(alldata_new, "seurat_clusters_new", reduction.use ="TSNE_on_CCA",color.panel = c("#d1001c", "#9900ff", "#ff6600", "#d5b60a", "#0881d1","#014600"),do.label = FALSE, labels.repel = FALSE,labels.highlight = TRUE,legend.show = TRUE, labels.size = 3, show.axes.numbers = TRUE, xlab=NULL,ylab=NULL,main=NULL)
dev.off()

svg("tsne_all_clusters.svg", width=7, height=7, res=300)
dittoDimPlot(alldata_new, "seurat_clusters_new", reduction.use ="TSNE_on_CCA",color.panel = c("#d1001c", "#9900ff", "#ff6600", "#d5b60a", "#0881d1","#014600"),do.label = TRUE, labels.repel = FALSE,labels.highlight = TRUE,legend.show = FALSE, labels.size = 3, show.axes.numbers = FALSE, xlab=NULL,ylab=NULL,main=NULL)
dev.off()


png("tsne_condition_1.png",units="in", width=7, height=7, res=300)
dittoDimPlot(alldata_new, "Condition", reduction.use ="TSNE_on_CCA",color.panel = c("#d1001c","#014600"),do.label = FALSE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = TRUE, labels.size = 3, show.axes.numbers = FALSE, xlab=NULL,ylab=NULL,main=NULL)
dev.off()


svg("tsne_condition_1.svg", width=7, height=7)
dittoDimPlot(alldata_new, "Condition", reduction.use ="TSNE_on_CCA",color.panel = c("#d1001c","#014600"),do.label = FALSE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = TRUE, labels.size = 3, show.axes.numbers = FALSE, xlab=NULL,ylab=NULL,main=NULL)
dev.off()



png("tsne_msi.png", units="in", width=7, height=7, res=300)
dittoDimPlot(alldata_new, "MSI_Status", reduction.use ="TSNE_on_CCA",color.panel = c("#d1001c","#014600"),do.label = FALSE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = FALSE, labels.size = 3, show.axes.numbers = FALSE, xlab=NULL,ylab=NULL,main=NULL)
dev.off()


svg("tsne_msi.svg", width=7, height=7)
dittoDimPlot(alldata_new, "MSI_Status", reduction.use ="TSNE_on_CCA",color.panel = c("#d1001c","#014600"),do.label = FALSE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = FALSE, labels.size = 3, show.axes.numbers = FALSE, xlab=NULL,ylab=NULL,main=NULL)
dev.off()


png("tsne_msi_1.png",units="in", width=7, height=7, res=300)
dittoDimPlot(alldata_new, "MSI_Status", reduction.use ="TSNE_on_CCA",color.panel = c("#d1001c","#014600"),do.label = FALSE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = TRUE, labels.size = 3, show.axes.numbers = FALSE, xlab=NULL,ylab=NULL,main=NULL)
dev.off()


svg("tsne_msi_1.svg", width=7, height=7)
dittoDimPlot(alldata_new, "MSI_Status", reduction.use ="TSNE_on_CCA",color.panel = c("#d1001c","#014600"),do.label = FALSE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = TRUE, labels.size = 3, show.axes.numbers = FALSE, xlab=NULL,ylab=NULL,main=NULL)
dev.off()


png("Alldata_markergenes_compartment.png", units="in", width=10, height=6, res=300)
cowplot::plot_grid(ncol = 6,
FeaturePlot(alldata,features = "KRT18" ,min.cutoff="q9", cols=c("lightgrey", "#014600"), label=FALSE,reduction="TSNE_on_CCA")+ labs(title = "Epithelial Cells"),
FeaturePlot(alldata,features = "CD3D" ,min.cutoff="q9", cols=c("lightgrey","#014600"), label=FALSE,reduction="TSNE_on_CCA")+ labs(title = "T Cells"),
FeaturePlot(alldata,features = "CD79A" ,min.cutoff="q9", cols=c("lightgrey","#014600"), label=FALSE,reduction="TSNE_on_CCA")+ labs(title = "B Cells"),
FeaturePlot(alldata,features = "LYZ" ,min.cutoff="q9", cols=c("lightgrey", "#014600"), label=FALSE,reduction="TSNE_on_CCA")+ labs(title = "Myeloid"),
FeaturePlot(alldata,features = "COL1A1" ,min.cutoff="q9", cols=c("lightgrey", "#014600"), label=FALSE,reduction="TSNE_on_CCA")+ labs(title = "Fibroblast"),
FeaturePlot(alldata,features = "CLDN5" ,min.cutoff="q9", cols=c("lightgrey","#014600"), label=FALSE,reduction="TSNE_on_CCA")+ labs(title = "Endothelial"))
dev.off()

png("Alldata_markergenes_compartment1.png", units="in", width=15, height=8, res=300)
cowplot::plot_grid(ncol = 3,
dittoDimPlot(alldata_new, "KRT18",reduction.use="TSNE_on_CCA",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = FALSE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = TRUE, labels.size = 3, show.axes.numbers = TRUE, xlab=NULL,ylab=NULL,main=NULL)+ labs(title = "Epithelial Cells: KRT18"),
dittoDimPlot(alldata_new, "CD3D",reduction.use="TSNE_on_CCA",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = FALSE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = TRUE, labels.size = 3, show.axes.numbers = TRUE, xlab=NULL,ylab=NULL,main=NULL)+ labs(title = "T Cells: CD3D"),
dittoDimPlot(alldata_new, "CD79A",reduction.use="TSNE_on_CCA",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = FALSE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = TRUE, labels.size = 3, show.axes.numbers = TRUE, xlab=NULL,ylab=NULL,main=NULL)+ labs(title = "B Cells: CD79A"),
dittoDimPlot(alldata_new, "LYZ",reduction.use="TSNE_on_CCA",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = FALSE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = TRUE, labels.size = 3, show.axes.numbers = TRUE, xlab=NULL,ylab=NULL,main=NULL)+ labs(title = "Myeloid Cells: LYZ"),
dittoDimPlot(alldata_new, "COL1A2",reduction.use="TSNE_on_CCA",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = FALSE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = TRUE, labels.size = 3, show.axes.numbers = TRUE, xlab=NULL,ylab=NULL,main=NULL)+ labs(title = "Fibroblast: COL1A2"),
dittoDimPlot(alldata_new, "CLDN5",reduction.use="TSNE_on_CCA",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = FALSE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = TRUE, labels.size = 3, show.axes.numbers = TRUE, xlab=NULL,ylab=NULL,main=NULL)+ labs(title = "Endothelial Cells: CLDN5"))
dev.off()





svg("Alldata_epi.svg", width=8, height=8)
dittoDimPlot(alldata_new, "KRT18",reduction.use="TSNE_on_CCA",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = FALSE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = TRUE, labels.size = 3, show.axes.numbers = TRUE, xlab=NULL,ylab=NULL,main=NULL)+ labs(title = "Epithelial Cells: KRT18")
dev.off()

svg("Alldata_tcell.svg", width=8, height=8)
dittoDimPlot(alldata_new, "CD3D",reduction.use="TSNE_on_CCA",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = FALSE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = TRUE, labels.size = 3, show.axes.numbers = TRUE, xlab=NULL,ylab=NULL,main=NULL)+ labs(title = "T Cells: CD3D")
dev.off()

svg("Alldata_bcell.svg", width=8, height=8)
dittoDimPlot(alldata_new, "CD79A",reduction.use="TSNE_on_CCA",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = FALSE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = TRUE, labels.size = 3, show.axes.numbers = TRUE, xlab=NULL,ylab=NULL,main=NULL)+ labs(title = "B Cells: CD79A")
dev.off()

svg("Alldata_myeloid.svg", width=8, height=8)
dittoDimPlot(alldata_new, "LYZ",reduction.use="TSNE_on_CCA",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = FALSE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = TRUE, labels.size = 3, show.axes.numbers = TRUE, xlab=NULL,ylab=NULL,main=NULL)+ labs(title = "Myeloid Cells: LYZ")
dev.off()

svg("Alldata_fibro.svg", width=8, height=8)
dittoDimPlot(alldata_new, "COL1A2",reduction.use="TSNE_on_CCA",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = FALSE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = TRUE, labels.size = 3, show.axes.numbers = TRUE, xlab=NULL,ylab=NULL,main=NULL)+ labs(title = "Fibroblast: COL1A2")
dev.off()

svg("Alldata_endo.svg", width=8, height=8)
dittoDimPlot(alldata_new, "CLDN5",reduction.use="TSNE_on_CCA",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = FALSE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = TRUE, labels.size = 3, show.axes.numbers = TRUE, xlab=NULL,ylab=NULL,main=NULL)+ labs(title = "Endothelial Cells: CLDN5")
dev.off()




png("Alldata_markergenes_compartment2.png", units="in", width=15, height=8, res=300)
cowplot::plot_grid(ncol = 3,
dittoDimPlot(alldata_new, "KRT18",reduction.use="TSNE_on_CCA",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = FALSE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = TRUE, labels.size = 3, show.axes.numbers = TRUE, xlab=NULL,ylab=NULL,main=NULL),
dittoDimPlot(alldata_new, "CD3D",reduction.use="TSNE_on_CCA",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = FALSE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = TRUE, labels.size = 3, show.axes.numbers = TRUE, xlab=NULL,ylab=NULL,main=NULL),
dittoDimPlot(alldata_new, "CD79A",reduction.use="TSNE_on_CCA",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = FALSE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = TRUE, labels.size = 3, show.axes.numbers = TRUE, xlab=NULL,ylab=NULL,main=NULL),
dittoDimPlot(alldata_new, "LYZ",reduction.use="TSNE_on_CCA",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = FALSE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = TRUE, labels.size = 3, show.axes.numbers = TRUE, xlab=NULL,ylab=NULL,main=NULL),
dittoDimPlot(alldata_new, "COL1A2",reduction.use="TSNE_on_CCA",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = FALSE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = TRUE, labels.size = 3, show.axes.numbers = TRUE, xlab=NULL,ylab=NULL,main=NULL),
dittoDimPlot(alldata_new, "CLDN5",reduction.use="TSNE_on_CCA",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = FALSE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = TRUE, labels.size = 3, show.axes.numbers = TRUE, xlab=NULL,ylab=NULL,main=NULL))
dev.off()



svg("Alldata_markergenes_compartment.svg", width=8, height=5)
cowplot::plot_grid(ncol = 3,
dittoDimPlot(alldata_new, "KRT18",reduction.use="TSNE_on_CCA",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = FALSE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = FALSE, labels.size = 3, show.axes.numbers = FALSE, xlab=NULL,ylab=NULL,main=NULL)+ labs(title = "Epithelial Cells: KRT18"),
dittoDimPlot(alldata_new, "CD3D",reduction.use="TSNE_on_CCA",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = FALSE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = FALSE, labels.size = 3, show.axes.numbers = FALSE, xlab=NULL,ylab=NULL,main=NULL)+ labs(title = "T Cells: CD3D"),
dittoDimPlot(alldata_new, "CD79A",reduction.use="TSNE_on_CCA",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = FALSE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = FALSE, labels.size = 3, show.axes.numbers = FALSE, xlab=NULL,ylab=NULL,main=NULL)+ labs(title = "B Cells: CD79A"),
dittoDimPlot(alldata_new, "LYZ",reduction.use="TSNE_on_CCA",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = FALSE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = FALSE, labels.size = 3, show.axes.numbers = FALSE, xlab=NULL,ylab=NULL,main=NULL)+ labs(title = "Myeloid Cells: LYZ"),
dittoDimPlot(alldata_new, "COL1A2",reduction.use="TSNE_on_CCA",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = FALSE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = FALSE, labels.size = 3, show.axes.numbers = FALSE, xlab=NULL,ylab=NULL,main=NULL)+ labs(title = "Fibroblast: COL1A2"),
dittoDimPlot(alldata_new, "CLDN5",reduction.use="TSNE_on_CCA",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = FALSE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = FALSE, labels.size = 3, show.axes.numbers = FALSE, xlab=NULL,ylab=NULL,main=NULL)+ labs(title = "Endothelial Cells: CLDN5"))
dev.off()

png("Alldata_barplots.png", units="in", width=15, height=7, res=300)
cowplot::plot_grid(ncol = 2,
dittoBarPlot(alldata_new, "seurat_clusters_new", group.by = "orig.ident",color.panel = c("#d1001c", "#9900ff", "#ff6600", "#d5b60a", "#0881d1","#014600"),main=NULL),
dittoBarPlot(alldata_new, "seurat_clusters_new", group.by = "bulk_prediction",color.panel = c("#d1001c", "#9900ff", "#ff6600", "#d5b60a", "#0881d1","#014600"),main=NULL))
dev.off()


svg("Alldata_barplots.svg", width=15, height=7,pointsize=12)
cowplot::plot_grid(ncol = 2,
dittoBarPlot(alldata_new, "seurat_clusters_new", group.by = "orig.ident",color.panel = c("#d1001c", "#9900ff", "#ff6600", "#d5b60a", "#0881d1","#014600"),main=NULL),
dittoBarPlot(alldata_new, "seurat_clusters_new", group.by = "bulk_prediction",color.panel = c("#d1001c", "#9900ff", "#ff6600", "#d5b60a", "#0881d1","#014600"),main=NULL))
dev.off()



#Lee Data 


svg("Lee_barplots_allcounts.svg", width=10, height=10)
dittoBarPlot(kor, "Cell_type", group.by = "orig.ident",color.panel = c("#d1001c", "#9900ff", "#ff6600", "#d5b60a", "#0881d1","#014600"),main=NULL)
dev.off()

immune/stromal --> T,B, myeloid, stromal and mast



png("kor_barplots.png", units="in", width=15, height=10, res=300)
cowplot::plot_grid(ncol = 3,
dittoDimPlot(kor, "seurat_clusters", reduction.use ="tsne",do.label = TRUE, labels.repel = FALSE,labels.highlight = TRUE,legend.show = TRUE, labels.size = 3, show.axes.numbers = TRUE, xlab=NULL,ylab=NULL,main=NULL),
dittoDimPlot(kor, "KRT18",reduction.use="tsne",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = FALSE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = FALSE, labels.size = 3, show.axes.numbers = FALSE, xlab=NULL,ylab=NULL,main=NULL)+ labs(title = "Epithelial Cells: KRT18"),
dittoDimPlot(kor, "CD3D",reduction.use="tsne",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = FALSE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = FALSE, labels.size = 3, show.axes.numbers = FALSE, xlab=NULL,ylab=NULL,main=NULL)+ labs(title = "T Cells: CD3D"),
dittoDimPlot(kor, "CD79A",reduction.use="tsne",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = FALSE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = FALSE, labels.size = 3, show.axes.numbers = FALSE, xlab=NULL,ylab=NULL,main=NULL)+ labs(title = "B Cells: CD79A"),
dittoDimPlot(kor, "LYZ",reduction.use="tsne",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = FALSE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = FALSE, labels.size = 3, show.axes.numbers = FALSE, xlab=NULL,ylab=NULL,main=NULL)+ labs(title = "Myeloid Cells: LYZ"),
dittoDimPlot(kor, "COL1A2",reduction.use="tsne",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = FALSE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = FALSE, labels.size = 3, show.axes.numbers = FALSE, xlab=NULL,ylab=NULL,main=NULL)+ labs(title = "Fibroblast: COL1A2"),
dittoDimPlot(kor, "CLDN5",reduction.use="tsne",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = FALSE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = FALSE, labels.size = 3, show.axes.numbers = FALSE, xlab=NULL,ylab=NULL,main=NULL)+ labs(title = "Endothelial Cells: CLDN5"),
dittoDimPlot(kor, "KIT",reduction.use="tsne",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = FALSE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = FALSE, labels.size = 3, show.axes.numbers = FALSE, xlab=NULL,ylab=NULL,main=NULL)+ labs(title = "MAST Cells: KIT"))
dev.off()

#Renaming the clusters of Lee-korean --> cell Types

old.cluster.ids <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25)
new.cluster.ids <- c("T cells","T cells", "B cells", "B cells","T cells","Myeloid", "T cells", "Epithelial" ,"T cells","T cells","B cells","T cells","Myeloid","Fibroblast","Fibroblast","B cells","Endothelial","T cells","Mast cells","T cells","Fibroblast","Myeloid","Mast cells","T cells","B cells","Endothelial")
kor$seurat_clusters <- plyr::mapvalues(x = kor$seurat_clusters, from = old.cluster.ids, to = new.cluster.ids)



kable(table(tcells_new@meta.data$seurat_clusters))

# subset the clusters

|Var1        |  Freq|
|:-----------|-----:|
|T cells     | 17150|
|B cells     |  7199|
|Myeloid     |  3113|
|Epithelial  |  1826|
|Fibroblast  |  1364|
|Endothelial |   388|
|Mast cells  |   343|



kor.immune <- SubsetData(kor, subset.name = "seurat_clusters", accept.value = c("B cells","T cells","Myeloid","Mast cells"))

kor.stromal<- SubsetData(kor, subset.name = "seurat_clusters", accept.value = c("Fibroblast","Endothelial"))

#8fce00 == light Green (Mast cells)
#014600 == Dark Green (T cells)
#d1001c == Dark Red (B cells)
#0881d1 == Dark Blue (Myeloid)
#d5b60a == cyan (methi color)(Epithelial)
#ff6600 == Peach (endothelial cells)
#9900ff == CAFs


c("#d1001c", "#9900ff", "#ff6600", "#d5b60a", "#0881d1","#014600")

svg("Kor_tsne_all_clusters.svg", width=7, height=7)
dittoDimPlot(kor, "seurat_clusters", reduction.use ="tsne",color.panel = c("#014600", "#d1001c", "#0881d1", "#d5b60a", "#9900ff","#ff6600","#8fce00"),do.label = FALSE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = TRUE, labels.size = 3, show.axes.numbers = FALSE, xlab=NULL,ylab=NULL,main=NULL)
dev.off()

svg("Kor_Barplots_all_clusters.svg", width=7, height=7)
dittoBarPlot(kor, "seurat_clusters", group.by = "orig.ident",color.panel = c("#d1001c", "#ff6600", "#d5b60a", "#0881d1", "#8fce00","#9900ff","#014600"),main=NULL)
dev.off()


svg("Kor__immune_Barplots.svg", width=7, height=7)
dittoBarPlot(kor.immune, "seurat_clusters", group.by = "orig.ident",color.panel = c("#d1001c", "#8fce00", "#0881d1","#014600"),main=NULL)
dev.off()


svg("Kor_stroma_Barplots.svg", width=7, height=7)
dittoBarPlot(kor.stromal, "seurat_clusters", group.by = "orig.ident",color.panel = c("#ff6600","#9900ff"),main=NULL)
dev.off()


# kor marker genes


svg("kor_markers.svg", width=10, height=10)
png("kor_markers.png", units="in", width=15, height=10, res=300)
cowplot::plot_grid(ncol = 3,
dittoDimPlot(kor, "KRT18",reduction.use="tsne",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = TRUE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = TRUE, labels.size = 3, show.axes.numbers = TRUE, xlab=NULL,ylab=NULL,main=NULL)+ labs(title = "Epithelial Cells: KRT18"),
dittoDimPlot(kor, "CD3D",reduction.use="tsne",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = TRUE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = TRUE, labels.size = 3, show.axes.numbers = TRUE, xlab=NULL,ylab=NULL,main=NULL)+ labs(title = "T Cells: CD3D"),
dittoDimPlot(kor, "CD79A",reduction.use="tsne",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = TRUE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = TRUE, labels.size = 3, show.axes.numbers = TRUE, xlab=NULL,ylab=NULL,main=NULL)+ labs(title = "B Cells: CD79A"),
dittoDimPlot(kor, "LYZ",reduction.use="tsne",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = TRUE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = TRUE, labels.size = 3, show.axes.numbers = TRUE, xlab=NULL,ylab=NULL,main=NULL)+ labs(title = "Myeloid Cells: LYZ"),
dittoDimPlot(kor, "COL1A2",reduction.use="tsne",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = TRUE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = TRUE, labels.size = 3, show.axes.numbers = TRUE, xlab=NULL,ylab=NULL,main=NULL)+ labs(title = "Fibroblast: COL1A2"),
dittoDimPlot(kor, "CLDN5",reduction.use="tsne",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = TRUE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = TRUE, labels.size = 3, show.axes.numbers = TRUE, xlab=NULL,ylab=NULL,main=NULL)+ labs(title = "Endothelial Cells: CLDN5"))
dev.off()

svg("kor_conditions.svg", width=7, height=5)
dittoBarPlot(kor, "seurat_clusters", group.by = "Condition",color.panel = c("#d1001c", "#ff6600", "#d5b60a", "#9900ff", "#8fce00","#0881d1","#014600"),main=NULL)
dev.off()

#nearestCMS.RF
svg("kor_cms.svg", width=7, height=5)
dittoBarPlot(kor, "seurat_clusters", group.by = "nearestCMS.RF",color.panel = c("#d1001c", "#ff6600", "#d5b60a", "#9900ff", "#8fce00","#0881d1","#014600"),main=NULL)
dev.off()

svg("kor_mast.svg", width=7, height=5)
dittoDimPlot(kor, "KIT",reduction.use="tsne",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = FALSE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = FALSE, labels.size = 3, show.axes.numbers = FALSE, xlab=NULL,ylab=NULL,main=NULL)+ labs(title = "MAST Cells: KIT")
dev.off()

# yuan

png("EPi_namf_markers.png", units="in", width=20, height=15, res=300)
multi_dittoDimPlot(epi_usage1, c("EPCAM","KRT8","KRT18","COL1A1","COL1A2","CLDN5","CD3D","CD79A"))
dev.off()

# yuan marker genes


svg("yuan_markers.svg", width=10, height=10)
#png("yaun_markers.png", units="in", width=15, height=10, res=300)
cowplot::plot_grid(ncol = 3,
dittoDimPlot(yuan_cr, "KRT18",reduction.use="tsne",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = TRUE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = TRUE, labels.size = 3, show.axes.numbers = TRUE, xlab=NULL,ylab=NULL,main=NULL)+ labs(title = "Epithelial Cells: KRT18"),
dittoDimPlot(yuan_cr, "CD3D",reduction.use="tsne",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = TRUE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = TRUE, labels.size = 3, show.axes.numbers = TRUE, xlab=NULL,ylab=NULL,main=NULL)+ labs(title = "T Cells: CD3D"),
dittoDimPlot(yuan_cr, "CD79A",reduction.use="tsne",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = TRUE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = TRUE, labels.size = 3, show.axes.numbers = TRUE, xlab=NULL,ylab=NULL,main=NULL)+ labs(title = "B Cells: CD79A"),
dittoDimPlot(yuan_cr, "LYZ",reduction.use="tsne",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = TRUE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = TRUE, labels.size = 3, show.axes.numbers = TRUE, xlab=NULL,ylab=NULL,main=NULL)+ labs(title = "Myeloid Cells: LYZ"),
dittoDimPlot(yuan_cr, "COL3A",reduction.use="tsne",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = TRUE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = TRUE, labels.size = 3, show.axes.numbers = TRUE, xlab=NULL,ylab=NULL,main=NULL)+ labs(title = "Fibroblast: COL3A1"),
dittoDimPlot(yuan_cr, "CLDN5",reduction.use="tsne",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = TRUE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = TRUE, labels.size = 3, show.axes.numbers = TRUE, xlab=NULL,ylab=NULL,main=NULL)+ labs(title = "Endothelial Cells: CLDN5"))
dev.off()

dittoDimPlot(yuan_cr, "COL3A1",reduction.use="tsne",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = TRUE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = TRUE, labels.size = 3, show.axes.numbers = TRUE, xlab=NULL,ylab=NULL,main=NULL)

png("yaun_dimplot.png", units="in", width=10, height=10, res=300)
dittoDimPlot(yuan_cr, "seurat_clusters", reduction.use ="tsne",do.label = TRUE, labels.repel = TRUE,labels.highlight = FALSE,legend.show = TRUE, labels.size = 3, show.axes.numbers = FALSE, xlab=NULL,ylab=NULL,main=NULL)
dev.off()



#Renaming the clusters of yuan--> cell Types

old.cluster.ids <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33)
new.cluster.ids <- c("T cells","B cells", "T cells", "T cells","T cells","Epithelial", "T cells", "Fibroblast" ,"T cells","Endothelial","Fibroblast","B cells","Fibroblast","Fibroblast","Epithelial","Epithelial","Endothelial","Endothelial","B cells","Myeloid","Epithelial","Fibroblast","Myeloid","Fibroblast","Epithelial","Fibroblast","Fibroblast","Myeloid","B cells","B cells","Fibroblast","Epithelial","Epithelial")

yuan_cr$seurat_clusters <- plyr::mapvalues(x = yuan_cr$seurat_clusters, from = old.cluster.ids, to = new.cluster.ids)

yuan.meta <- read.csv("yuanmetadata.csv", header=TRUE, sep = ',',row.names=1)

# add metadata to the seurat object

yuan_cr <- AddMetaData(yuan_cr, yuan.meta)


yuan.immune <- SubsetData(yuan_cr, subset.name = "seurat_clusters", accept.value = c("B cells","T cells","Myeloid"))

yuan.stromal<- SubsetData(yuan_cr, subset.name = "seurat_clusters", accept.value = c("Fibroblast","Endothelial"))


#8fce00 == light Green (Mast cells)
#014600 == Dark Green (T cells)
#d1001c == Dark Red (B cells)
#0881d1 == Dark Blue (Myeloid)
#d5b60a == cyan (methi color)(Epithelial)
#ff6600 == Peach (endothelial cells)
#9900ff == CAFs


svg("yuan_tsne_all_clusters.svg", width=7, height=7)
dittoDimPlot(yuan_cr, "seurat_clusters", reduction.use ="tsne",color.panel = c("#014600", "#d1001c", "#d5b60a", "#9900ff", "#ff6600","#0881d1"),do.label = FALSE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = TRUE, labels.size = 3, show.axes.numbers = FALSE, xlab=NULL,ylab=NULL,main=NULL)
dev.off()

svg("yuan_Barplots_all_clusters.svg", width=7, height=7)
dittoBarPlot(yuan_cr, "seurat_clusters", group.by = "orig.ident",color.panel = c("#d1001c", "#ff6600", "#d5b60a", "#9900ff", "#0881d1","#014600"),main=NULL)
dev.off()

#####
svg("kor_Barplots_all_clusters.svg", width=7, height=7)
dittoBarPlot(kor_clean_comb, "celltype_new", group.by = "orig.ident",main=NULL, color.panel = c("#d1001c", "#d5b60a", "#9900ff", "#ff6600","#0881d1","#014600"))
dev.off()


svg("kor_immune_Barplots.svg", width=7, height=7)
dittoBarPlot(kor_clean_comb.immune, "celltype_new", group.by = "orig.ident",color.panel = c("#d1001c", "#0881d1", "#014600"),main=NULL)
dev.off()

svg("kor_stroma_Barplots.svg", width=7, height=7)
dittoBarPlot(kor.fibro, "seurat_clusters", group.by = "orig.ident",color.panel = c("#ff6600","#9900ff"),main=NULL)
dev.off()

svg("kor_conditions.svg", width=7, height=5)
dittoBarPlot(kor_clean_comb, "celltype_new", group.by = "Condition",color.panel = c("#d1001c", "#d5b60a", "#9900ff", "#ff6600","#0881d1","#014600"),main=NULL)
dev.off()

#nearestCMS.RF

svg("kor_cms.svg", width=7, height=5)
dittoBarPlot(kor_clean_comb, "celltype_new", group.by = "nearestCMS.RF",color.panel = ,main=NULL)
dev.off()

svg("kor_mast.svg", width=7, height=5)
dittoDimPlot(kor, "KIT",reduction.use="tsne",assay="RNA",max.color = "#FF0000", min.color = "lightgrey",do.label = FALSE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = FALSE, labels.size = 3, show.axes.numbers = FALSE, xlab=NULL,ylab=NULL,main=NULL)+ labs(title = "MAST Cells: KIT")
dev.off()
#######


svg("yuan_immune_Barplots.svg", width=7, height=7)
dittoBarPlot(yuan.immune, "seurat_clusters", group.by = "orig.ident",color.panel = c("#d1001c", "#0881d1", "#014600"),main=NULL)
dev.off()


svg("yuan_stroma_Barplots.svg", width=7, height=7)
dittoBarPlot(yuan.stromal, "seurat_clusters", group.by = "orig.ident",color.panel = c("#ff6600","#9900ff"),main=NULL)
dev.off()



svg("yuan_conditions.svg", width=7, height=5)
dittoBarPlot(yuan_cr, "seurat_clusters", group.by = "condition_new",color.panel = c("#d1001c", "#ff6600", "#d5b60a", "#9900ff","#0881d1","#014600"),main=NULL)
dev.off()


#Epi

svg("tsne_epi_clusters.svg", width=7, height=7, res=300)
dittoDimPlot(epi, "seurat_clusters", reduction.use ="tsne",color.panel = c("#d1001c", "#9900ff", "#ff6600", "#d5b60a", "#0881d1","#014600","#ff8000","#40ff00","#0040ff","  darkgreen","#3399ff","#9900ff","#993399","#660033","#993300","#e6e600","#999966"),do.label = TRUE, labels.repel = FALSE,labels.highlight = TRUE,legend.show = FALSE, labels.size = 3, show.axes.numbers = FALSE, xlab=NULL,ylab=NULL,main=NULL)
dev.off()


## New Figures for menuscript 10/11/2021

svg("immune_stromal_total.svg", width=7, height=7)
cowplot::plot_grid(ncol = 3,
dittoBarPlot(alldata_new, "seurat_clusters_new", group.by = "orig.ident",color.panel = c("#d1001c", "#9900ff", "#ff6600", "#d5b60a", "#0881d1","#014600"),main=NULL),
dittoBarPlot(immune.comp, "seurat_clusters_new", group.by = "orig.ident",color.panel = c("#d1001c", "#9900ff", "#ff6600", "#d5b60a", "#0881d1","#014600"),main=NULL),
dittoBarPlot(stromal.comp, "seurat_clusters_new", group.by = "orig.ident",color.panel = c("#d1001c", "#9900ff", "#ff6600", "#d5b60a", "#0881d1","#014600"),main=NULL))
dev.off()


svg("immune_comp.svg", width=10, height=10)
dittoBarPlot(immune.comp, "seurat_clusters_new", group.by = "orig.ident",color.panel = c("#d1001c", "#9900ff", "#ff6600", "#d5b60a", "#0881d1","#014600"),main=NULL)
dev.off()


svg("total_comp.svg", width=10, height=10)
dittoBarPlot(alldata_new, "seurat_clusters_new", group.by = "orig.ident",color.panel = c("#d1001c", "#9900ff", "#ff6600", "#d5b60a", "#0881d1","#014600"),main=NULL)
dev.off()

svg("stromal_comp.svg", width=10, height=10)
dittoBarPlot(stromal.comp, "seurat_clusters_new", group.by = "orig.ident",color.panel = c("#d1001c", "#9900ff", "#ff6600", "#d5b60a", "#0881d1","#014600"),main=NULL)
dev.off()

svg("cms_total.svg", width=10, height=10)
dittoBarPlot(alldata_new, "seurat_clusters_new", group.by = "bulk_prediction",color.panel = c("#d1001c", "#9900ff", "#ff6600", "#d5b60a","#0881d1","#014600"),main=NULL)
dev.off()

svg("normal_tumor_total.svg", width=10, height=10)
dittoBarPlot(alldata_new, "seurat_clusters_new", group.by = "Condition",color.panel = c("#d1001c", "#9900ff", "#ff6600", "#d5b60a","#0881d1","#014600"),main=NULL)
dev.off()


# Fibro


genes <- c("CXCL13","LUM","DCN","PDLIM1","VCAN","AEBP1","FBLN1","MFAP4","COL4A1","COL18A1","SPARCL1","ITGA1","CD248","PDGFRB","EPAS1","NOTCH3","MCAM","RBPMS","CREG1","ANGPTL2","NR2F2","COL5A3","NID2","MEF2C","ESAM","PDGFA","KCNJ8")

svg("fibro_genes_exp.svg", width=10, height=10)
dittoDotPlot(fibro, vars = genes, group.by = "seurat_clusters",min.color = "#2eafff",max.col="#ff4d00")
dev.off()


svg("tsne_fibro_clusters.svg", width=4, height=4)
dittoDimPlot(fibro, "seurat_clusters", reduction.use ="tsne",color.panel = c("#d1001c","#338000","#9900ff"),do.label = TRUE, labels.repel = FALSE,labels.highlight = TRUE,legend.show = FALSE, labels.size = 3, show.axes.numbers = FALSE, xlab=NULL,ylab=NULL,main=NULL)
dev.off()

svg("tsne_int_bc_crc.svg", width=10, height=15)
FeaturePlot(int_bc_crc, features = c("GJB2", "SCARA5", "ADH1B", "SEMA3C","TGFB1"), reduction="tsne",split.by = "disease_type", max.cutoff = 3, cols = c("grey", "red"),label.size = 4)
dev.off()

("ANTXR1", "SDC1", "LAMP5", "CD9", "GPC3","DLK1")


png("tsne_int_bc_crc_2.png", units="in",width=10, height=15,res=300)
FeaturePlot(int_bc_crc, features = c("ANTXR1", "SDC1", "LAMP5", "CD9", "GPC3","DLK1"), reduction="tsne",split.by = "disease_type", max.cutoff = 3, cols = c("grey", "red"),label.size = 4)
dev.off()



png("tsne_int_bc_crc.png", units="in",width=10, height=15,res=300)
FeaturePlot(int_bc_crc, features = c("GJB2", "SCARA5", "ADH1B", "SEMA3C","TGFB1"), reduction="tsne",split.by = "disease_type", max.cutoff = 3, cols = c("grey", "red"),label.size = 4)
dev.off()

png("tsne_int_lee_Ourcrc.png", units="in",width=10, height=15,res=300)
FeaturePlot(immune.combined_0.5, features = c("GJB2", "SCARA5", "ADH1B", "SEMA3C","CST1","TGFB1"), reduction="tsne",split.by = "data_type", max.cutoff = 5, cols = c("darkgrey", "red"),label.size = 4,pt.size=1)
dev.off()

svg("tsne_int_lee_Ourcrc.svg", width=10, height=15)
FeaturePlot(immune.combined_0.5, features = c("GJB2", "SCARA5", "ADH1B", "SEMA3C","CST1","TGFB1"), reduction="tsne",split.by = "data_type", max.cutoff = 5, cols = c("darkgrey", "red"),label.size = 4,pt.size=1)
dev.off()

genes<-c("FAP","PDPN","PDGFRA","CXCL14","MFAP5","FBLN1","MFAP4","LOXL1","DCN","LUM","CREG1","SPARCL1","RGS5","CSPG4","COL5A3","NID2","COL4A1","COL4A2","KCNJ8","RBPMS","MCAM","NOTCH3","EPAS1","COL18A1","NR2F2","PDGFA","ESAM","ITGA1","PDGFRB")

CAFS1<-c("PDPN","CXCL14","MFAP5","FBLN1","MFAP4","LOXL1","DCN","LUM","CREG1","SPARCL1")



svg("fibro_genes_exp.svg", width=10, height=10)
dittoDotPlot(fibro, vars = genes, group.by = "seurat_clusters",min.color = "#2eafff",max.col="#ff4d00")
dev.off()




svg("tsne_fibro_markers.svg", width=7, height=7, res=300)
dittoDimPlot(fibro, "seurat_clusters", reduction.use ="tsne",color.panel = c("#d1001c", "#9900ff", "#ff6600", "#d5b60a", "#0881d1","#014600","#ff8000","#40ff00","#0040ff","  darkgreen","#3399ff","#9900ff","#993399","#660033","#993300","#e6e600","#999966"),do.label = TRUE, labels.repel = FALSE,labels.highlight = TRUE,legend.show = FALSE, labels.size = 3, show.axes.numbers = FALSE, xlab=NULL,ylab=NULL,main=NULL)
dev.off()


gene1 <- c("FAP", "PDPN", "PDGFRA","COL1A1","COL1A2","CXCL12", "CXCL14","RGS5", "CSPG4","PDGFRB","CD248","EPAS1")
dittoHeatmap(fibro, genes, annot.by = "seurat_clusters", highlight.features = genes[1:3], complex = TRUE,scaled.to.max = TRUE, heatmap.colors = colorRampPalette(c("blue", "white", "yellow"))(50))



d1=DoHeatmap(object = fibro,genes) + scale_fill_viridis(option="magma")
ggsave(file="heatmap_fibro_markers.svg", plot=d1, width=10, height=8)
dev.off()

png("heatmap_fibro_markers.png", units="in",width=10, height=15,res=300)
DoHeatmap(object = fibro,genes) + scale_fill_viridis(option="magma")
 dev.off()



 ### fibro-caf-s1 subset and reclustering

 obj<- SubsetData(fibro, subset.name = "seurat_clusters", accept.value = "CAF-S1")



obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, verbose = FALSE, features = VariableFeatures(object = obj))
obj <- RunTSNE(object = obj, dims.use = 1:30)
obj <- RunUMAP(obj, dims = 1:30)
obj <- FindNeighbors(obj, dims = 1:30, force.recalc = T)

obj<- FindClusters(obj, resolution = 0.6, algorithm = 1, force.recalc=T)
obj<- FindClusters(obj, resolution = 1.5, algorithm = 1, force.recalc=T)
obj$seurat_clusters <- as.factor(as.numeric(as.character(obj$seurat_clusters)) + 1)
Idents(obj) <- "seurat_clusters"



table(obj$seurat_clusters)
png("tsne_obj.png", units="in", width=10, height=10, res=300)
DimPlot(obj, reduction = "tsne", label=T, group.by="seurat_clusters")
dev.off()

markers_genes_obj <- FindAllMarkers(obj, logfc.threshold = 0.2, test.use = "wilcox", only.pos = TRUE,assay = "RNA")
write.table(markers_genes_obj, file="marker_genes_obj.1.txt", sep="\t",append=F, row.names = F)
cbind(table(markers_genes_obj$cluster))

top25 <- markers_genes_obj  %>% group_by(cluster) %>% top_n(-25, p_val_adj)

obj <- ScaleData(obj, features = as.character(unique(top5$gene)), assay = "RNA")
DoHeatmap(obj, features = as.character(unique(top10$gene)), group.by = "seurat_clusters", assay = "RNA")

multi_dittoPlot(obj, mycaf.genes[1:9], group.by = "seurat_clusters",vlnplot.lineweight = 0.2, jitter.size = 0.3)

mycaf.genes <- c("TAGLN", "ACTA2", "MMP11", "PDGFRB", "HOPX", "POSTN","MYL9","TPM1","TPM2","FAP","PDPN","COL12A1")
icaf.genes <- c("IL6", "PDGFRA", "CXCL12", "CFD", "DPT", "LMNA", "AGTR1", "HAS1", "CXCL1", "CXCL2", "CCL2", "IL8","PLA2G2A","CLU","EMP1")
 

current.cluster.ids <- c(1,2,3,4,5,6,7,8,9)
new.cluster.ids <- c("iCAF","myCAF","myCAF","myCAF","myCAF","iCAF","iCAF","iCAF","myCAF")

obj$seurat_clusters <- plyr::mapvalues(x = obj$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)

kable(table(obj@meta.data$seurat_clusters))


png("caf_subtype.png", width=5, height=5, res=300)
#svg("caf_subtype.svg", width=5, height=5)
DimPlot(obj, reduction = "tsne", label=T, group.by="seurat_clusters",pt.size=4,cols=c("#008000","#D40000"))
dev.off()



FeaturePlot(obj, reduction = "tsne",features = mycaf.genes ,min.cutoff="q9", cols=c("lightgrey", "red"), label=TRUE) + labs(title = "myCAF")


caf.genes<-c("CXCL12","CFD","DPT","CCL2","TAGLN","ACTA2","MMP11","PDGFRB","FAP","PDPN","COL12A1")


svg("cafs_markers.svg", width=7, height=15)
StackedVlnPlot(obj = obj, features = caf.genes)
dev.off()

svg("meyloid_markers_mycaf.svg", width=15, height=10)
png("meyloid_markers_mycaf.png", units="in",width=8, height=10,res=300)
multi_dittoDimPlot(obj, mycaf.genes,reduction.use ="tsne",axes.labels.show=FALSE,legend.show = TRUE,min.color="lightgrey", max.color="red")
dev.off()
 

svg("meyloid_markers_icaf.svg", width=15, height=10)
png("meyloid_markers_icaf.png", units="in",width=8, height=10,res=300) 
multi_dittoDimPlot(obj, icaf.genes,reduction.use ="tsne",axes.labels.show=FALSE,legend.show = TRUE,min.color="lightgrey", max.color="red")
dev.off()


DotPlot(object = obj, features = c("CXCL12","PDGFD","CD34","FGF10","CXCL1","ID2","DPT","CFD","CLU","MCAM","PDGFRA","SOD2","TAGLN","COL1A2","ANTXR1","GPC3","LAMP5","CD9","SDC1","TGFB1"))+theme(axis.text.x = element_text(angle = 90))

 ####

ggsave(file="test.svg", plot=d1, width=10, height=8)



d2=DoHeatmap(object = fibro,genes)
ggsave(file="test.svg", plot=d2, width=10, height=8)


my_genes <- c("CD1C","C1QA","APOE", "SEPP1", "CD163","TAM","MRC1","S100A8","1L1B","S100A9","FCGR3B")


svg("myeloid_heatmap.svg", width=10, height=10)
DoHeatmap(object = myeloid,my_genes,group.by="seurat_clusters_new") + scale_fill_viridis(option="magma")
dev.off()

library(viridis)
ggsave(file="test1.svg", plot=DoHeatmap(object = myeloid,my_genes,group.by="seurat_clusters_new") + scale_fill_viridis(option="magma"), width=10, height=10)


svg("tsne_myeloid_markers.svg", width=7, height=7, res=300)
dittoDimPlot(myeloid, "seurat_clusters_new", reduction.use ="tsne",color.panel = c("#d1001c", "#9900ff", "#ff6600", "#d5b60a", "#0881d1","#014600","#ff8000","#40ff00","#0040ff","  darkgreen","#3399ff","#9900ff","#993399","#660033","#993300","#e6e600","#999966"),do.label = TRUE, labels.repel = FALSE,labels.highlight = TRUE,legend.show = FALSE, labels.size = 3, show.axes.numbers = FALSE, xlab=NULL,ylab=NULL,main=NULL)
dev.off()




myeloid <- NormalizeData(myeloid, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE,force.recalc = T)
myeloid <- FindVariableFeatures(myeloid, selection.method = "vst", nfeatures = 2000, verbose = FALSE,force.recalc = T)
myeloid <- ScaleData(myeloid, verbose = FALSE,force.recalc = T)
myeloid <- RunPCA(myeloid, verbose = FALSE, features = VariableFeatures(object = myeloid),force.recalc = T)
myeloid <- RunTSNE(object = myeloid, dims.use = 1:30,force.recalc = T)
myeloid <- RunUMAP(myeloid, dims = 1:30,force.recalc = T)

myeloid <- FindNeighbors(myeloid, dims = 1:30, force.recalc = T)
myeloid <- FindClusters(myeloid, force.recalc = T)



dittoHeatmap(myeloid, my_genes, annot.by = "seurat_clusters_new", scaled.to.max = TRUE, heatmap.colors = colorRampPalette(c("blue", "white", "yellow"))(50))

### Myeloid markers

Macrophages --> CD163 [phenotypic moulding], PTGS2[cancer cell]
Granulocytes --> S100A12 [phenotypic moulding]
Dendritic cells --> CLEC9A, FCER1A, XCR1, CDC1, LAMP3 [phenotypic moulding]
Monocytes --> LYBC2, CHIL3, CYBB [Cancer cell], S100A8, S100A9, FGER3B, CD163
TAM-C1Q+  --> C1QA, C1QB,C1QC,CD63,PLTP , HLA-DQB2, CD74, HLA-DQB1 [Cancer cell]

mygenes <- c("CD163", "S100A12", "CLEC9A", "FCER1A", "XCR1", "CDC1", "LYBC2", "CHIL3", "CYBB", "C1QA", "C1QB","C1QC", "CD63", "PLTP", "HLADQB2", "CD74", "HLADQB1" )
mygenes <- c("S100A12", "CLEC9A","FCER1A","LAMP3","S100A8","S100A9","FCGR3B","CD163","C1QA","C1QB","C1QC","CD63","PLTP","CD74","APOE","MRC1","MARCO" )
DoHeatmap(object = yuan_myeloid_cr,mygenes,group.by="seurat_clusters") + scale_fill_viridis(option="magma")

#dittoHeatmap(kor_myeloid, my_genes, annot.by = "seurat_clusters", scaled.to.max = TRUE, heatmap.colors = colorRampPalette(c("blue", "white", "yellow"))(50))



mygenes <- c("S100A12", "CLEC9A","FCER1A","LAMP3","S100A8","S100A9","FCGR3B","CD163","C1QA","C1QB","C1QC","PLTP","APOE","MRC1","MARCO" )

svg("myloid_markers.svg", width=7, height=7)
multi_dittoDimPlot(myeloid, mygenes,reduction.use ="tsne",axes.labels.show=FALSE,legend.show = TRUE)
dev.off()


c("CXCL12","PDGFD","CD34","FGF10","CXCL1","ID2","DPT","CFD","CLU","MCAM","PDGFRA","SOD2","TAGLN","COL1A2","ANTXR1","GPC3","LAMP5","CD9","SDC1","TGFB1")

dittoDotPlot(object = fibro, features = c("CXCL12","PDGFD","CD34","FGF10","CXCL1","ID2","DPT","CFD","CLU","MCAM","PDGFRA","SOD2","TAGLN","COL1A2","ANTXR1","GPC3","LAMP5","CD9","SDC1","TGFB1"),min.color = "#2eafff",max.col="#ff4d00")

dittoDotPlot(fibro, vars = genes, group.by = "seurat_clusters",min.color = "#2eafff",max.col="#ff4d00")



####
dittoHeatmap(kor_myeloid, my_genes, annot.by = "seurat_clusters", scaled.to.max = TRUE, heatmap.colors = colorRampPalette(c("blue", "white", "yellow"))(50))


DotPlot(object = myeloid, features = my_genes)+ scale_colour_gradient2(low="steelblue", mid="lightgrey", high="darkgoldenrod1")


png("heatmap_yuan_myeloid_markers.png", units="in",width=15, height=10,res=300)
DoHeatmap(yuan_myeloid,my_genes) + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + guides(color=FALSE)
dev.off()

png("heatmap_korean_myeloid_markers.png", units="in",width=15, height=10,res=300)
DoHeatmap(object = kor_myeloid,my_genes) + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + guides(color=FALSE)
dev.off()


ggsave(file="test3.svg", plot=DoHeatmap(myeloid,my_genes) + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + guides(color=FALSE))


# markers fibro

svg("fibro_caf_markers.svg", width=7, height=7)
multi_dittoDimPlot(fibro, c("CXCL12","CFD","FAP","PDPN","RGS5","MCAM"),reduction.use ="tsne",axes.labels.show=TRUE,legend.show = TRUE)
dev.off()



#Heatmap

library(Seurat)
library(dplyr)
library(ComplexHeatmap)
library(circlize) # for the colorRamp2() function

#seurat_object <- readRDS(...)

# Annotations

meta_data <- myeloid@meta.data

# order of annotations/colors are defined here
ordered_meta_data <- meta_data[order(meta_data$seurat_clusters_new), ]

# OPTIONAL: YOU CAN PICK COLORS FOR EACH LEVEL OF ANNOTATION
# HERE I PROVIDE COLORS FOR ONLY TWO OF THE THREE LEVELS, COLORS FOR THE
# REMAINING LEVEL IS TAKEN CARE OF BY THE PACKAGE
annotation_colors <- list("seurat_clusters_new"= c("CD1C-DENDRITIC CELLS" = "#d1001c",
                                         "GRANULOCYTES" = "#9900ff",
                                         "MONOCYTES" = "#ff6600",
                                         "TAM" = "#d5b60a"),
                          "bulk_prediction" = c("CMS1" = "blue",
                                              "CMS2" = "yellow",
                                              "CMS3" = "maroon",
                                              "CMS4" = "darkgreen"))
ha = HeatmapAnnotation(df = ordered_meta_data,
                       show_annotation_name = TRUE,
                       col = annotation_colors)

# Expression data

genes_to_use <- c("CD1C","C1QA","APOE", "SEPP1", "CD163","TAM","MRC1","S100A8","1L1B","S100A9","FCGR3B")

seurat_object <- ScaleData(myeloid,
                           genes.use = genes_to_use)

my_data <- seurat_object@assays$RNA@scale.data

# COLUMN ORDER OF THE EXPRESSION DATA SHOULD MATCH THE ROW ORDER OF THE
# ANNOTATION TABLE
my_data <- my_data[, rownames(ordered_meta_data)]

# Heatmap

col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
png(filename = "example_heatmap.png",
    width = 1000,
    height = 1000)
Heatmap(
  my_data,
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_order = NULL,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  use_raster = TRUE,
  raster_device = c("png"),
  bottom_annotation = NULL,
  top_annotation = ha
)
dev.off()

#myeloid


svg("tsne_myeloid_dim.svg", width=7, height=7)
dittoDimPlot(myeloid, "seurat_clusters_new", reduction.use ="tsne",color.panel = c("#d1001c", "#9900ff", "#ff6600", "#d5b60a", "#0881d1","#014600","#ff8000","#40ff00","#0040ff","  darkgreen","#3399ff","#9900ff","#993399","#660033","#993300","#e6e600","#999966"),do.label = TRUE, labels.repel = FALSE,labels.highlight = TRUE,legend.show = FALSE, labels.size = 3, show.axes.numbers = FALSE, xlab=NULL,ylab=NULL,main=NULL,size=2)
dev.off()

ggsave(file="test3.png", plot=DoHeatmap(myeloid,my_genes, group.by="seurat_clusters_new") + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + guides(color=FALSE))



ggsave(file="test4.png", plot=DoHeatmap(myeloid,my_genes, group.by="bulk_prediction") + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + guides(color=FALSE))



svg("tsne_myeloid_dim.svg", width=7, height=7)
dittoDimPlot(myeloid, "seurat_clusters_new", reduction.use ="tsne",color.panel = c("#d1001c", "#1100ff", "#06540e", "#d5b60a", "#0881d1","#014600","#ff8000","#40ff00","#0040ff","  darkgreen","#3399ff","#9900ff","#993399","#660033","#993300","#e6e600","#999966"),do.label = TRUE, labels.repel = FALSE,labels.highlight = TRUE,legend.show = FALSE, labels.size = 3, show.axes.numbers = FALSE, xlab=NULL,ylab=NULL,main=NULL,size=5)
dev.off()

#  Tcells Rename




current.cluster.ids <- c("CD4-ANXA1","CD4-CLUSTER-6","CD4-CXCL13","CD4-HSP","CD4-MEMORY","CD4-TH17","CD4-TREG","CD8-EFFECTOR 1","CD8-EFFECTOR 2","CD8-EFFECTOR 3","CD8-MAIT CELLS","CD8-MEMORY 1","CD8-MEMORY 2","CD8-TRM","INNATE LYMPHOID CELLS","NK-CELLS")
new.cluster.ids <- c("CD4+ Resting cells","CD4+ CLUSTER 6","CD4+ CXCL13","CD4+ HSP","CD4+ Memory","CD4+ TH17","CD4+ Tregs","CD8+ GZMK","CD8+ GZMA","CD8+ GZMB","CD8+ MAIT Cells","CD8+ Memory 1","CD8+ Memory 2","CD8+ TRM","ILC","NK cells")
tcells$seurat_clusters <- plyr::mapvalues(x = tcells$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)

kable(table(tcells_new@meta.data$seurat_clusters))


#####


#Tcells



svg("tsne_myeloid_dim.svg", width=7, height=7)
dittoDimPlot(tcells, "seurat_clusters", reduction.use ="tsne",color.panel = c("#d1001c", "#9900ff", "#ff6600", "#d5b60a", "#0881d1","#014600","#ff8000","#40ff00","#0040ff","  darkgreen","#3399ff","#9900ff","#993399","#660033","#993300","#e6e600","#999966"),do.label = TRUE, labels.repel = FALSE,labels.highlight = TRUE,legend.show = FALSE, labels.size = 3, show.axes.numbers = FALSE, xlab=NULL,ylab=NULL,main=NULL,size=2)
dev.off()

svg("tsne_tcells_dim.svg", width=7, height=7)
dittoDimPlot(tcells, "seurat_clusters", reduction.use ="tsne",color.panel = c("#d1001c", "#9900ff", "#ff6600", "#d5b60a", "#0881d1","#014600","#ff8000","#40ff00","#0040ff","  darkgreen","#3399ff","#9900ff","#993399","#660033","#993300","#e6e600","#999966"),do.label = FALSE, labels.repel = FALSE,labels.highlight = FALSE,legend.show = TRUE, labels.size = 3, show.axes.numbers = FALSE, xlab=NULL,ylab=NULL,main=NULL,size=1)
dev.off()



tcells_marker <- c("FOXP3", "CTLA4", "PD1", "PDL1","LAG3","CCR7", "SELL","TCF7","ANXA1", "IL7R","LMNA","CXCL13","GZMA","GZMB","GZMH","IFNG","PRF1","HAVCR2","CD96","EOMES","KLRG1")


svg("tcells_marker.svg", width=10, height=10)
DoHeatmap(tcells,tcells_marker,group.by="seurat_clusters") + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + guides(color=FALSE)
dev.off()

ggsave(file="test4.svg", plot=DoHeatmap(tcells,tcells_marker,group.by="seurat_clusters") + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + guides(color=FALSE))

png("tcells_marker.png",units="in", width=10, height=10, res=300)
DoHeatmap(tcells,tcells_marker,group.by="seurat_clusters") + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) + guides(color=FALSE)
dev.off()



multi_dittoDimPlot(tcells,tcells_marker[1:18])

FeaturePlot(tcells,tcells_marker)

svg("tcells_marker.svg", width=10, height=10)

svg("tcells_marker_tsne.svg", width=10, height=10)
FeaturePlot(tcells, reduction = "tsne",features = tcells_marker ,min.cutoff="q9", cols=c("lightgrey", "red"))
dev.off()


png("tcells_marker_tsne.png",units="in", width=10, height=10, res=300)
FeaturePlot(tcells, reduction = "tsne",features = tcells_marker ,min.cutoff="q9", cols=c("lightgrey", "darkred"))
dev.off()



svg("tcells_bar_samples.svg", width=10, height=10)
dittoBarPlot(tcells, "seurat_clusters", group.by = "orig.ident",color.panel = c("#d1001c", "#9900ff", "#ff6600", "#d5b60a", "#0881d1","#014600","#ff8000","#40ff00","#0040ff","  darkgreen","#3399ff","#9900ff","#993399","#660033","#993300","#e6e600","#999966"),main=NULL)
dev.off()
#B-cells


svg("tsne_bcells_dim.svg", width=7, height=7)
dittoDimPlot(bcells, "seurat_clusters", reduction.use ="umap",color.panel = c("#d1001c", "#9900ff", "#00ff77", "#d5b60a", "#0881d1","#014600","#ff8000","#40ff00","#0040ff","darkgreen","#3399ff","#9900ff","#993399","#660033","#993300","#e6e600","#999966"),do.label = TRUE, labels.repel = FALSE,labels.highlight = TRUE,legend.show = FALSE, labels.size = 3, show.axes.numbers = FALSE, xlab=NULL,ylab=NULL,main=NULL,size=2)
dev.off()

/Users/akhaliq/Desktop/Figures/Final Folder Cancer Cell/References/My Library.bib
#Lee

# for traj epi lee --> load("/Users/akhaliq/Desktop/All_Clusters/Trajectory_analysis_Epi_Korean_02022021.RData")

svg("tsne_yuan_dim.svg", width=7, height=7)
dittoDimPlot(yuan_epi, "seurat_clusters", reduction.use ="umap",color.panel = c("#d1001c", "#9900ff", "#00ff77", "#d5b60a", "#0881d1","#014600","#ff8000","#40ff00","#0040ff","  darkgreen","#3399ff","#9900ff","#993399","#660033","#993300","#e6e600","#999966"),do.label = TRUE, labels.repel = FALSE,labels.highlight = TRUE,legend.show = FALSE, labels.size = 3, show.axes.numbers = FALSE, xlab=NULL,ylab=NULL,main=NULL,size=2)
dev.off()


svg("trajectory_yuan.svg", width=10, height=7)
cowplot::plot_grid(ncol = 3,
plot_cell_trajectory(cds1, color_by="State"),
plot_cell_trajectory(cds1, color_by="sc_cms"),
plot_cell_trajectory(cds1, color_by="nearestCMS.RF"),
plot_cell_trajectory(cds1, color_by="Condition"),
plot_cell_trajectory(cds1, color_by="Location"),
plot_cell_trajectory(cds1, color_by="Pseudotime"))
dev.off()


svg("tsne_yuan_dim.svg", width=7, height=7)
cowplot::plot_grid(ncol = 1,
dittoDimPlot(yuan_epi, "seurat_clusters", reduction.use ="umap",color.panel = c("#d1001c", "#9900ff", "#00ff77", "#d5b60a", "#0881d1","#014600","#ff8000","#40ff00","#0040ff","  darkgreen","#3399ff","#9900ff","#993399","#660033","#993300","#e6e600","#999966"),do.label = TRUE, labels.repel = FALSE,labels.highlight = TRUE,legend.show = FALSE, labels.size = 3, show.axes.numbers = FALSE, xlab=NULL,ylab=NULL,main=NULL,size=2),
dittoDimPlot(yuan_epi, "Condition", reduction.use ="umap",color.panel = c("#d1001c", "#9900ff", "#00ff77", "#d5b60a", "#0881d1","#014600","#ff8000","#40ff00","#0040ff","  darkgreen","#3399ff","#9900ff","#993399","#660033","#993300","#e6e600","#999966"),do.label = TRUE, labels.repel = FALSE,labels.highlight = TRUE,legend.show = FALSE, labels.size = 3, show.axes.numbers = FALSE, xlab=NULL,ylab=NULL,main=NULL,size=2),
dittoDimPlot(yuan_epi, "MSI_Status", reduction.use ="umap",color.panel = c("#d1001c", "#9900ff", "#00ff77", "#d5b60a", "#0881d1","#014600","#ff8000","#40ff00","#0040ff","  darkgreen","#3399ff","#9900ff","#993399","#660033","#993300","#e6e600","#999966"),do.label = TRUE, labels.repel = FALSE,labels.highlight = TRUE,legend.show = FALSE, labels.size = 3, show.axes.numbers = FALSE, xlab=NULL,ylab=NULL,main=NULL,size=2))
dev.off()

## Korean re 

png("korean_dim.png",units="in", width=15, height=10, res=300)
#svg("korean_dim.svg", width=15, height=7)
cowplot::plot_grid(ncol = 2, 
dittoDimPlot(alldata, "Cell_type",do.label = TRUE, labels.repel = FALSE,reduction="tsne")+ labs(title = "After QC"))
dev.off()


svg("korean_myloid_conta_epi.svg", width=15, height=7)
cowplot::plot_grid(ncol = 2,
FeaturePlot(koreandata,features = "LYZ",reduction="tsne" ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE)+ labs(title = "Myloid Contamination in Epithelial cell clusters"),
FeaturePlot(alldata,features = "LYZ",reduction="tsne" ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE))
dev.off()



png("Epi_Genes.png", units="in", width=15, height=6, res=300)
cowplot::plot_grid(ncol = 2,
FeaturePlot(epi.tum,features = c("IL6","OSMR") ,reduction="umap" ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE),
FeaturePlot(epi.tum,features = c("IL6","OSMR") ,reduction="tsne" ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE))
dev.off()


png("korean_Data.png", units="in", width=40, height=10, res=300)
cowplot::plot_grid(ncol = 6,
FeaturePlot(koreandata,features = "KRT18" ,reduction="tsne" ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE)+ labs(title = "Epithelial Cells"),
FeaturePlot(koreandata,features = "CD3D" ,reduction="tsne" ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE)+ labs(title = "T Cells"),
FeaturePlot(koreandata,features = "CD79A" ,reduction="tsne" ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE)+ labs(title = "B Cells"),
FeaturePlot(koreandata,features = "LYZ" ,reduction="tsne" ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE)+ labs(title = "Myeloid"),
FeaturePlot(koreandata,features = "COL1A1",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE)+ labs(title = "Fibroblast"),
FeaturePlot(koreandata,features = "CLDN5" ,reduction="tsne" ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE)+ labs(title = "Endothelial"))
dev.off()

# Epi

svg("korean_epicontamination_1.svg", width=15, height=7)
cowplot::plot_grid(ncol = 6,
FeaturePlot(kor.epi,features = "JCHAIN",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=4),
FeaturePlot(kor.epi,features = "IGHG1",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=4),
FeaturePlot(kor.epi,features = "IGHG3",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=4),
FeaturePlot(kor.epi,features = "CD3D",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=4),
FeaturePlot(kor.epi,features = "COL1A1",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=4),
FeaturePlot(kor.epi,features = "CD68",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=4))
dev.off()


svg("korean_Cleaned_epi_contamination_2.svg", width=15, height=7)
cowplot::plot_grid(ncol = 5,
FeaturePlot(kor.epi,features = "JCHAIN",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=4),
FeaturePlot(kor.epi,features = "IGHG1",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=4),
FeaturePlot(kor.epi,features = "IGHG3",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=4),
FeaturePlot(kor.epi,features = "CD3D",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=4),
FeaturePlot(kor.epi,features = "COL1A1",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=4))
dev.off()

svg("korean_Cleaned_epi_contamination_3.svg", width=15, height=7)
cowplot::plot_grid(ncol = 6,
FeaturePlot(kor.epi,features = "JCHAIN",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=4),
FeaturePlot(kor.epi,features = "IGHG1",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=4),
FeaturePlot(kor.epi,features = "IGHG3",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=4),
FeaturePlot(kor.epi,features = "CD3D",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=4),
FeaturePlot(kor.epi,features = "COL1A1",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=4),
FeaturePlot(kor.epi,features = "CD68",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=4))
dev.off()


svg("our_Cleaned_epi_contamination_4.svg", width=15, height=7)
cowplot::plot_grid(ncol = 6,
FeaturePlot(epi,features = "JCHAIN",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=4),
FeaturePlot(epi,features = "IGHG1",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=4),
FeaturePlot(epi,features = "IGHG3",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=4),
FeaturePlot(epi,features = "CD3D",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=4),
FeaturePlot(epi,features = "COL1A1",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=4),
FeaturePlot(epi,features = "CD68",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=4))
dev.off()

#Myeloid

svg("korean_myeloid_contamination_1.svg", width=15, height=7)
cowplot::plot_grid(ncol = 5,
FeaturePlot(kor.myloid,features = "JCHAIN",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5),
FeaturePlot(kor.myloid,features = "IGHA1",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5),
FeaturePlot(kor.myloid,features = "IGHG3",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5),
FeaturePlot(kor.myloid,features = "KRT8",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5),
FeaturePlot(kor.myloid,features = "CD3D",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5))
dev.off()

svg("korean_myeloid_contamination_2.svg", width=15, height=7)
cowplot::plot_grid(ncol = 5,
FeaturePlot(kor.myloid,features = "JCHAIN",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5),
FeaturePlot(kor.myloid,features = "IGHA1",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5),
FeaturePlot(kor.myloid,features = "IGHG3",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5),
FeaturePlot(kor.myloid,features = "KRT8",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5),
FeaturePlot(kor.myloid,features = "CD3D",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5))
dev.off()


svg("korean_myeloid_cleaned_3.svg", width=15, height=7)
cowplot::plot_grid(ncol = 5,
FeaturePlot(kor.myeloid,features = "JCHAIN",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5),
FeaturePlot(kor.myeloid,features = "IGHA1",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5),
FeaturePlot(kor.myeloid,features = "IGHG3",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5),
FeaturePlot(kor.myeloid,features = "KRT8",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5),
FeaturePlot(kor.myeloid,features = "CD3D",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5))
dev.off()



svg("our_myeloid_cleaned_4.svg", width=15, height=5)
cowplot::plot_grid(ncol = 5,
FeaturePlot(myeloid,features = "JCHAIN",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=3),
FeaturePlot(myeloid,features = "IGHA1",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=3),
FeaturePlot(myeloid,features = "IGHG3",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=3),
FeaturePlot(myeloid,features = "KRT8",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=3),
FeaturePlot(myeloid,features = "CD3D",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=3))
dev.off()


# Fibro

svg("korean_fibro_contamination_1.svg", width=15, height=7)
cowplot::plot_grid(ncol = 6,
FeaturePlot(kor.stroma,features = "JCHAIN",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=2), 
FeaturePlot(kor.stroma,features = "IGHA1",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=2),
FeaturePlot(kor.stroma,features = "IGHG3",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=2),
FeaturePlot(kor.stroma,features = "KRT8",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=2),
FeaturePlot(kor.stroma,features = "CD3D",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=2),
FeaturePlot(kor.stroma,features = "CD68",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=2))
dev.off()

svg("korean_fibro_contamination_2.svg", width=15, height=7)
cowplot::plot_grid(ncol = 6,
FeaturePlot(kor.stroma,features = "JCHAIN",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=4), 
FeaturePlot(kor.stroma,features = "IGHA1",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=4),
FeaturePlot(kor.stroma,features = "IGHG3",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=4),
FeaturePlot(kor.stroma,features = "KRT8",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=4),
FeaturePlot(kor.stroma,features = "CD3D",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=4),
FeaturePlot(kor.stroma,features = "CD68",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=4))
dev.off()

png("korean_Data.png", units="in", width=20, height=10, res=300)
#svg("korean_fibro_cleaned_3.svg", width=15, height=4)
cowplot::plot_grid(ncol = 6,
FeaturePlot(kor.fibro,features = "JCHAIN",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=2), 
FeaturePlot(kor.fibro,features = "IGHA1",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=2),
FeaturePlot(kor.fibro,features = "IGHG3",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=2),
FeaturePlot(kor.fibro,features = "KRT8",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=2),
FeaturePlot(kor.fibro,features = "CD3D",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=2),
FeaturePlot(kor.fibro,features = "CD68",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=2))
dev.off()

svg("our_fibro_cleaned_4.svg", width=15, height=4)
cowplot::plot_grid(ncol = 6,
FeaturePlot(fibro,features = "JCHAIN",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=4), 
FeaturePlot(fibro,features = "IGHA1",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=4),
FeaturePlot(fibro,features = "IGHG3",reduction="tsne"  , cols=c("lightgrey", "red"), label=FALSE,pt.size=4),
FeaturePlot(fibro,features = "KRT8",reduction="tsne"  , cols=c("lightgrey", "red"), label=FALSE,pt.size=4),
FeaturePlot(fibro,features = "CD3D",reduction="tsne"  , cols=c("lightgrey", "red"), label=FALSE,pt.size=4),
FeaturePlot(fibro,features = "CD68",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=4))
dev.off()

#Tcells

svg("korean_Tcell_contamination_1.svg", width=15, height=7)
cowplot::plot_grid(ncol = 7,
FeaturePlot(kor.t,features = "JCHAIN",reduction="tsne", cols=c("lightgrey", "red"), label=FALSE,pt.size=5),
FeaturePlot(kor.t,features = "IGHA1",reduction="tsne", cols=c("lightgrey", "red"), label=FALSE,pt.size=5),
FeaturePlot(kor.t,features = "IGHG3",reduction="tsne", cols=c("lightgrey", "red"), label=FALSE,pt.size=5),
FeaturePlot(kor.t,features = "KRT8",reduction="tsne", cols=c("lightgrey", "red"), label=FALSE,pt.size=5),
FeaturePlot(kor.t,features = "COL1A1",reduction="tsne", cols=c("lightgrey", "red"), label=FALSE,pt.size=5),
FeaturePlot(kor.t,features = "COL1A1",reduction="tsne", cols=c("lightgrey", "red"), label=FALSE,pt.size=5),
FeaturePlot(kor.t,features = "CD68",reduction="tsne", cols=c("lightgrey", "red"), label=FALSE,pt.size=5))
dev.off()


svg("korean_Tcell_contamination_2.svg", width=15, height=7)
cowplot::plot_grid(ncol = 7,
FeaturePlot(kor.t,features = "JCHAIN",reduction="tsne", cols=c("lightgrey", "red"), label=FALSE,pt.size=5),
FeaturePlot(kor.t,features = "IGHA1",reduction="tsne", cols=c("lightgrey", "red"), label=FALSE,pt.size=5),
FeaturePlot(kor.t,features = "IGHG3",reduction="tsne", cols=c("lightgrey", "red"), label=FALSE,pt.size=5),
FeaturePlot(kor.t,features = "KRT8",reduction="tsne", cols=c("lightgrey", "red"), label=FALSE,pt.size=5),
FeaturePlot(kor.t,features = "COL1A1",reduction="tsne", cols=c("lightgrey", "red"), label=FALSE,pt.size=5),
FeaturePlot(kor.t,features = "COL1A1",reduction="tsne", cols=c("lightgrey", "red"), label=FALSE,pt.size=5),
FeaturePlot(kor.t,features = "CD68",reduction="tsne", cols=c("lightgrey", "red"), label=FALSE,pt.size=5))
dev.off()


svg("korean_Tcell_cleaned_3.svg", width=15, height=7)
cowplot::plot_grid(ncol = 6,
FeaturePlot(kor.tcells,features = "JCHAIN",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5), 
FeaturePlot(kor.tcells,features = "IGHA1",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5),
FeaturePlot(kor.tcells,features = "IGHG3",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5),
FeaturePlot(kor.tcells,features = "KRT8",reduction="tsne",min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5),
FeaturePlot(kor.tcells,features = "COL1A1",reduction="tsne",min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5),
FeaturePlot(kor.tcells,features = "CD68",reduction="tsne",min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5))
dev.off()


svg("our_Tcell_cleaned_4.svg", width=15, height=7)
cowplot::plot_grid(ncol = 6,
FeaturePlot(tcells,features = "JCHAIN",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5), 
FeaturePlot(tcells,features = "IGHA1",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5),
FeaturePlot(tcells,features = "IGHG3",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5),
FeaturePlot(tcells,features = "KRT8",reduction="tsne",min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5),
FeaturePlot(tcells,features = "COL1A1",reduction="tsne",min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5),
FeaturePlot(tcells,features = "CD68",reduction="tsne",min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5))
dev.off()


#B cells

svg("korean_Bcell_contamination.svg", width=15, height=7)
cowplot::plot_grid(ncol = 2,
FeaturePlot(kor.b,features = "LYZ",reduction="tsne"  , cols=c("lightgrey", "red"), label=FALSE,pt.size=5),
FeaturePlot(kor.b,features = "CXCL8",reduction="tsne", cols=c("lightgrey", "red"), label=FALSE,pt.size=5),
FeaturePlot(kor.b,features = "CD3D",reduction="tsne", cols=c("lightgrey", "red"), label=FALSE,pt.size=5),
FeaturePlot(kor.b,features = "COL1A1",reduction="tsne", cols=c("lightgrey", "red"), label=FALSE,pt.size=5))
dev.off()


svg("korean_Bcell_cleaned.svg", width=15, height=7)
cowplot::plot_grid(ncol = 2,
FeaturePlot(kor.bcells,features = "LYZ",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5), 
FeaturePlot(kor.bcells,features = "CXCL8",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5),
FeaturePlot(kor.bcells,features = "CD3D",reduction="tsne",min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5),
FeaturePlot(kor.bcells,features = "COL1A1",reduction="tsne",min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5))
dev.off()



svg("korean_EPI_contamination.svg", width=15, height=7)
cowplot::plot_grid(ncol = 2,
FeaturePlot(kor.epi,features = "LYZ",reduction="tsne"  , cols=c("lightgrey", "red"), label=FALSE),
FeaturePlot(kor.epi,features = "CXCL8",reduction="tsne", cols=c("lightgrey", "red"), label=FALSE),
FeaturePlot(kor.epi,features = "CD3D",reduction="tsne", cols=c("lightgrey", "red"), label=FALSE),
FeaturePlot(kor.epi,features = "GZMB",reduction="tsne", cols=c("lightgrey", "red"), label=FALSE),
FeaturePlot(kor.epi,features = "GZMA",reduction="tsne", cols=c("lightgrey", "red"), label=FALSE),
FeaturePlot(kor.epi,features = "MALAT1",reduction="tsne",min.cutoff="q9",cols=c("lightgrey", "red"), label=FALSE),
FeaturePlot(kor.epi,features = "PLA2G2A",reduction="tsne", cols=c("lightgrey", "red"), label=FALSE),
FeaturePlot(kor.epi,features = "NKG7",reduction="tsne", cols=c("lightgrey", "red"), label=FALSE))
dev.off()

GZMB
GZMA
MALAT1 --> Tcells
PLA2G2A --> fibro
NKG7 --> Tcells
GNLY --> Tcells


/Users/akhaliq/Desktop/merged_epi/msih_removed/Graphs/GSVA_emt.png
Acknowledgements
This study was supported by the startup fund provided to A.M. by the Rush University Medical Center; the OCM grant to A.M by Rush University Cancer Center. Part of A.H.’s time was supported by a Merit Review Award (I01 BX000545) 

svg("korean_EPI_cleaned.svg", width=15, height=7)
cowplot::plot_grid(ncol = 2,
FeaturePlot(kor.epi,features = "LYZ",reduction="tsne"  , cols=c("lightgrey", "red"), label=FALSE),
FeaturePlot(kor.epi,features = "CXCL8",reduction="tsne", cols=c("lightgrey", "red"), label=FALSE),
FeaturePlot(kor.epi,features = "CD3D",reduction="tsne", cols=c("lightgrey", "red"), label=FALSE),
FeaturePlot(kor.epi,features = "GZMB",reduction="tsne", cols=c("lightgrey", "red"), label=FALSE),
FeaturePlot(kor.epi,features = "GZMA",reduction="tsne", cols=c("lightgrey", "red"), label=FALSE),
FeaturePlot(kor.epi,features = "MALAT1",reduction="tsne",min.cutoff="q9",cols=c("lightgrey", "red"), label=FALSE),
FeaturePlot(kor.epi,features = "PLA2G2A",reduction="tsne", cols=c("lightgrey", "red"), label=FALSE),
FeaturePlot(kor.epi,features = "NKG7",reduction="tsne", cols=c("lightgrey", "red"), label=FALSE))
dev.off()




svg("korean_Bcell_contamination.svg", width=15, height=7)
cowplot::plot_grid(ncol = 2,
FeaturePlot(kor.b,features = "LYZ",reduction="tsne"  , cols=c("lightgrey", "red"), label=FALSE,pt.size=5),
FeaturePlot(kor.b,features = "CXCL8",reduction="tsne", cols=c("lightgrey", "red"), label=FALSE,pt.size=5),
FeaturePlot(kor.b,features = "CD3D",reduction="tsne", cols=c("lightgrey", "red"), label=FALSE,pt.size=5),
FeaturePlot(kor.b,features = "COL1A1",reduction="tsne", cols=c("lightgrey", "red"), label=FALSE,pt.size=5))
dev.off()


svg("korean_Bcell_cleaned.svg", width=15, height=7)
cowplot::plot_grid(ncol = 2,
FeaturePlot(kor.bcells,features = "LYZ",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5), 
FeaturePlot(kor.bcells,features = "CXCL8",reduction="tsne"  ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5),
FeaturePlot(kor.bcells,features = "CD3D",reduction="tsne",min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5),
FeaturePlot(kor.bcells,features = "COL1A1",reduction="tsne",min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5))
dev.off()



#epithelial cell markers Smiliy

png("korean_fibro_contamination_smilie.png", units="in", width=40, height=20, res=300)
FeaturePlot(kor.stroma, reduction = "tsne",features = c("ASCL2","LGR5","LEFTY1","CDC25C","SLC26A2","FABP1","CA1","EPCAM","RBP2","SLC26A3","CCL20","TNFRSF11A","BEST4","TFF1","MUC2","HCK","TRPM5","SCGN","CHGA","CCL7","CCL13"),min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5)
dev.off()

png("korean_myloid_contamination_smilie.png", units="in", width=40, height=20, res=300)
FeaturePlot(kor.myloid, reduction = "tsne",features = c("ASCL2","LGR5","LEFTY1","CDC25C","SLC26A2","FABP1","CA1","EPCAM","RBP2","SLC26A3","CCL20","TNFRSF11A","BEST4","TFF1","MUC2","HCK","TRPM5","SCGN","CHGA","CCL7","CCL13"),min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5)
dev.off()

png("korean_T_contamination_smilie.png", units="in", width=40, height=20, res=300)
FeaturePlot(kor.t, reduction = "tsne",features = c("ASCL2","LGR5","LEFTY1","CDC25C","SLC26A2","FABP1","CA1","EPCAM","RBP2","SLC26A3","CCL20","TNFRSF11A","BEST4","TFF1","MUC2","HCK","TRPM5","SCGN","CHGA","CCL7","CCL13"),min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5)
dev.off()


#Myeloid markers smiliey 

png("korean_fibro_contamination_smilie_with_myeloid.png", units="in", width=40, height=40, res=300)
FeaturePlot(kor.stroma, reduction = "tsne",features = c("CLEC10A", "CD1C", "S100A8", "S100A9", "TPSAB1", "OSM", "NKG7", "KLRC1"),min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5)
dev.off()

png("korean_T_contamination_smilie_myeloid.png", units="in", width=40, height=20, res=300)
FeaturePlot(kor.t, reduction = "tsne",features = c("CLEC10A", "CD1C", "S100A8", "S100A9", "TPSAB1", "OSM", "NKG7", "KLRC1"),min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5)
dev.off()

png("korean_epi_contamination_smilie_with_myeloid.png", units="in", width=40, height=20, res=300)
FeaturePlot(kor.epi, reduction = "tsne",features = c("CLEC10A", "CD1C", "S100A8", "S100A9", "TPSAB1", "OSM", "NKG7", "KLRC1"),min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5)
dev.off()



#Fibro markers smiliey 

png("korean_myeloid_contamination_smilie_with_fibro.png", units="in", width=40, height=40, res=300)
FeaturePlot(kor.myloid, reduction = "tsne",features = c("RSPO3",  "CCL11", "WNT5B", "BMP4", "CHI3L1", "TAGLN", "ACTA2", "WNT2B", "SOX17", "RBP7", "VWA1", "PLVAP", "CDH5", "CD36", "MADCAM1", "DARC", "RGS5", "SOX10", "S100B", "CD68", "XCR1", "CLEC9A"),min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5)
dev.off()

png("korean_T_contamination_smilie_with_fibro.png", units="in", width=40, height=20, res=300)
FeaturePlot(kor.t, reduction = "tsne",features = c("RSPO3", "CCL11", "WNT5B", "BMP4", "CHI3L1", "TAGLN", "ACTA2", "WNT2B", "SOX17", "RBP7", "VWA1", "PLVAP", "CDH5", "CD36", "MADCAM1", "DARC", "RGS5", "SOX10", "S100B", "CD68", "XCR1", "CLEC9A"),min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5)
dev.off()

png("korean_epi_contamination_smilie_with_fibro.png", units="in", width=40, height=20, res=300)
FeaturePlot(kor.epi, reduction = "tsne",features = c("RSPO3", "CCL11", "WNT5B", "BMP4", "CHI3L1", "TAGLN", "ACTA2", "WNT2B", "SOX17", "RBP7", "VWA1", "PLVAP", "CDH5", "CD36", "MADCAM1", "DARC", "RGS5", "SOX10", "S100B", "CD68", "XCR1", "CLEC9A"),min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5)
dev.off()



#T markers smilliey 


png("korean_myeloid_contamination_smilie_with_T.png", units="in", width=40, height=40, res=300)
FeaturePlot(kor.myloid, reduction = "tsne",features = c("FOXP3", "CTLA4", "CD8B", "CXCR6", "CD3D"),min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5)
dev.off()

png("korean_Fibro_contamination_smilie_with_T.png", units="in", width=40, height=20, res=300)
FeaturePlot(kor.stroma, reduction = "tsne",features = c("FOXP3", "CTLA4", "CD8B", "CXCR6", "CD3D"),min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5)
dev.off()

png("korean_epi_contamination_smilie_with_T.png", units="in", width=40, height=20, res=300)
FeaturePlot(kor.epi, reduction = "tsne",features = c("FOXP3", "CTLA4", "CD8B", "CXCR6", "CD3D"),min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5)
dev.off()


#B markers smiliey 

png("korean_fibro_contamination_smilie_with_Bcells.png", units="in", width=40, height=20, res=300)
FeaturePlot(kor.stroma, reduction = "tsne",features = c("MZB1", "IGHA1", "SELL", "CD19", "AICDA"),min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5)
dev.off()

png("korean_myloid_contamination_smilie_with_bcells.png", units="in", width=40, height=20, res=300)
FeaturePlot(kor.myloid, reduction = "tsne",features = c("MZB1", "IGHA1", "SELL", "CD19", "AICDA"),min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5)
dev.off()

png("korean_T_contamination_smilie_with_Bcells.png", units="in", width=40, height=20, res=300)
FeaturePlot(kor.t, reduction = "tsne",features = c("MZB1", "IGHA1", "SELL", "CD19", "AICDA"),min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5)
dev.off()

png("korean_epi_contamination_smilie_with_Bcells.png", units="in", width=40, height=20, res=300)
FeaturePlot(kor.epi, reduction = "tsne",features = c("MZB1", "IGHA1", "SELL", "CD19", "AICDA"),min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,pt.size=5)
dev.off()



png("epi_contamination_compartment.png", units="in", width=20, height=6, res=300)
cowplot::plot_grid(ncol = 6,
FeaturePlot(kor.epi,features = "KRT18" ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,reduction="tsne",pt.size=5)+ labs(title = "Epithelial Cells:KRT18"),
FeaturePlot(kor.epi,features = "CD3D" ,min.cutoff="q9", cols=c("lightgrey","red"), label=FALSE,reduction="tsne",pt.size=5)+ labs(title = "T Cells:CD3D"),
FeaturePlot(kor.epi,features = "CD79A" ,min.cutoff="q9", cols=c("lightgrey","red"), label=FALSE,reduction="tsne",pt.size=5)+ labs(title = "B Cells:CD79A"),
FeaturePlot(kor.epi,features = "LYZ" ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,reduction="tsne",pt.size=5)+ labs(title = "Myeloid:LYZ"),
FeaturePlot(kor.epi,features = "COL1A1" ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,reduction="tsne",pt.size=5)+ labs(title = "Fibroblast:COL1A1"),
FeaturePlot(kor.epi,features = "CLDN5" ,min.cutoff="q9", cols=c("lightgrey","red"), label=FALSE,reduction="tsne",pt.size=5)+ labs(title = "Endothelial:CLDN5"))
dev.off()



png("fibro_contamination_compartment.png", units="in", width=20, height=6, res=300)
cowplot::plot_grid(ncol = 6,
FeaturePlot(kor.stroma,features = "KRT18" ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,reduction="tsne",pt.size=5)+ labs(title = "Epithelial Cells:KRT18"),
FeaturePlot(kor.stroma,features = "CD3D" ,min.cutoff="q9", cols=c("lightgrey","red"), label=FALSE,reduction="tsne",pt.size=5)+ labs(title = "T Cells:CD3D"),
FeaturePlot(kor.stroma,features = "CD79A" ,min.cutoff="q9", cols=c("lightgrey","red"), label=FALSE,reduction="tsne",pt.size=5)+ labs(title = "B Cells:CD79A"),
FeaturePlot(kor.stroma,features = "LYZ" ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,reduction="tsne",pt.size=5)+ labs(title = "Myeloid:LYZ"),
FeaturePlot(kor.stroma,features = "COL1A1" ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,reduction="tsne",pt.size=5)+ labs(title = "Fibroblast:COL1A1"),
FeaturePlot(kor.stroma,features = "CLDN5" ,min.cutoff="q9", cols=c("lightgrey","red"), label=FALSE,reduction="tsne",pt.size=5)+ labs(title = "Endothelial:CLDN5"))
dev.off()



png("myeloid_contamination_compartment.png", units="in", width=20, height=6, res=300)
cowplot::plot_grid(ncol = 6,
FeaturePlot(kor.myloid,features = "KRT18" ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,reduction="tsne",pt.size=5)+ labs(title = "Epithelial Cells:KRT18"),
FeaturePlot(kor.myloid,features = "CD3D" ,min.cutoff="q9", cols=c("lightgrey","red"), label=FALSE,reduction="tsne",pt.size=5)+ labs(title = "T Cells:CD3D"),
FeaturePlot(kor.myloid,features = "CD79A" ,min.cutoff="q9", cols=c("lightgrey","red"), label=FALSE,reduction="tsne",pt.size=5)+ labs(title = "B Cells:CD79A"),
FeaturePlot(kor.myloid,features = "LYZ" ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,reduction="tsne",pt.size=5)+ labs(title = "Myeloid:LYZ"),
FeaturePlot(kor.myloid,features = "COL1A1" ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,reduction="tsne",pt.size=5)+ labs(title = "Fibroblast:COL1A1"),
FeaturePlot(kor.myloid,features = "CLDN5" ,min.cutoff="q9", cols=c("lightgrey","red"), label=FALSE,reduction="tsne",pt.size=5)+ labs(title = "Endothelial:CLDN5"))
dev.off()



png("tcell_contamination_compartment.png", units="in", width=20, height=6, res=300)
cowplot::plot_grid(ncol = 6,
FeaturePlot(kor.t,features = "KRT18" ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,reduction="tsne",pt.size=5)+ labs(title = "Epithelial Cells:KRT18"),
FeaturePlot(kor.t,features = "CD3D" ,min.cutoff="q9", cols=c("lightgrey","red"), label=FALSE,reduction="tsne",pt.size=5)+ labs(title = "T Cells:CD3D"),
FeaturePlot(kor.t,features = "CD79A" ,min.cutoff="q9", cols=c("lightgrey","red"), label=FALSE,reduction="tsne",pt.size=5)+ labs(title = "B Cells:CD79A"),
FeaturePlot(kor.t,features = "LYZ" ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,reduction="tsne",pt.size=5)+ labs(title = "Myeloid:LYZ"),
FeaturePlot(kor.t,features = "COL1A1" ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,reduction="tsne",pt.size=5)+ labs(title = "Fibroblast:COL1A1"),
FeaturePlot(kor.t,features = "CLDN5" ,min.cutoff="q9", cols=c("lightgrey","red"), label=FALSE,reduction="tsne",pt.size=5)+ labs(title = "Endothelial:CLDN5"))
dev.off()


png("bcell_contamination_compartment.png", units="in", width=20, height=6, res=300)
cowplot::plot_grid(ncol = 6,
FeaturePlot(kor.b,features = "KRT18" ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,reduction="tsne",pt.size=5)+ labs(title = "Epithelial Cells:KRT18"),
FeaturePlot(kor.b,features = "CD3D" ,min.cutoff="q9", cols=c("lightgrey","red"), label=FALSE,reduction="tsne",pt.size=5)+ labs(title = "T Cells:CD3D"),
FeaturePlot(kor.b,features = "CD79A" ,min.cutoff="q9", cols=c("lightgrey","red"), label=FALSE,reduction="tsne",pt.size=5)+ labs(title = "B Cells:CD79A"),
FeaturePlot(kor.b,features = "LYZ" ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,reduction="tsne",pt.size=5)+ labs(title = "Myeloid:LYZ"),
FeaturePlot(kor.b,features = "COL1A1" ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE,reduction="tsne",pt.size=5)+ labs(title = "Fibroblast:COL1A1"),
FeaturePlot(kor.b,features = "CLDN5" ,min.cutoff="q9", cols=c("lightgrey","red"), label=FALSE,reduction="tsne",pt.size=5)+ labs(title = "Endothelial:CLDN5"))
dev.off()



#######


rsvg-convert /Users/akhaliq/Desktop/new_plots/Manuscript/Figures/Main_figures/mf3/mf3.svg -w 210 -h 297 -f png -o /Users/akhaliq/Desktop/new_plots/Manuscript/Figures/Main_figures/mf3/test-skewed.png


#######

#merging 2 seurat obj


#####



epi.combined <- merge(epi.tumor, y = kor.epi.tumor, add.cell.ids = c("RUSH-CRC", "LEE-CRC"), project = "CRC")

Removed samples which are less than 10 cells

alldata_subset<- subset(x = alldata, subset = orig.ident==c("SMC01-T","SMC02-T","SMC03-T","SMC04-T","SMC07-T","SMC08-T","SMC09-T","SMC10-T","SMC11-T","SMC14-T","SMC15-T","SMC16-T","SMC17-T","SMC18-T","SMC19-T","SMC20-T","SMC21-T","SMC22-T","SMC23-T","T_cac1","T_cac2","T_cac3","T_cac4","T_cac6","T_cac7","T_cac8","T_cac9","T_cac10","T_cac11","T_cac12","T_cac13","T_cac14","T_cac15","T_cac16"))



#### 
Merge PDFs 
####

in Linux


pdfunite 1_SUPPLEMENTARY_FIGURE_DISCRIPTIONS.pdf 1.1_sf1_new.pdf 2.1_sf2.pdf 2.2_sf2_desc.pdf 3_merged_suppli_figures.pdf 3.1_sf3.pdf 3.2_sf3_desc.pdf 4.1_sf4.pdf 4.2_sf4_desc.pdf 5.1_sf5.pdf 5.2_sf5_desc.pdf 6.1_sf6.pdf 6.2_sf6_desc.pdf 7.1_sf7.pdf 7.2_sf7_desc.pdf 8.1_sf8.pdf 8.2_sf8_desc.pdf 9.1_sf9.pdf 9.2_sf9_desc.pdf 1.2_sf1_desc.pdf 10.1_sf10.pdf 10.2_sf10_desc.pdf 11.1_sf11.pdf 11.2_sf11_desc.pdf 12.1_sf12.pdf 12.2_sf12_desc.pdf 13.1_sf13_new.pdf 13.2_sf13_desc.pdf 14.1_sf14_new.pdf 14.2_sf14_desc.pdf 15.1_sf15_new.pdf 15.2_sf15_desc.pdf 16.1_sf16_new.pdf 16.2_sf16_desc.pdf 17.1_sf17.pdf 17.2_sf17_desc.pdf 18.1_sf18.pdf 18.2_sf18_desc.pdf 19.1_sf19.pdf 19.2_sf19_desc.pdf merged_suppli_figures_01312022.pdf


pdfunite mf1_new.sevgi.pdf Fig1_desc.pdf mf2.pdf Fig2_desc.pdf mf3.pdf Fig3_desc.pdf mf4.pdf Fig4_desc.pdf mf5.pdf Fig5_desc.pdf mf6.pdf Fig6_desc.pdf mf7.pdf Fig7_desc.pdf /Users/akhaliq/Desktop/genome_biology/suppli_figures/merged_suppli_figures_01312022.pdf Figures_merged.all.pdf

######
Conda create new enviroment and add new version of python in it
######
conda create -n infercnv python=3.7.6





####

Reading h5ad (anndata) file in python


adata = sc.read("/Users/akhaliq/Desktop/drugs_org/epi.tumor.h5ad")


# Convert Seurat clusters to anndata

# lets first subset the data 
allexcp.epi <- SubsetData(alldata_new, subset.name = "seurat_clusters_new", accept.value = c("B-cells","CAF","Endothelial","Myeloid","T-cells"))
epi.tumor <- SubsetData(epi, subset.name = "Condition", accept.value = "Tumor")

# merge allexcp.epi and epi.tumor

allcomb.excep.epinor<- merge(allexcp.epi1, y = epi.tumor1, add.cell.ids = c("all", "epi"),project = "CRC")

# now run all the tsnes and umaps
alldata_new <- NormalizeData(alldata_new, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE,assay="RNA")
alldata_new <- FindVariableFeatures(alldata_new, selection.method = "vst", nfeatures = 2000, verbose = FALSE,assay="RNA")
alldata_new <- ScaleData(alldata_new, verbose = FALSE)

alldata_new <- RunPCA(alldata_new, verbose = FALSE, features = VariableFeatures(object = alldata_new))
alldata_new <- RunTSNE(object = alldata_new, dims.use = 1:30)
alldata_new <- RunUMAP(alldata_new, dims = 1:30)  
alldata_new <- FindNeighbors(alldata_new, dims = 1:30, force.recalc = T)

alldata_new <- FindClusters(alldata_new, resolution = 0.8, algorithm = 1, force.recalc=T)

alldata_new$seurat_clusters <- as.factor(as.numeric(as.character(alldata_new$seurat_clusters)) + 1)
Idents(alldata_new) <- "seurat_clusters"


# now convert seurat to anndata format(h5ad)
#now switch on the conda activate infercnv
# R terminal

library(Seurat)

alldata_new <- readRDS("alldata_new.rds")

library(sceasy)

sceasy::convertFormat(alldata_new, from="seurat", to="anndata",outFile='alldata.h5ad')

# now read the anndata into python 
# open terminal and conda activate infercnv

python
import scanpy as sc
import infercnvpy as cnv
import matplotlib.pyplot as plt

sc.settings.set_figure_params(figsize=(20, 30))


sc.logging.print_header()

adata = sc.read("/Users/akhaliq/Desktop/drugs_org/alldata.h5ad")

sc.pl.umap(adata, color="seurat_clusters_new")

#Genomic positions need to be added. There need to be `chromosome`, `start`, and `end` columns in `adata.var` in adata.

infercnvpy.io.genomic_position_from_gtf(/Users/akhaliq/Desktop/drugs_org/gencode.v39.annotation.gtf, adata=None, *, gtf_gene_id='gene_name', adata_gene_id=None, inplace=True)



cnv.tl.infercnv(
    adata,
    reference_key="seurat_clusters_new",
    reference_cat=[
        "B-cells",
        "CAF",
        "Endothelial",
        "Myeloid",
        "T-cells"
    ],
    window_size=250,
)


#######
Conda activate
#######


source /home/akhaliq/miniconda3/bin/activate
conda create --name seurat
conda activate seurat

######
CellphoneDB analysis
######


cpdb.combined <- merge(tcells, y = c(fibro_5subtypes, myeloid))
write.csv(cpdb.combined@meta.data,"combined_meta.csv")
# modify the csv to include a new metadata coloumn comb_clusters
# add that colomn in the seurat obj

meta <- read.csv("combined_meta.csv",header=TRUE, sep=',', row.names=1)
cpdb.combined<- AddMetaData(cpdb.combined, meta)
kable(table(cpdb.combined@meta.data$comb_clusters))


####
 "/Users/akhaliq/Desktop/genome_biology/cellphonedb_analysis"
conda activate cpdb
cellphonedb method statistical_analysis alldata_meta.txt log_counts.txt --counts-data gene_name --threads 8 
#[ ][CORE][06/06/21-03:13:45][INFO] [Cluster Statistical Analysis] Threshold:0.1 Iterations:1000 Debug-seed:-1 Threads:8 Precision:3
cellphonedb plot dot_plot --rows rows.txt --columns columns.txt

#Converting Seurat object to LOG 2 - Log normalised count data
log_counts= log2(exp(as.matrix(GetAssayData(object = cpdb.combined, slot = "data"))))
### just another way to get log counts
#log_counts= as.matrix(epi.subset@assays$RNA@data)

library(tibble)
log_counts1 <- tibble::rownames_to_column(as.data.frame(log_counts), "Gene")
write.table(log_counts1, file="log_counts.txt", sep="\t",append=F, row.names = F)

write.table(cpdb.combined@meta.data, file="metadata.txt", sep="\t",append=F)

conda activate cpdb
cellphonedb method statistical_analysis metadata.txt log_counts.txt --counts-data gene_name --threads 8 
#[ ][CORE][06/06/21-03:13:45][INFO] [Cluster Statistical Analysis] Threshold:0.1 Iterations:1000 Debug-seed:-1 Threads:8 Precision:3
cellphonedb plot dot_plot --rows rows.txt --columns columns.txt

cellphonedb plot dot_plot

cellphonedb plot heatmap_plot metadata.txt

cellphonedb plot dot_plot --rows rows_col.txt --columns columns_caf.txt --output-name cafs_1.pdf

cellphonedb plot dot_plot --rows rows_col.txt --columns columns_tam.txt --output-name tam_1.pdf

cellphonedb plot dot_plot --rows rows_col.txt --columns columns_cd8.txt --output-name cd8_1.pdf

--output-name

significant_means.txt

/Users/akhaliq/Desktop/genome_biology/cellphonedb_analysis/out/significant_means.txt


sig <- read.csv("/Users/akhaliq/Desktop/genome_biology/cellphonedb_analysis/out/significant_means_test.csv",header=TRUE, sep=',', row.names=1)



#######

library(Seurat)
library(ktplots)
cpdb = readRDS("/Users/akhaliq/Desktop/genome_biology/cellphonedb_analysis/input_cpdb.rds")
sce <- Seurat::as.SingleCellExperiment(cpdb)

pvals <- read.delim("out/pvalues.txt", check.names = FALSE)
means <- read.delim("out/means.txt", check.names = FALSE)
#chemokines


#Keeps Significant Interactions

# TAM|CAF-S4

TGFB-myCAF|CAF-S4
detox-IL-iCAF|CAF-S4
ecm-myCAF|CAF-S4
wound-myCAF|CAF-S4


p1<-plot_cpdb(cell_type1 = 'CAF-S4', cell_type2 = 'TAM|TGFB-myCAF|detox-IL-iCAF|ecm-myCAF|wound-myCAF|CD1C-DENDRITIC CELLS|CD4-CXCL13', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE) + 
small_axis(fontsize = 15) + small_grid() + small_guide() + small_legend(fontsize = 5) # some helper functions included in ktplots to help with the plotting


png("CAFS4.png", units="in", width=20, height=45, res=300)
p1
dev.off()

pdf("CAFS4_TAM_5subtypes.pdf")
p1
dev.off()

# CD8

p2<-plot_cpdb(cell_type1 = 'CD8', cell_type2 = 'TAM|TGFB-myCAF|detox-IL-iCAF|ecm-myCAF|wound-myCAF|CAF-S4|CD1C-DENDRITIC CELLS|CXCL13|TIGIT', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE) + 
small_axis(fontsize = 25) + small_grid() + small_guide() + small_legend(fontsize = 10) # some helper functions included in ktplots to help with the plotting


png("CD8_Cellphonedb.png", units="in", width=50, height=60, res=300)
p2
dev.off()


# TAM

p3<-plot_cpdb(cell_type1 = 'TAM', cell_type2 = 'TGFB-myCAF|detox-IL-iCAF|ecm-myCAF|wound-myCAF|CAF-S4|CD1C-DENDRITIC CELLS|CXCL13|TIGIT', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE) + 
small_axis(fontsize = 25) + small_grid() + small_guide() + small_legend(fontsize = 10) # some helper functions included in ktplots to help with the plotting


png("TAM_Cellphonedb.png", units="in", width=20, height=50, res=300)
p3
dev.off()


####for future reference, genes or gene.family can be specified, not both.
gene.family can be one of the following:
[1] "chemokines"    "Th1"           "Th2"           "Th17"         
[5] "Treg"          "costimulatory" "coinhibitory"  "niche" 

########

p4<-plot_cpdb(cell_type1 = 'CAF-S4', cell_type2 = 'TAM|TGFB-myCAF|detox-IL-iCAF|ecm-myCAF|wound-myCAF|CD1C-DENDRITIC CELLS|CD4-CXCL13', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE,gene.family='chemokines') + 
small_axis(fontsize = 15) + small_grid() + small_guide() + small_legend(fontsize = 5) # some helper functions included in ktplots to help with the plotting


png("CAFS4.chemokines.png", units="in", width=10, height=15, res=300)
p4
dev.off()


p5<-plot_cpdb(cell_type1 = 'CAF-S4', cell_type2 = 'TAM|TGFB-myCAF|detox-IL-iCAF|ecm-myCAF|wound-myCAF|CD1C-DENDRITIC CELLS|CD4-CXCL13', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE,gene.family='costimulatory') + 
small_axis(fontsize = 15) + small_grid() + small_guide() + small_legend(fontsize = 5) # some helper functions included in ktplots to help with the plotting


png("CAFS4.costimulatory.png", units="in", width=10, height=15, res=300)
p5
dev.off()



p6<-plot_cpdb(cell_type1 = 'CAF-S4', cell_type2 = 'TAM|TGFB-myCAF|detox-IL-iCAF|ecm-myCAF|wound-myCAF|CD1C-DENDRITIC CELLS|CD4-CXCL13', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE,gene.family='coinhibitory') + 
small_axis(fontsize = 15) + small_grid() + small_guide() + small_legend(fontsize = 5) # some helper functions included in ktplots to help with the plotting


png("CAFS4.coinhibitory.png", units="in", width=10, height=15, res=300)
p6
dev.off()



# CD8

p7<-plot_cpdb(cell_type1 = 'CD8', cell_type2 = 'TAM|TGFB-myCAF|detox-IL-iCAF|ecm-myCAF|wound-myCAF|CAF-S4|CD1C-DENDRITIC CELLS|CXCL13|TIGIT', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE,gene.family='chemokines') + 
small_axis(fontsize = 25) + small_grid() + small_guide() + small_legend(fontsize = 10) # some helper functions included in ktplots to help with the plotting


png("CD8.chemokines.png", units="in", width=55, height=20, res=300)
p7
dev.off()


p8<-plot_cpdb(cell_type1 = 'CD8', cell_type2 = 'TAM|TGFB-myCAF|detox-IL-iCAF|ecm-myCAF|wound-myCAF|CAF-S4|CD1C-DENDRITIC CELLS|CXCL13|TIGIT', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE,gene.family='costimulatory') + 
small_axis(fontsize = 25) + small_grid() + small_guide() + small_legend(fontsize = 10) # some helper functions included in ktplots to help with the plotting


png("CD8.costimulatory.png", units="in", width=55, height=30, res=300)
p8
dev.off()


p9<-plot_cpdb(cell_type1 = 'CD8', cell_type2 = 'TAM|TGFB-myCAF|detox-IL-iCAF|ecm-myCAF|wound-myCAF|CAF-S4|CD1C-DENDRITIC CELLS|CXCL13|TIGIT', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE,gene.family='coinhibitory') + 
small_axis(fontsize = 25) + small_grid() + small_guide() + small_legend(fontsize = 10) # some helper functions included in ktplots to help with the plotting


png("CD8.coinhibitory.png", units="in", width=55, height=30, res=300)
p9
dev.off()



# TAM

p10<-plot_cpdb(cell_type1 = 'TAM', cell_type2 = 'TGFB-myCAF|detox-IL-iCAF|ecm-myCAF|wound-myCAF|CAF-S4|CD1C-DENDRITIC CELLS|CXCL13|TIGIT', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE,gene.family='chemokines') + 
small_axis(fontsize = 25) + small_grid() + small_guide() + small_legend(fontsize = 10) # some helper functions included in ktplots to help with the plotting


png("TAM.chemokines.png", units="in", width=20, height=10, res=300)
p10
dev.off()



p10<-plot_cpdb(cell_type1 = 'TAM', cell_type2 = 'TGFB-myCAF|detox-IL-iCAF|ecm-myCAF|wound-myCAF|CAF-S4|CD1C-DENDRITIC CELLS|CXCL13|TIGIT', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE,gene.family='costimulatory') + 
small_axis(fontsize = 25) + small_grid() + small_guide() + small_legend(fontsize = 10) # some helper functions included in ktplots to help with the plotting


png("TAM.costimulatory.png", units="in", width=20, height=25, res=300)
p10
dev.off()



p11<-plot_cpdb(cell_type1 = 'TAM', cell_type2 = 'TGFB-myCAF|detox-IL-iCAF|ecm-myCAF|wound-myCAF|CAF-S4|CD1C-DENDRITIC CELLS|CXCL13|TIGIT', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE,gene.family='coinhibitory') + 
small_axis(fontsize = 25) + small_grid() + small_guide() + small_legend(fontsize = 10) # some helper functions included in ktplots to help with the plotting


png("TAM.coinhibitory.png", units="in", width=20, height=25, res=300)
p11
dev.off()


CD8_Cellphonedb.pdf
CAFS4.pdf
TAM_Cellphonedb.pdf
CAFS4.chemokines.pdf
CAFS4.costimulatory.pdf
CAFS4.coinhibitory.pdf
CD8.chemokines.pdf
CD8.costimulatory.pdf
CD8.coinhibitory.pdf
TAM.chemokines.pdf
TAM.costimulatory.pdf
TAM.coinhibitory.pdf

pdfunite CAFS4.pdf CAFS4.coinhibitory.pdf CAFS4.costimulatory.pdf CAFS4.chemokines.pdf TAM_Cellphonedb.pdf TAM.coinhibitory.pdf TAM.costimulatory.pdf TAM.chemokines.pdf CD8_Cellphonedb.pdf CD8.coinhibitory.pdf CD8.costimulatory.pdf CD8.chemokines.pdf CellphoneDB_results.pdf

# convert .png to pdf 
library(magick)

all_images_1 <- purrr::reduce(
    purrr::map(all_images,image_read),
    c
)

image_write( "CAFS4.png", "CAFS4.coinhibitory.png", "CAFS4.costimulatory.png", "CAFS4.chemokines.png", "TAM_Cellphonedb.png", "TAM.coinhibitory.png", "TAM.costimulatory.png", "TAM.chemokines.png", "CD8_Cellphonedb.png", "CD8.coinhibitory.png", "CD8.costimulatory.png", "CD8.chemokines.png", format = "pdf", "check.pdf")

#####




p12<-plot_cpdb(cell_type1 = 'TAM', cell_type2 = 'TGFB-myCAF|detox-IL-iCAF|ecm-myCAF|wound-myCAF|CD1C-DENDRITIC CELLS|CD4-CXCL13', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE) + 
small_axis(fontsize = 15) + small_grid() + small_guide() + small_legend(fontsize = 5) # some helper functions included in ktplots to help with the plotting


png("TAM-CAFS1.png", units="in", width=20, height=45, res=300)
p12
dev.off()



p13<-plot_cpdb(cell_type1 = 'TAM', cell_type2 = 'CD8', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE) + 
small_axis(fontsize = 15) + small_grid() + small_guide() + small_legend(fontsize = 5) # some helper functions included in ktplots to help with the plotting


png("TAM-CD8.png", units="in", width=20, height=45, res=300)
p13
dev.off()


p14<-plot_cpdb(cell_type1 = 'TGFB-myCAF|detox-IL-iCAF|ecm-myCAF|wound-myCAF', cell_type2 = 'CD8', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE) + 
small_axis(fontsize = 15) + small_grid() + small_guide() + small_legend(fontsize = 5) # some helper functions included in ktplots to help with the plotting


png("TAM-CD8.png", units="in", width=40, height=45, res=300)
p14
dev.off()



#TGFB-CD8
p15<-plot_cpdb(cell_type1 = 'TGFB-myCAF', cell_type2 = 'CD8', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE) + 
small_axis(fontsize = 15) + small_grid() + small_guide() + small_legend(fontsize = 5) # some helper functions included in ktplots to help with the plotting


png("TGFB-CD8.png", units="in", width=20, height=25, res=300)
p15
dev.off()



#detox-IL-iCAF-CD8
p16<-plot_cpdb(cell_type1 = 'detox-IL-iCAF', cell_type2 = 'CD8', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE) + 
small_axis(fontsize = 15) + small_grid() + small_guide() + small_legend(fontsize = 5) # some helper functions included in ktplots to help with the plotting


png("detox-IL-iCAF-CD8.png", units="in", width=20, height=25, res=300)
p16
dev.off()


#ecm-myCAF-CD8
p17<-plot_cpdb(cell_type1 = 'ecm-myCAF', cell_type2 = 'CD8', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE) + 
small_axis(fontsize = 15) + small_grid() + small_guide() + small_legend(fontsize = 5) # some helper functions included in ktplots to help with the plotting


png("ecm-myCAF-CD8.png", units="in", width=20, height=25, res=300)
p17
dev.off()


#wound-myCAF-CD8
p18<-plot_cpdb(cell_type1 = 'wound-myCAF', cell_type2 = 'CD8', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE) + 
small_axis(fontsize = 15) + small_grid() + small_guide() + small_legend(fontsize = 5) # some helper functions included in ktplots to help with the plotting


png("wound-myCAF-CD8.png", units="in", width=20, height=25, res=300)
p18
dev.off()



#TGFB-TAM
p19<-plot_cpdb(cell_type1 = 'TGFB-myCAF', cell_type2 = 'TAM', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE) + 
small_axis(fontsize = 15) + small_grid() + small_guide() + small_legend(fontsize = 5) # some helper functions included in ktplots to help with the plotting


png("TGFB-TAM.png", units="in", width=10, height=25, res=300)
p19
dev.off()



#detox-IL-iCAF-TAM

p20<-plot_cpdb(cell_type1 = 'detox-IL-iCAF', cell_type2 = 'TAM', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE) + 
small_axis(fontsize = 15) + small_grid() + small_guide() + small_legend(fontsize = 5) # some helper functions included in ktplots to help with the plotting


png("detox-IL-iCAF-TAM.png", units="in", width=10, height=25, res=300)
p20
dev.off()


#ecm-myCAF-TAM
p21<-plot_cpdb(cell_type1 = 'ecm-myCAF', cell_type2 = 'TAM', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE) + 
small_axis(fontsize = 15) + small_grid() + small_guide() + small_legend(fontsize = 5) # some helper functions included in ktplots to help with the plotting


png("ecm-myCAF-TAM.png", units="in", width=10, height=25, res=300)
p21
dev.off()


#wound-myCAF-TAM
p22<-plot_cpdb(cell_type1 = 'wound-myCAF', cell_type2 = 'TAM', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE) + 
small_axis(fontsize = 15) + small_grid() + small_guide() + small_legend(fontsize = 5) # some helper functions included in ktplots to help with the plotting


png("wound-myCAF-TAM.png", units="in", width=10, height=25, res=300)
p22
dev.off()


#TAM-CD8
p23<-plot_cpdb(cell_type1 = 'TAM', cell_type2 = 'CD8', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE) + 
small_axis(fontsize = 15) + small_grid() + small_guide() + small_legend(fontsize = 5) # some helper functions included in ktplots to help with the plotting


png("TAM-CD8.png", units="in", width=20, height=25, res=300)
p23
dev.off()


##
#TAM-CAFS1

p24<-plot_cpdb(cell_type1 = 'TAM', cell_type2 = 'TGFB-myCAF|detox-IL-iCAF|ecm-myCAF|wound-myCAF|CD1C-DENDRITIC CELLS|CD4-CXCL13', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE,gene.family='chemokines') + 
small_axis(fontsize = 15) + small_grid() + small_guide() + small_legend(fontsize = 5) # some helper functions included in ktplots to help with the plotting


png("TAM-CAFS1.chemokines.png", units="in", width=10, height=15, res=300)
p24
dev.off()



p25<-plot_cpdb(cell_type1 = 'TAM', cell_type2 = 'TGFB-myCAF|detox-IL-iCAF|ecm-myCAF|wound-myCAF|CD1C-DENDRITIC CELLS|CD4-CXCL13', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE,gene.family='costimulatory') + 
small_axis(fontsize = 15) + small_grid() + small_guide() + small_legend(fontsize = 5) # some helper functions included in ktplots to help with the plotting


png("TAM-CAFS1.costimulatory.png", units="in", width=10, height=15, res=300)
p25
dev.off()



p26<-plot_cpdb(cell_type1 = 'TAM', cell_type2 = 'TGFB-myCAF|detox-IL-iCAF|ecm-myCAF|wound-myCAF|CD1C-DENDRITIC CELLS|CD4-CXCL13', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE,gene.family='coinhibitory') + 
small_axis(fontsize = 15) + small_grid() + small_guide() + small_legend(fontsize = 5) # some helper functions included in ktplots to help with the plotting


png("TAM-CAFS1.coinhibitory.png", units="in", width=10, height=15, res=300)
p26
dev.off()

pdfunite TAM-CAFS1.pdf TAM-CAFS1.chemokines.pdf TAM-CAFS1.costimulatory.pdf TAM-CAFS1.coinhibitory.pdf TAM-CAFS1.Cellphonedb.pdf

#
#TAM-CD8


p27<-plot_cpdb(cell_type1 = 'TAM', cell_type2 = 'CD8', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE,gene.family='chemokines') + 
small_axis(fontsize = 15) + small_grid() + small_guide() + small_legend(fontsize = 5) # some helper functions included in ktplots to help with the plotting


png("TAM-CD8.chemokines.png", units="in", width=10, height=15, res=300)
p27
dev.off()


p28<-plot_cpdb(cell_type1 = 'TAM', cell_type2 = 'CD8', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE,gene.family='costimulatory') + 
small_axis(fontsize = 15) + small_grid() + small_guide() + small_legend(fontsize = 5) # some helper functions included in ktplots to help with the plotting


png("TAM-CD8.costimulatory.png", units="in", width=10, height=15, res=300)
p28
dev.off()


p29<-plot_cpdb(cell_type1 = 'TAM', cell_type2 = 'CD8', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE,gene.family='coinhibitory') + 
small_axis(fontsize = 15) + small_grid() + small_guide() + small_legend(fontsize = 5) # some helper functions included in ktplots to help with the plotting


png("TAM-CD8.coinhibitory.png", units="in", width=10, height=15, res=300)
p29
dev.off()


TAM-CD8.chemokines.pdf
TAM-CD8.costimulatory.pdf
TAM-CD8.coinhibitory.pdf
TAM-CD8.pdf

pdfunite TAM-CD8.pdf TAM-CD8.chemokines.pdf TAM-CD8.costimulatory.pdf TAM-CD8.coinhibitory.pdf TAM-CD8.Cellphonedb.pdf


TGFB-CD8.pdf
detox-IL-iCAF-CD8.pdf
ecm-myCAF-CD8.pdf
wound-myCAF-CD8.pdf

pdfunite TGFB-CD8.pdf detox-IL-iCAF-CD8.pdf ecm-myCAF-CD8.pdf wound-myCAF-CD8.pdf CAFS1-CD8.Cellphonedb.pdf


p29<-plot_cpdb(cell_type1 = 'TAM', cell_type2 = 'CD8', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE,genes=c("SELL","TIGIT","NECTIN","CD86","CD40","CTLA")) + 
small_axis(fontsize = 15) + small_grid() + small_guide() + small_legend(fontsize = 5) # some helper functions included in ktplots to help with the plotting



plot_cpdb(cell_type1 = 'TAM', cell_type2 = 'TGFB-myCAF|detox-IL-iCAF|ecm-myCAF|wound-myCAF', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE,genes=c("TIGIT","NECTIN","CD226","HLA","TNF","CHI3L1")) + 
small_axis(fontsize = 15) + small_grid() + small_guide() + small_legend(fontsize = 5) #





plot_cpdb(cell_type1 = 'TAM', cell_type2 = 'CD8-EFFECTOR 1|CD8-EFFECTOR 2|CD8-EFFECTOR 3', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE,gene.family='coinhibitory'|'costimulatory'|'chemokines') + 
small_axis(fontsize = 15) + small_grid() + small_guide() + small_legend(fontsize = 5) #



png("TAM-CD8.coinhibitory.png", units="in", width=10, height=15, res=300)

pl1 <- plot_cpdb(cell_type1 = 'TAM', cell_type2 = 'CD8-EFFECTOR 1|CD8-EFFECTOR 2|CD8-EFFECTOR 3', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE,gene.family='coinhibitory') + 
small_axis(fontsize = 15) + small_grid() + small_guide() + small_legend(fontsize = 5) #
dev.off()

png("TAM-CD8.costimulatory.png", units="in", width=10, height=15, res=300)
pl2 <- plot_cpdb(cell_type1 = 'TAM', cell_type2 = 'CD8-EFFECTOR 1|CD8-EFFECTOR 2|CD8-EFFECTOR 3', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE,gene.family='costimulatory') + 
small_axis(fontsize = 15) + small_grid() + small_guide() + small_legend(fontsize = 5) #
dev.off()

png("TAM-CD8.chemokines.png", units="in", width=10, height=15, res=300)
pl3 <- plot_cpdb(cell_type1 = 'TAM', cell_type2 = 'CD8-EFFECTOR 1|CD8-EFFECTOR 2|CD8-EFFECTOR 3', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE,gene.family='chemokines') + 
small_axis(fontsize = 15) + small_grid() + small_guide() + small_legend(fontsize = 5) #
dev.off()


png("TAM-CD8.cytokines.png", units="in", width=10, height=15, res=300)
pl4 <- plot_cpdb(cell_type1 = 'TAM', cell_type2 = 'CD8-EFFECTOR 1|CD8-EFFECTOR 2|CD8-EFFECTOR 3', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE,gene.family='Th2') + 
small_axis(fontsize = 15) + small_grid() + small_guide() + small_legend(fontsize = 5) #
dev.off()


png("TAM-CD8.all.png", units="in", width=10, height=15, res=300)
pl05 <- plot_cpdb(cell_type1 = 'TAM', cell_type2 = 'CD8-EFFECTOR 1|CD8-EFFECTOR 2|CD8-EFFECTOR 3', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE,genes=c("SIRPA","CD47","CD58","CD2","MIF","TNFRSF14","SPP1","CD44","LGALS9","CD44","TIGIT","NECTIN2","CCR1","CCR4","SELL","ICAM1","CD74","CCL3","CCL4","SELL","SELPGL","PVR","CTLA4","PDCD1","BTLA"),col_option = viridis::cividis(50)) + 
small_axis(fontsize = 20) + small_grid() + small_guide() + small_legend(fontsize = 10) #
dev.off()


svg("tam_cd8.1.svg", width=15, height=10)
pl1
dev.off()


svg("tam_cd8.2.svg", width=15, height=10)
pl2
dev.off()

svg("tam_cd8.3.svg", width=15, height=10)
pl3
dev.off()


svg("tam_cd8.4.svg", width=15, height=10)
pl4
dev.off()


svg("tam_cd8.05.svg", width=15, height=10)
plot_cpdb(cell_type1 = 'TAM', cell_type2 = 'CD8-EFFECTOR 1|CD8-EFFECTOR 2|CD8-EFFECTOR 3', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE,genes=c("SIRPA","CD47","CD58","CD2","MIF","TNFRSF14","SPP1","CD44","LGALS9","CD44","TIGIT","NECTIN2","CCR1","CCR4","SELL","ICAM1","CD74","CCL3","CCL4","SELL","SELPGL","PVR","CTLA4","PDCD1","BTLA"),col_option = viridis::cividis(50))
dev.off()



svg("cafs1_cd8.05.svg", width=15, height=10)
plot_cpdb(cell_type1 = 'TGFB-myCAF|detox-IL-iCAF|ecm-myCAF|wound-myCAF', cell_type2 = 'CD8-EFFECTOR 1|CD8-EFFECTOR 2|CD8-EFFECTOR 3', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE,genes=c("PDCD1","LAG3","TIGIT","CTLA4","PDGFRB","THY1","COL1A1","ITGB1","S100A4","CXCL13","CD274","CCL21","MCP1","SDF1","MCSF1","CCL2","CHI3L1","CCL5"),col_option = viridis::cividis(50))
dev.off()

png("cafs-CD8.costimulatory.png", units="in", width=10, height=15, res=300)
plot_cpdb(cell_type1 = 'TGFB-myCAF|detox-IL-iCAF|ecm-myCAF|wound-myCAF', cell_type2 = 'CD8-EFFECTOR 1|CD8-EFFECTOR 2|CD8-EFFECTOR 3', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE,gene.family='costimulatory',col_option = viridis::cividis(50))
dev.off()


png("cafs-CD8.costimulatory.png", units="in", width=10, height=15, res=300)
plot_cpdb(cell_type1 = 'TGFB-myCAF|detox-IL-iCAF|ecm-myCAF|wound-myCAF', cell_type2 = 'CD8-EFFECTOR 1|CD8-EFFECTOR 2|CD8-EFFECTOR 3', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE,gene.family='costimulatory',col_option = viridis::cividis(50))
dev.off()


png("cafs-CD8.coinhibitory.png", units="in", width=10, height=15, res=300)
plot_cpdb(cell_type1 = 'TGFB-myCAF|detox-IL-iCAF|ecm-myCAF|wound-myCAF', cell_type2 = 'CD8-EFFECTOR 1|CD8-EFFECTOR 2|CD8-EFFECTOR 3', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE,gene.family='coinhibitory',col_option = viridis::cividis(50))
dev.off()


svg("cafs1_cd8.svg", width=15, height=10)
plot_cpdb(cell_type1 = 'TGFB-myCAF|detox-IL-iCAF|ecm-myCAF|wound-myCAF', cell_type2 = 'CD8-EFFECTOR 1|CD8-EFFECTOR 2|CD8-EFFECTOR 3', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE,genes=c("PDCD1","LAG3","TIGIT","TNFRSF1B","MIF","CD46","TNF","CD226","HLA","SELL","TNFSF9","CD47","NECTIN2","LGALS9"),col_option = viridis::cividis(50))
dev.off()


plot_cpdb(cell_type1 = 'TGFB-myCAF|detox-IL-iCAF|ecm-myCAF|wound-myCAF', cell_type2 = 'CD8-EFFECTOR 1|CD8-EFFECTOR 2|CD8-EFFECTOR 3', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE,genes=c("PDCD1","LAG3","TIGIT","TNFRSF1B","MIF","CD46","TNF","CD226","HLA","SELL","TNFSF9","CD47","NECTIN2","LGALS9"),col_option = viridis::cividis(50))


###

png("cafs-TAM.costimulatory.png", units="in", width=10, height=15, res=300)
plot_cpdb(cell_type1 = 'TGFB-myCAF|detox-IL-iCAF|ecm-myCAF|wound-myCAF', cell_type2 = 'TAM', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE,gene.family='costimulatory',col_option = viridis::cividis(50))
dev.off()


png("cafs-TAM.coinhibitory.png", units="in", width=10, height=15, res=300)
plot_cpdb(cell_type1 = 'TGFB-myCAF|detox-IL-iCAF|ecm-myCAF|wound-myCAF', cell_type2 = 'TAM', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE,gene.family='coinhibitory',col_option = viridis::cividis(50))
dev.off()

plot_cpdb(cell_type1 = 'TGFB-myCAF|detox-IL-iCAF|ecm-myCAF|wound-myCAF', cell_type2 = 'CD8-EFFECTOR 1|CD8-EFFECTOR 2|CD8-EFFECTOR 3', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE,gene.family='chemokines',col_option = viridis::cividis(50))


plot_cpdb(cell_type1 = 'TGFB-myCAF|detox-IL-iCAF|ecm-myCAF|wound-myCAF', cell_type2 = 'TAM', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE,gene.family='chemokines',col_option = viridis::cividis(50))

chemokines


chemokines

png("cafs-TAM.chemokines.png", units="in", width=10, height=15, res=300)
plot_cpdb(cell_type1 = 'TGFB-myCAF|detox-IL-iCAF|ecm-myCAF|wound-myCAF', cell_type2 = 'TAM', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE,gene.family='chemokines',col_option = viridis::cividis(50))
dev.off()


png("cafs1_tam.svg", units="in", width=10, height=15, res=300)
plot_cpdb(cell_type1 = 'TGFB-myCAF|detox-IL-iCAF|ecm-myCAF|wound-myCAF', cell_type2 = 'TAM', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE,genes=c("FOLRA","SPP1","MT1G","CD163","CD40LD","CD40","CD28","SIRPA","CD86","NECTIN1","IFNG","IL15","CSF1R","CXCL12","CXCR6","PDCD1","TNFSF13","TNF","LGALS9","CD46","NOTCH","HLA","CD55","CD99","CD74","PGRMC2","VS1R"),col_option = viridis::cividis(50))
dev.off()




#Cytokines

"CCL2","CXCL2","CXCL12","CCL1","C3","LGALS9","CCL2","CSF1R","CF1","CXCL9","CXCL8","CXCL16","CXCL10","CXCL1","CCL5","CCL4L2","CCL22","CCL21","CCL17","CCL13","CCL11","DPPR","ACKR1","CXCR6","CXCR4","CXCR3","ACKR3","CCR5","CCR1"

#growth factors 

"VEGFB","VEGFA","PGF","PDGFC","IGF1","IGF2","HGF","FGF9","FGF7","FGF2","FGF10","NRP1","FLT1","NRP2","ADRB2","KDR","FLT4","IGF2R","IGF1R","MET","CD44","FGFR4","FGFR3","FGFR2","FGFR1","CD44"



png("tam_cd8.cytokines.png", , units="in", width=7, height=10, res=300)
plot_cpdb(cell_type1 = 'TAM', cell_type2 = 'CD8-EFFECTOR 1|CD8-EFFECTOR 2|CD8-EFFECTOR 3', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE,genes=c("CCL2","CXCL2","CXCL12","CCL1","C3","LGALS9","CCL2","CSF1R","CF1","CXCL9","CXCL8","CXCL16","CXCL10","CXCL1","CCL5","CCL4L2","CCL22","CCL21","CCL17","CCL13","CCL11","DPPR","ACKR1","CXCR6","CXCR4","CXCR3","ACKR3","CCR5","CCR1"),col_option = viridis::cividis(50))
dev.off()



png("tam_cd8.growthFactors.png", , units="in", width=6, height=10, res=300)
plot_cpdb(cell_type1 = 'TAM', cell_type2 = 'CD8-EFFECTOR 1|CD8-EFFECTOR 2|CD8-EFFECTOR 3', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE,genes=c("VEGFB","VEGFA","PGF","PDGFC","IGF1","IGF2","HGF","FGF9","FGF7","FGF2","FGF10","NRP1","FLT1","NRP2","ADRB2","KDR","FLT4","IGF2R","IGF1R","MET","CD44","FGFR4","FGFR3","FGFR2","FGFR1","CD44"),col_option = viridis::cividis(50))
dev.off()



png("cafs_cd8_cytokines.png", units="in", width=10, height=10, res=300)
plot_cpdb(cell_type1 = 'TGFB-myCAF|detox-IL-iCAF|ecm-myCAF|wound-myCAF', cell_type2 = 'CD8-EFFECTOR 1|CD8-EFFECTOR 2|CD8-EFFECTOR 3', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE,genes=c("CCL2","CXCL2","CXCL12","CCL1","C3","LGALS9","CCL2","CSF1R","CF1","CXCL9","CXCL8","CXCL16","CXCL10","CXCL1","CCL5","CCL4L2","CCL22","CCL21","CCL17","CCL13","CCL11","DPPR","ACKR1","CXCR6","CXCR4","CXCR3","ACKR3","CCR5","CCR1"
),col_option = viridis::cividis(50))
dev.off()


png("cafs_cd8_growthFactors.png", units="in", width=10, height=10, res=300)
plot_cpdb(cell_type1 = 'TGFB-myCAF|detox-IL-iCAF|ecm-myCAF|wound-myCAF', cell_type2 = 'CD8-EFFECTOR 1|CD8-EFFECTOR 2|CD8-EFFECTOR 3', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE,genes=c("VEGFB","VEGFA","PGF","PDGFC","IGF1","IGF2","HGF","FGF9","FGF7","FGF2","FGF10","NRP1","FLT1","NRP2","ADRB2","KDR","FLT4","IGF2R","IGF1R","MET","CD44","FGFR4","FGFR3","FGFR2","FGFR1","CD44"),col_option = viridis::cividis(50))
dev.off()



png("tam_cafs.cytokines.png", , units="in", width=10, height=10, res=300)
plot_cpdb(cell_type1 = 'TAM', cell_type2 = 'TGFB-myCAF|detox-IL-iCAF|ecm-myCAF|wound-myCAF', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE,genes=c("CCL2","CXCL2","CXCL12","CCL1","C3","LGALS9","CCL2","CSF1R","CF1","CXCL9","CXCL8","CXCL16","CXCL10","CXCL1","CCL5","CCL4L2","CCL22","CCL21","CCL17","CCL13","CCL11","DPPR","ACKR1","CXCR6","CXCR4","CXCR3","ACKR3","CCR5","CCR1"),col_option = viridis::cividis(50))
dev.off()


png("tam_cafs.growthFactors.png", , units="in", width=10, height=10, res=300)
plot_cpdb(cell_type1 = 'TAM', cell_type2 = 'TGFB-myCAF|detox-IL-iCAF|ecm-myCAF|wound-myCAF', scdata = sce,
  idents = 'comb_clusters', # column name where the cell ids are located in the metadata
  means = means, pvals = pvals,keep_significant_only = TRUE,genes=c("VEGFB","VEGFA","PGF","PDGFC","IGF1","IGF2","HGF","FGF9","FGF7","FGF2","FGF10","NRP1","FLT1","NRP2","ADRB2","KDR","FLT4","IGF2R","IGF1R","MET","CD44","FGFR4","FGFR3","FGFR2","FGFR1","CD44"),col_option = viridis::cividis(50))
dev.off()


####
library(ktplots)
cell_type1 = 'TAM', cell_type2 = 'CD8-EFFECTOR 1|CD8-EFFECTOR 2|CD8-EFFECTOR 3'

#sce <- Seurat::as.SingleCellExperiment(kidneyimmune)

p <- plot_cpdb2(cell_type1 = 'TAM', cell_type2 = 'CD8-EFFECTOR 1|CD8-EFFECTOR 2|CD8-EFFECTOR 3',
    scdata = sce,
    idents = 'comb_clusters', # column name where the cell ids are located in the metadata
    means = means,
    pvals = pvals,
    deconvoluted = decon2, # new options from here on specific to plot_cpdb2
    desiredInteractions = list(
        c('TAM'),
        c('CD8-EFFECTOR 1', 'CD8-EFFECTOR 2' ,'CD8-EFFECTOR 3')),
    interaction_grouping = interaction_annotation,
    edge_group_colors = c(
        "Activating" = "#e15759",
        "Chemotaxis" = "#59a14f",
        "Inhibitory" = "#4e79a7",
        "Intracellular trafficking" = "#9c755f",
        "DC_development" = "#B07aa1",
        "Unknown" = "#e7e7e7"
        ),
    node_group_colors = c(
        "TAM" = "red",
        "CD8-EFFECTOR 1" = "blue"),
    keep_significant_only = TRUE,
    standard_scale = TRUE,
    remove_self = TRUE
    )
p







#### Pelka Analysis

#home/akhaliq/data/pelka

library(Seurat)
library(knitr)

pelka.h5<- Seurat::Read10X_h5(filename = "GSE178341_crc10x_full_c295v4_submit.h5", use.names =T)

pelka.obj <- CreateSeuratObject(pelka.h5, project = "pelka")
pelka.meta <- read.csv("GSE178341_crc10x_full_c295v4_submit_metatables.csv",header=TRUE,sep=',', row.names=1)
pelka.alldata<- AddMetaData(pelka.obj, pelka.meta)
pelka.clust <- read.csv("GSE178341_crc10x_full_c295v4_submit_cluster.csv",header=TRUE,sep=',', row.names=1)
pelka.alldata<- AddMetaData(pelka.alldata, pelka.clust)

head(pelka.alldata@meta.data)
saveRDS(pelka.alldata,"pelka.alldata.rds")

# since its Seurat V 4 the subset function is bit different

pelka.epi <- subset(x = pelka.alldata, subset = clMidwayPr == "Epi")


kable(table(pelka.alldata@meta.data$clMidwayPr))
|Var1         |   Freq|
|:------------|------:|
|B            |  25660|
|DC           |   5549|
|Endo         |   7520|
|Epi          | 168295|
|Fibro        |   5231|
|Granulo      |   2043|
|ILC          |    832|
|Macro        |  20280|
|Mast         |   3834|
|Mono         |  14242|
|NK           |   3924|
|Peri         |   1525|
|Plasma       |  37809|
|Schwann      |    281|
|SmoothMuscle |    881|
|TCD4         |  34598|
|TCD8         |  23486|
|Tgd          |   9383|
|TZBTB16      |   4742|



> kable(table(pelka.alldata@meta.data$SPECIMEN_TYPE))


|Var1 |   Freq|
|:----|------:|
|N    | 112864|
|T    | 257251|

> kable(table(pelka.epi@meta.data$SPECIMEN_TYPE))


|Var1 |   Freq|
|:----|------:|
|N    |  55079|
|T    | 113216|



pelka.epi.t <- subset(x = pelka.epi, subset = SPECIMEN_TYPE == "T")


####
Convert seurat obj to Anndata 

####
# https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html
library(Seurat)
library(SeuratData)
library(SeuratDisk)

SaveH5Seurat(epi.tumor, filename = "epi.tumor.h5Seurat")

Convert("epi.tumor.h5Seurat", dest = "h5ad")



#####
Cellrank
#####


Conda activate seurat

# in python terminal
import sys

if "google.colab" in sys.modules:
    !pip install -q git+https://github.com/theislab/cellrank@dev
    !pip install python-igraph


##SlingShot


suppressPackageStartupMessages({
    # Single cell libraries
    library(Seurat)
    library(scran)
    library(scater)

    # Plotting
    library(rafalib)
    library(cowplot)
    library(rgl)
    library(plotly)
    options(rgl.printRglwidget = TRUE)

    # Sparse Matrix manipulation tools
    library(Matrix)
    library(sparseMatrixStats) # this will be installed only in R v 4

    # Trajectory analysis
    library(slingshot)
    library(tradeSeq)

    # Dimensionality reduction library(destiny)
    library(fastICA)
})

# Define some color palette
pal <- c((scales::hue_pal())(8), RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8,
    "Set2"))
set.seed(1)
pal <- rep(sample(pal, length(pal)), 200)




# Add graph to the base R graphics plot
draw_graph <- function(layout, graph, lwd = 0.2, col = "grey") {
    res <- rep(x = 1:(length(graph@p) - 1), times = (graph@p[-1] - graph@p[-length(graph@p)]))
    segments(x0 = layout[graph@i + 1, 1], x1 = layout[res, 1], y0 = layout[graph@i +
        1, 2], y1 = layout[res, 2], lwd = lwd, col = col)
}


obj = belgian.tumor

vars <- c("Condition","Cell_type","Cell_subtype","MSI.Status","Location")


pl <- list()


png("png", , units="in", width=10, height=10, res=300)
for (i in vars) {
    pl[[i]] <- DimPlot(obj, group.by = i, label = T) + theme_void() + NoLegend()
}
plot_grid(plotlist = pl)
dev.off()

vars <- c("CMS1","CMS2","CMS3","CMS4")


pl <- list()

# SEURAT
pl <- list(DimPlot(obj, group.by = "Cell_subtype", label = T) + theme_void() + NoLegend())
for (i in vars) {
    pl[[i]] <- FeaturePlot(obj, features = i, order = T) + theme_void() + NoLegend()
}
plot_grid(plotlist = pl)



### Slim down the seurat obj

obj <- DietSeurat(obj,count=TRUE,data=TRUE)

###################### 
# SlingShot Analysis Pipeline
######################

ref: https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_03_integration.html


suppressPackageStartupMessages({
    # Single cell libraries
    library(Seurat)
    library(scran)
    library(scater)

    # Plotting
    library(rafalib)
    library(cowplot)
    library(rgl)
    library(plotly)
    options(rgl.printRglwidget = TRUE)

    # Sparse Matrix manipulation tools
    library(Matrix)
    library(sparseMatrixStats)

    # Trajectory analysis
    library(slingshot)
    library(tradeSeq)

    # Dimensionality reduction library(destiny)
    library(fastICA)
})

# Define some color palette
pal <- c((scales::hue_pal())(8), RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8,
    "Set2"))
set.seed(1)
pal <- rep(sample(pal, length(pal)), 200)


# Add graph to the base R graphics plot
draw_graph <- function(layout, graph, lwd = 0.2, col = "grey") {
    res <- rep(x = 1:(length(graph@p) - 1), times = (graph@p[-1] - graph@p[-length(graph@p)]))
    segments(x0 = layout[graph@i + 1, 1], x1 = layout[res, 1], y0 = layout[graph@i +
        1, 2], y1 = layout[res, 2], lwd = lwd, col = col)
}

#### Usual Seurat clustering pipeline

obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE,assay="RNA")
obj <- ScaleData(obj, verbose = FALSE)

obj <- RunPCA(obj, verbose = FALSE, features = VariableFeatures(object = obj))
obj <- RunTSNE(object = obj, dims.use = 1:30)
obj <- RunUMAP(obj, dims = 1:30) 
obj <- RunUMAP(obj, dims = 1:30,n.components = 3,reduction.name = "umap3d") 
obj <- FindNeighbors(obj, dims = 1:30, force.recalc = T)

obj <- FindClusters(obj, resolution = 0.8, algorithm = 1, force.recalc=T)

obj$seurat_clusters <- as.factor(as.numeric(as.character(obj$seurat_clusters)) + 1)
Idents(obj) <- "seurat_clusters"



obj.harmony <- RunHarmony(obj, group.by.vars = "orig.ident", reduction = "pca", dims.use = 1:50, assay.use = "RNA")

###########


table(obj$seurat_clusters)

NORM_COUNTS <- obj@assays$RNA@data
UMAP2 <- obj@reductions$umap@cell.embeddings
UMAP3 <- obj@reductions$umap3d@cell.embeddings
HARMONY <- obj@reductions$harmony@cell.embeddings
PCA <- obj@reductions$pca@cell.embeddings
PCA_loadings <- obj@reductions$pca@feature.loadings
clustering <- factor(obj$seurat_clusters)
KNN <- obj@graphs$CCA_nn # check with this before running


mm <- sparse.model.matrix(~0 + factor(clustering))
colnames(mm) <- levels(factor(clustering))
centroids3d <- as.matrix(t(t(UMAP3) %*% mm)/Matrix::colSums(mm))
centroids2d <- as.matrix(t(t(UMAP2) %*% mm)/Matrix::colSums(mm))



rgl::open3d()

points3d(x = UMAP3[, 1], y = UMAP3[, 2], z = UMAP3[, 3], col = pal[factor(clustering)])
text3d((centroids3d[, 1:3]), texts = rownames(centroids3d), cex = 1)


# calculate the distance from every cell to its neighbors
expected_U2d <- t(t(UMAP2) %*% KNN)/colSums2(KNN)
d <- rowSums((expected_U2d - UMAP2)^2)^(1/2)

# Define a distance cutoff
hist(d, breaks = 400)
cutoff <- mean(d) + 5 * sd(d)
abline(v = (cutoff), col = "red", xpd = F)

cutoff



to_keep <- (d < cutoff)

mypar()
plot(UMAP2, type = "n")
draw_graph(layout = UMAP2, graph = KNN)
points(UMAP2, cex = ifelse(!to_keep, 1, 0.3), lwd = ifelse(!to_keep, 2, 0), bg = pal[clustering],
    pch = 21)
text(centroids2d, labels = rownames(centroids2d), cex = 1, font = 2)


new_UMAP2 <- UMAP2
res <- as.matrix(t(t(UMAP2) %*% KNN)/colSums2(KNN))
new_UMAP2[!to_keep, ] <- res[!to_keep, ]

new_centroids2d <- as.matrix(t(t(new_UMAP2) %*% mm)/Matrix::colSums(mm))


# Check the UMAP in 2D
mypar()
plot(new_UMAP2, type = "n")
draw_graph(layout = new_UMAP2, graph = KNN)
points(new_UMAP2, cex = ifelse(!to_keep, 1, 0.3), lwd = ifelse(!to_keep, 2, 0), bg = pal[clustering],
    pch = 21)
text(new_centroids2d, labels = rownames(new_centroids2d), cex = 1, font = 2)




#Trajectory inference with Slingshot


# Computing ICA
ICA_object <- fastICA::fastICA(X = HARMONY, n.comp = 20, method = "C", row.norm = T)


ICA <- ICA_object$S
colnames(ICA) <- paste0("ICA_", 1:ncol(ICA))


png("diffusion_maps_lee_crc.png",units="in" , width=10, height=10,res=300)
mypar(3, 3)
plot(new_UMAP2, pch = 16, col = pal[clustering])
text(new_centroids2d, labels = rownames(new_centroids2d), cex = 1, font = 2)

for (i in 1:8) {
    cc <- t(t(ICA[, c(i * 2 - 1, i * 2)]) %*% mm)/Matrix::colSums(mm)
    plot(ICA[, c(i * 2 - 1, i * 2)], pch = 16, col = pal[clustering])
    text(cc[, 1], cc[, 2], labels = rownames(cc), cex = 1, font = 2)
}

dev.off()



# Run Slingshot on UMAP3d
set.seed(1)
lineages <- as.SlingshotDataSet(getLineages(data = new_UMAP2, clusterLabels = clustering))

# Change the reduction (FOR VISUALISATION ONLY, in case you use another
# dimension for calculations)
lineages


lineages@reducedDim <- new_UMAP2


png("lineages_lee_crc.png",units="in" , width=10, height=7,res=300)
# Plot the lineages
mypar(1, 2)
plot(new_UMAP2, col = pal[clustering], cex = 0.5, pch = 16)
lines(lineages, lwd = 2, col = "black", cex = 3)
text(new_centroids2d, labels = rownames(new_centroids2d), cex = 0.8, font = 2, col = "white")

# Check the UMAP in 2D
plot(new_UMAP2, type = "n")
draw_graph(layout = new_UMAP2, graph = KNN)
points(new_UMAP2, cex = ifelse(!to_keep, 1, 0.3), lwd = ifelse(!to_keep, 2, 0), bg = pal[clustering],
    pch = 21)
text(new_centroids2d, labels = rownames(new_centroids2d), cex = 0.8, font = 2)
dev.off()

png("dim_LEE_crc.png",units="in" , width=7, height=7,res=300)
DimPlot(obj, group.by = "Cell_subtype", label = F,pt.size=1,cols = c("#dcc134", "#008000", "#37c8ab","#ff5555"))
dev.off()




svg("lineages_lee_crc.svg", width=10, height=10)
# Plot the lineages
mypar(1, 2)
plot(new_UMAP2, col = pal[clustering], cex = 0.5, pch = 16)
lines(lineages, lwd = 2, col = "black", cex = 3)
text(new_centroids2d, labels = rownames(new_centroids2d), cex = 0.8, font = 2, col = "white")

# Check the UMAP in 2D
plot(new_UMAP2, type = "n")
draw_graph(layout = new_UMAP2, graph = KNN)
points(new_UMAP2, cex = ifelse(!to_keep, 1, 0.3), lwd = ifelse(!to_keep, 2, 0), bg = pal[clustering],
    pch = 21)
text(new_centroids2d, labels = rownames(new_centroids2d), cex = 0.8, font = 2)
dev.off()


svg("dim_leecrc.svg", width=7, height=7)
DimPlot(obj, group.by = "Cell_subtype", label = F,pt.size=1,cols = c("#dcc134", "#008000", "#37c8ab","#ff5555"))
dev.off()


##################################################################


###
Intercellular
###


/Users/akhaliq/Desktop/genome_biology/cellphonedb_analysis/cpdb_intercellar_results



#############################
Python script for running the samples effectively


import os
a="/home/akhaliq/data/sr_new/"
#os.mkdir(a)
for i in open("area_new.txt"):
  j=i.strip().split("\t")
  #b=a+sample+"/"
  sample=j[0]
  area=j[1]
  b=a+sample+"/"
  os.mkdir(b)
  os.chdir(b)
  os.system('/home/akhaliq/data/spaceranger/spaceranger-1.3.1/spaceranger count --id="'+sample+'" --description="PDAC Spatial" --transcriptome=/home/akhaliq/data/spatial_anal_AK/refdata-gex-GRCh38-2020-A --fastqs=/home/akhaliq/data/spatialdata_new/'+sample+' --localcores=32 --localmem=128 --image=/home/akhaliq/data/spatialdata_new/'+sample+'/spatial/detected_tissue_image.jpg --slide=V11L12-005 --area='+area+' --probe-set=/home/akhaliq/data/spatial_anal_AK/Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv  --reorient-images=')
  #os.system('/home/akhaliq/data/spaceranger/spaceranger-1.3.1/spaceranger count --id="'+sample+'" --description="PDAC Spatial" --transcriptome=/home/akhaliq/data/spatial_anal_AK/refdata-gex-GRCh38-2020-A --fastqs=/home/akhaliq/data/spatialdata_new/'+sample+' --localcores=32 --localmem=128 --image=/home/akhaliq/data/spatialdata_new/'+sample+'/spatial/detected_tissue_image.jpg --slide=V11L12-005 --area='+area+' --probe-set=/home/akhaliq/data/spatial_anal_AK/Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv  --reorient-images=')
  os.chdir(a)
  os.getcwd()

---
NOTE:save as .py
Note: area_new.txt looks as follows

Sample_ID area
M016  A1
M017  B1
M018  C1
M019  D1
M020  A1
M021  B1
M022  C1
M023  D1
M024  A1
M025  B1
M026  C1

____


############
HPC access

ssh bigred3.uits.iu.edu


# to transfer files...
scp addu.jpeg bigred3.uits.iu.edu:/N/u/akhaliq/BigRed3/test

# for copying whole directory

scp -r sr_new bigred3.uits.iu.edu:/N/u/akhaliq/BigRed3

# carbonate

scp /Users/akhaliq/Desktop/pdac_mets/CLR__20220520_161157.tiff carbonate.uits.iu.edu:/N/slate/akhaliq
/N/slate/akhaliq

ssh carbonate.uits.iu.edu carbonate.uits.iu.edu:/N/slate/akhaliq

## for copying files form slate to local hard drive
scp -r carbonate.uits.iu.edu carbonate.uits.iu.edu:/N/slate/akhaliq/sc_ac /Volumes/Seagate\ Bac/bk_GCP

#### SLATE 

/N/slate/akhaliq

scp -r crc_dataset carbonate.uits.iu.edu:/N/slate/akhaliq

# to check the modules available
module avail

# conda activation
module load anaconda/python3.8/2020.07
conda activate seurat

# to get a tree on files
#conda install -c conda-forge tree

tree -d
#########
# Convert h5ad file to seurat object
#########
conda activate seurat

library(Seurat)
library(SeuratDisk)

Convert("VUMC_HTAN_DIS_EPI_V2.h5ad",dest = "h5seurat", overwrite = TRUE)

htan_62 <- LoadH5Seurat("VUMC_HTAN_DIS_EPI_V2.h5seurat")


##########

#Run in Slrum

#Submit a Job
sbatch s18.c9.script

#to check the job running in username
squeue -u akhaliq  -p dl -t PENDING
#killjob Scancel <jobid>

# scp from carbonate to local 

# in a new local terminal

scp -r carbonate.uits.iu.edu carbonate.uits.iu.edu:/N/slate/akhaliq/sc_analysis/pdac_mets/spatials_mets/analysis/spatial_counts2.png /Users/akhaliq/Desktop/pdac_mets/new


########
Cytospace
########


module load anaconda/python3.8/2020.07
conda activate cytospace

#

sbatch my_job.script

/Volumes/Seagate Bac/bkup_mac_0430/gb/Additional file 1/Fig. S5.pdf

#



SpatialFeaturePlot(st, features = c("EPCAM","KRT8","KRT18"), pt.size.factor = 2, alpha = c(0.1, 1),images="s18_a1_Liver_Mets.3")


####################
Getting Count martix for huge counts.txt
This is a workaround using Python. I managed to export a matrix that was too big in R.
########

ref= https://stackoverflow.com/questions/54712621/writing-a-very-large-sparse-matrix-to-file-in-r


Export data in R as sparse matrix:

library(Matrix)
counts <- epi@assays$RNA@counts

write(colnames(counts), file = "colnames.txt")
write(rownames(counts), file = "rownames.txt")
writeMM(counts, file = "sparsematrix.txt")

#Read then convert in Python:

from scipy import sparse, io
import pandas as pd
import numpy as np

sparsematrix = io.mmread('sparsematrix.txt')

m_dense = sparsematrix.toarray()

var_names = np.genfromtxt('rownames.txt', dtype=str)
col_names = np.genfromtxt('colnames.txt', dtype=str)

# Export to txt:
df = pd.DataFrame(m_dense, columns=col_names, index=var_names)
df.to_csv('export_sparsematrix.txt', sep='\t', header=True, index=True, index_label='Genes')

######

#Count number of cells expressing a gene
Ref: 
sum(GetAssayData(object = my_object, slot = "data")[my_gene,my_cells]>0)
sum(GetAssayData(object = my_object, slot = "data")[my_gene,my_cells]>0)/length(my_cells)


#####
Subsetting the cluster-sampl
######

At times subset will not work so  first set the ident to which wver thing you want to subst and then subset, for example

Idents(peng) <- "seurat_clusters"

peng.sub <- subset(peng, idents = c("Epithelial cells","Stem cells","Fibroblasts","Endocrine cells"))
kable(table(peng.sub$seurat_clusters))

### All R packages are installed in the following Locations

‘/Users/akhaliq/Library/R/x86_64/4.2/library’

919902922617

zabiadv79@gmail.com
+919972120340
919108956799


#

## replacing teh orig.idents 
# This code replaces the patterns in old_patterns with the corresponding ones in new_patterns using the gsub() function. Note that the match() function is used inside the replacement argument to find the index of the matching pattern in old_patterns, and use it to get the corresponding pattern in new_patterns. The resulting orig.ident column is updated in samp@meta.data.

# get the old and new patterns
old_patterns <- c("IU_PDA_NH3", "IU_PDA_HM4", "IU_PDA_HM5", "IU_PDA_HM6", "IU_PDA_LNM7", "IU_PDA_LNM8", "IU_PDA_LNM9", "IU_PDA_LNM11", "IU_PDA_NP11", "IU_PDA_NP12")
new_patterns <- c("IU_PDA_NH2", "IU_PDA_HM4", "IU_PDA_HM5", "IU_PDA_HM6", "IU_PDA_LNM6", "IU_PDA_LNM8", "IU_PDA_T9", "IU_PDA_LNM10", "IU_PDA_NP10", "IU_PDA_LNM12")

# replace the old patterns with the new ones
samp@meta.data$orig.ident <- gsub(
  pattern = paste0("(", paste(old_patterns, collapse = "|"), ")"),
  replacement = function(x) new_patterns[match(x, old_patterns)],
  x = samp@meta.data$orig.ident
)

###Changing the image name in seurat spatial data

old_names <- c("IU_PDA_T1", "IU_PDA_NP2", "IU_PDA_T2", "IU_PDA_NH3", "IU_PDA_HM3", "IU_PDA_HM4", "IU_PDA_T4", "IU_PDA_HM5", "IU_PDA_T5", "IU_PDA_HM6", "IU_PDA_HM7", "IU_PDA_LNM7", "IU_PDA_T7", "IU_PDA_LNM8", "IU_PDA_HM9", "IU_PDA_LNM9", "IU_PDA_T9", "IU_PDA_T10", "IU_PDA_HM10", "IU_PDA_HM11", "IU_PDA_LNM11", "IU_PDA_NP11", "IU_PDA_T11", "IU_PDA_HM12", "IU_PDA_T12", "IU_PDA_NP12", "IU_PDA_HM13", "IU_PDA_LNM13", "IU_PDA_T13", "IU_PDA_HM14")

new_names <- c("IU_PDA_T1", "IU_PDA_NP2", "IU_PDA_T2", "IU_PDA_NH2", "IU_PDA_HM2", "IU_PDA_HM3", "IU_PDA_T3", "IU_PDA_HM4", "IU_PDA_T4", "IU_PDA_HM5", "IU_PDA_HM6", "IU_PDA_LNM6", "IU_PDA_T6", "IU_PDA_LNM7", "IU_PDA_HM8", "IU_PDA_LNM8", "IU_PDA_T8", "IU_PDA_T9", "IU_PDA_HM9", "IU_PDA_HM10", "IU_PDA_LNM10", "IU_PDA_NP10", "IU_PDA_T10", "IU_PDA_HM11", "IU_PDA_T11", "IU_PDA_NP11", "IU_PDA_HM12", "IU_PDA_LNM12", "IU_PDA_T12", "IU_PDA_HM13")

names(samp@images) <- new_names

####

# Changing the name of the images in St object 

# Define the old and new image name mapping
image_mapping <- c(
  "IU_PDA_T1" = "IU_PDA_HM9_t",
  "IU_PDA_NP2" = "IU_PDA_HM10_t",
  "IU_PDA_T2" = "IU_PDA_HM11_t",
  "IU_PDA_NH2" = "IU_PDA_HM12_t",
  "IU_PDA_HM2" = "IU_PDA_HM13_t",
  "IU_PDA_HM3" = "IU_PDA_HM2_t",
  "IU_PDA_T3" = "IU_PDA_HM3_t",
  "IU_PDA_HM4" = "IU_PDA_HM4_t",
  "IU_PDA_T4" = "IU_PDA_HM5_t",
  "IU_PDA_HM5" = "IU_PDA_HM6_t",
  "IU_PDA_HM6" = "IU_PDA_HM8_t",
  "IU_PDA_LNM6" = "IU_PDA_LNM10_t",
  "IU_PDA_T6" = "IU_PDA_LNM12_t",
  "IU_PDA_LNM7" = "IU_PDA_LNM6_t",
  "IU_PDA_HM8" = "IU_PDA_LNM7_t",
  "IU_PDA_LNM8" = "IU_PDA_LNM8_t",
  "IU_PDA_T8" = "IU_PDA_NH2_t",
  "IU_PDA_T9" = "IU_PDA_NP10_t",
  "IU_PDA_HM9" = "IU_PDA_NP11_t",
  "IU_PDA_HM10" = "IU_PDA_NP2_t",
  "IU_PDA_LNM10" = "IU_PDA_T1_t",
  "IU_PDA_NP10" = "IU_PDA_T9_t",
  "IU_PDA_T10" = "IU_PDA_T10_t",
  "IU_PDA_HM11" = "IU_PDA_T11_t",
  "IU_PDA_T11" = "IU_PDA_T12_t",
  "IU_PDA_NP11" = "IU_PDA_T2_t",
  "IU_PDA_HM12" = "IU_PDA_T3_t",
  "IU_PDA_LNM12" = "IU_PDA_T4_t",
  "IU_PDA_T12" = "IU_PDA_T6_t",
  "IU_PDA_HM13" = "IU_PDA_T8_t"
)

# Rename the image names in pdac_allres
names(pdac_allres@images) <- image_mapping[names(pdac_allres@images)]

# now remove the _t from the image to the following

# Define the old and new image name mapping
image_mapping <- c(
  "IU_PDA_HM9_t" = "IU_PDA_HM9",
  "IU_PDA_HM10_t" = "IU_PDA_HM10",
  "IU_PDA_HM11_t" = "IU_PDA_HM11",
  "IU_PDA_HM12_t" = "IU_PDA_HM12",
  "IU_PDA_HM13_t" = "IU_PDA_HM13",
  "IU_PDA_HM2_t" = "IU_PDA_HM2",
  "IU_PDA_HM3_t" = "IU_PDA_HM3",
  "IU_PDA_HM4_t" = "IU_PDA_HM4",
  "IU_PDA_HM5_t" = "IU_PDA_HM5",
  "IU_PDA_HM6_t" = "IU_PDA_HM6",
  "IU_PDA_HM8_t" = "IU_PDA_HM8",
  "IU_PDA_LNM10_t" = "IU_PDA_LNM10",
  "IU_PDA_LNM12_t" = "IU_PDA_LNM12",
  "IU_PDA_LNM6_t" = "IU_PDA_LNM6",
  "IU_PDA_LNM7_t" = "IU_PDA_LNM7",
  "IU_PDA_LNM8_t" = "IU_PDA_LNM8",
  "IU_PDA_NH2_t" = "IU_PDA_NH2",
  "IU_PDA_NP10_t" = "IU_PDA_NP10",
  "IU_PDA_NP11_t" = "IU_PDA_NP11",
  "IU_PDA_NP2_t" = "IU_PDA_NP2",
  "IU_PDA_T1_t" = "IU_PDA_T1",
  "IU_PDA_T9_t" = "IU_PDA_T9",
  "IU_PDA_T10_t" = "IU_PDA_T10",
  "IU_PDA_T11_t" = "IU_PDA_T11",
  "IU_PDA_T12_t" = "IU_PDA_T12",
  "IU_PDA_T2_t" = "IU_PDA_T2",
  "IU_PDA_T3_t" = "IU_PDA_T3",
  "IU_PDA_T4_t" = "IU_PDA_T4",
  "IU_PDA_T6_t" = "IU_PDA_T6",
  "IU_PDA_T8_t" = "IU_PDA_T8"
)

# Rename the image names in pdac_allres
names(pdac_allres@images) <- image_mapping[names(pdac_allres@images)]


####

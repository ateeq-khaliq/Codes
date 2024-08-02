

new <- NormalizeData(new, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
new <- FindVariableFeatures(new, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
new <- ScaleData(new, verbose = FALSE)

new <- RunPCA(new, verbose = FALSE, features = VariableFeatures(object = new),npcs=60)
new <- RunTSNE(object = new, dims.use = 1:60)
new <- RunUMAP(new, dims = 1:60)
new <- FindNeighbors(new, dims = 1:60, force.recalc = T)

#new_allres <- FindClusters(new, resolution = c(0,0.2,0.5,0.8,1.5,2,2.5,3,3.5,4) , algorithm = 1, force.recalc=T)

new <- FindClusters(new, resolution = 2, algorithm = 1, force.recalc=T)


#new <- FindClusters(new, resolution = 0.8, algorithm = 1, force.recalc=T)

new$seurat_clusters <- as.factor(as.numeric(as.character(new$seurat_clusters)) + 1)
Idents(new) <- "seurat_clusters"

table(new$seurat_clusters)
png("UMAP_Orig_ident_new_2_celltypes_All_merged.png", units="in", width=10, height=10, res=300)
DimPlot(new, reduction = "umap", label=T, group.by="new_celltype")
dev.off()

markers_genes_new <- FindAllMarkers(new, logfc.threshold = 0.2, test.use = "wilcox", only.pos = TRUE,assay = "RNA")
write.table(markers_genes_new, file="marker_genes_new_2.txt", sep="\t",append=F, row.names = F)
cbind(table(markers_genes_new$cluster))

saveRDS(new,"new.rds")

## Merge 3 objects

new <- merge(new,c(anu,dip))

### ProjectTils
library(ProjecTILs)
Idents(new)<- "new_celltype"
tcells <- subset(new,ident="T cells")

ref <- load.reference.map()
query.projected <- Run.ProjecTILs(tcells, ref = ref)

png("ProjectTILs_Projections.png", units="in", width=5, height=5, res=300)
plot.projection(ref, query.projected, linesize = 0.5, pointsize = 0.5)
dev.off()


png("ProjectTILs_composition.png", units="in", width=10, height=10, res=300)
plot.statepred.composition(ref, query.projected, metric = "Percent")
dev.off()


genes4radar = c("Foxp3", "Cd4", "Cd8a", "Tcf7", "Ccr7", "Sell", "Gzmb", "Gzmk", "Pdcd1",
    "Havcr2", "Tox", "Mki67")

png("ProjectTILs_radar.png", units="in", width=20, height=10, res=300)
plot.states.radar(ref, query = query.projected, genes4radar = genes4radar, min.cells = 20)
dev.off()

png("UMAP_ProjectIL_Tcells.png", units="in", width=8, height=5, res=300)
DimPlot(query.projected, reduction = "umap", label=T, group.by="functional.cluster")
dev.off()

png("Barplot_ProjectIL_Tcells_treatment.png", units="in", width=10, height=10, res=300)
dittoBarPlot(query.projected, "functional.cluster", group.by = "treatment_group")
dev.off()


#annotate
#All markers
python /Users/akhaliq/Desktop/mouse/annotations/annotate_markers_Script/annotate.py marker_genes_new_2.txt marker_genes_new_alltype_2.xlsx /Users/akhaliq/Desktop/mouse/annotations/annotate_markers_Script/markers/master_markerListMouse.txt

#basicMarkers

python /Users/akhaliq/Desktop/mouse/annotations/annotate_markers_Script/annotate.py marker_genes_new_2.txt marker_genes_new_BASICtype_2.xlsx /Users/akhaliq/Desktop/mouse/annotations/annotate_markers_Script/markers/basic_markers.txt


python /Users/akhaliq/Desktop/mouse/annotations/annotate_markers_Script/annotate.py /Users/akhaliq/Downloads/markers_genes_cc10.1.txt /Users/akhaliq/Downloads/markers_genes_cc10.1.xlsx /Users/akhaliq/Desktop/MarkerGenesListCode/markers_master.txt

top50 <- markers_genes_new %>% group_by(cluster) %>% top_n(50, avg_log2FC)


FeaturePlot(new, "Anxa2", cols.use = c("lightgrey", "blue") )

FeaturePlot(new, features = c("Cd4","Cd8a"), max.cutoff = 3, cols = c("grey", "red"),label.size = 4)

DimPlot(new, reduction = "umap", label=F, group.by="copykat_pred")


#rename Clusters

current.cluster.ids <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40", "41", "42", "43", "44", "45", "46", "47", "48", "49", "50")



new.cluster.ids <- c("Fibroblasts", "Fibroblasts", "Fibroblasts", "T cells", "Myeloids", "Myeloids", "Myeloids", "Myeloids", "Myeloids", "B Cells", "T cells", "Myeloids", "Fibroblasts", "Fibroblasts", "Myeloids", "Myeloids", "Myeloids", "Myeloids", "Fibroblasts", "Myeloids", "T cells", "Myeloids", "Myeloids", "Fibroblasts", "Fibroblasts", "Fibroblasts", "Myeloids", "Fibroblasts", "Myeloids", "T cells", "Myeloids", "Endothelial Cells", "Myeloids", "Myeloids", "B Cells", "Fibroblasts", "DC", "T cells", "Fibroblasts", "Fibroblasts", "T cells", "Endothelial Cells", "Fibroblasts", "Myeloids", "Fibroblasts", "Fibroblasts", "Myeloids", "Epithelial cells", "Myeloids", "Myeloids")



new$seurat_clusters <- plyr::mapvalues(x = new$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)




###
# subset otherwise fibro and stroma

new <- subset(ct26,ident=c("Fibroblasts","Stromal cells"),invert=TRUE)

####

Cks2=Quasi_Mesanchymal_Adel
Anxa1: Tumor LiDing
Krt20:ClassicalB_Adel : Epithelial
Col3A1:newblast
Oasl2: eCAF_Raghu_Kalluri
C4b: new Raghu



png("Fibro_new.png", units="in", width=10, height=10, res=300)
DimPlot(new, reduction = "umap", label=F, group.by="copykat_pred")
FeaturePlot(new, reduction = "umap",features = c("Lgals1"), order = T, min.cutoff="q9", cols=c("lightgrey", "red"), label=F) # Tumor
FeaturePlot(new, reduction = "umap",features = c("Lgals1"), order = T, min.cutoff="q9", cols=c("lightgrey", "red"), label=F) # Tumor
FeaturePlot(new, reduction = "umap",features = c("Col1a1"), order = T, min.cutoff="q9", cols=c("lightgrey", "red"), label=F) # Fibro
FeaturePlot(new, reduction = "umap",features = c("Dcn"), order = T, min.cutoff="q9", cols=c("lightgrey", "red"), label=F) # Fibro
FeaturePlot(new, reduction = "umap",features = c("Dcn"), order = T, min.cutoff="q9", cols=c("lightgrey", "red"), label=F) # Fibro
FeaturePlot(new, reduction = "umap",features = c("Anxa2"), order = T, min.cutoff="q9", cols=c("lightgrey", "red"), label=F)
dev.off()




png("dimplot_immegen_ct263.png", units="in", width=20, height=10, res=300)
dittoDimPlot(ct26, "seurat_clusters", split.by = "singlr_labels_Immgen")
dev.off()

singlr_mouseRNA



png("dimplot_mouse_ct263.png", units="in", width=20, height=10, res=300)
dittoDimPlot(ct26, "seurat_clusters", split.by = "singlr_mouseRNA")
dev.off()

png("dimplot_immegen_ct263_2.png", units="in", width=15, height=15, res=300)
dittoDimPlot(ct26, "singlr_labels_Immgen",do.label=T)
dev.off()


png("dimplot_singlr_mouseRNA_ct263_2.png", units="in", width=15, height=15, res=300)
dittoDimPlot(ct26, "singlr_mouseRNA",do.label=T)
dev.off()



DimPlot(fibro, reduction = "umap", label=F, group.by="copykat_pred")



#####
#MIA


library(Seurat)
library(dplyr)
library(ggplot2)
library(reshape2)
library(scales)



SC <- readRDS("/N/project/cytassist/Integration/All_cells_merged.rds")
DefaultAssay(SC) <- "RNA"
SC <- DietSeurat(SC, assay = "RNA")
Idents(SC) <- SC$Celltype3
SC <- subset(SC, idents=c("Tumor_epithelium", "Normal_epithelium"), invert=TRUE)
SC <- ScaleData(SC)


ST <- readRDS("/N/project/cytassist/Integration/RCTD/Tier6/k10/ISCHIA_k10.rds")
#Idents(a) <- a$orig.ident
#samples <- c('S20-2878-C5', 'S20-3289-A3', 'S20-5578-A4', 'S20-6967-D7', 'S20-7543-A7', 'S21_10041_B1', 'S21_10041_B6', 'S21_11953_A1', 'S21_11953_A6', 'S21_13051_5', 'S21_33367_A4', 'S21_38552_7', 'S22_12298_A7', 'S22-21049-E8', 'S22-21554-B8', 'S22-28140-A4', 'S22-29157-A4', 'S22-32932-A7', 'S22-34391-C9', 'S22_3666_4', 'S22-40774-A5', 'S22_7205_A11', 'S22_8339_A8', 'S22_9021_A5')






#ST <- subset(a, idents=c(sample))
Idents(ST) <- ST$CompositionCluster_CC





# SC markers
Idents(SC) <- SC$functional.cluster
sc.markers <- FindAllMarkers(SC,assay = 'RNA',  only.pos = TRUE, test.use = 'wilcox', min.pct = 0.25, logfc.threshold = 0.5, verbose = F)
sc.markers['cluster'] %>% summary(maxsum=50) 






# ST markers
st.markers <- FindAllMarkers(ST, assay = 'RNA', only.pos = TRUE, test.use = 'wilcox', min.pct = 0.25, logfc.threshold = 0.5, verbose = F)
st.markers['cluster'] %>% summary()






# Create a list object containing the marker genes for each ST region:
st.clusts <- Idents(ST) %>% levels()
N <- length(st.clusts)


st.marker.list <- vector(mode = 'list', length = N)
names(st.marker.list) <- st.clusts
for(i in st.clusts) {
    st.marker.list[[i]] <- st.markers[st.markers$cluster == i,'gene']
}


# Create a list object containing the marker genes for each cell type:
sc.clusts <- Idents(SC) %>% levels()
M <- length(sc.clusts)


sc.marker.list <- vector(mode = 'list', length = M)
names(sc.marker.list) <- sc.clusts


for (i in sc.clusts) {
  sc.marker.list[[i]] <- sc.markers[sc.markers$cluster == i,'gene']
}


# Initialize a dataframe for us to store values in:
N <- length(st.clusts) ; M <- length(sc.clusts)
MIA.results <- matrix(0,nrow = M, ncol = N)
row.names(MIA.results) <- sc.clusts
colnames(MIA.results) <- st.clusts


# Gene universe
gene.universe <- intersect(ST@assays$Spatial %>% rownames(), SC@assays$RNA %>% rownames()) %>% length()


  # Loop over ST clusters
for (i in 1:N) {
  # Then loop over SC clusters
  for (j in 1:M) {
    genes1 <- st.marker.list[[st.clusts[i]]]
    genes2 <- sc.marker.list[[sc.clusts[j]]]
    
    # Hypergeometric    
    A <- length(intersect(genes1,genes2))
    B <- length(genes1)
    C <- length(genes2)
    enr <- -log10(phyper(A, B, gene.universe-B, C, lower.tail = FALSE))
    dep <- -log10(1-phyper(A, B, gene.universe-B, C, lower.tail = FALSE))
    if (enr < dep) {
      MIA.results[j,i] = -dep
    } else {
      MIA.results[j,i] = enr
    MIA.results[is.infinite(MIA.results)] <- 0
    }
  }
}
        


outfile <- "All_cells_MIA_results.rds"


saveRDS(MIA.results, outfile)


####
#Label_transfer


DefaultAssay(new) <- "RNA"
DefaultAssay(fibro) <- "RNA"

ifnb.list <- lapply(X = c(fibro,new), FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features1 <- SelectIntegrationFeatures(object.list = ifnb.list)
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features1)
# this command creates an 'integrated' data assay
immune.combined <- IntegrateData(anchorset = immune.anchors)
# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunTSNE(object = immune.combined, dims.use = 1:30)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined_0.5 <- FindClusters(immune.combined, resolution = 0.5)

#modify the metadata include CRC and BC
kefeer_meta1 <- read.csv("Meta_data_combined.csv",row.names=1)
immune.combined_0.5<- AddMetaData(immune.combined_0.5, kefeer_meta1)

DimPlot(immune.combined_0.5, reduction = "umap", split.by = "disease_type")


DimPlot(immune.combined, reduction = "umap", group.by = "cellline", label=T)



png("integrated_analysis4.png", units="in", width=20, height=20, res=400)
cowplot::plot_grid(ncol = 2,
DimPlot(fibro_sub, reduction = "umap", group.by = "cellline"),
DimPlot(fibro_sub, reduction = "umap", group.by = "celltypes_combined"),
DimPlot(fibro_sub, reduction = "umap", group.by = "orig.ident"),
DimPlot(fibro_sub, reduction = "umap", group.by = "copykat_pred"),
DimPlot(fibro_sub, reduction = "umap", label = FALSE, repel = TRUE))
dev.off()


png("Integration_analysis_fibroGenes.png", units="in", width=15, height=20, res=400)
FeaturePlot(fibro_sub, features = fib_genes, split.by = "cellline", max.cutoff = 3, cols = c("grey", "red"),label.size = 4)
dev.off()

# For performing differential expression after integration, we switch back to the original data
DefaultAssay(immune.combined) <- "RNA"
nk.markers <- FindConservedMarkers(immune.combined_0.5, ident.1 = 6, grouping.var = "disease_type", verbose = FALSE)
head(nk.markers)
dim(nk.markers)


####
#ref https://biocellgen-public.svi.edu.au/mig_2019_scrnaseq-workshop/comparing-and-combining-scrna-seq-datasets.html
library(scater)
library(SingleCellExperiment)
library(scmap)


emt6.sub.sce <- as.SingleCellExperiment(emt6)
rowData(emt6.sub.sce)$feature_symbol <- rownames(emt6.sub.sce)
#emt6.sub.sce <- emt6.sub.sce[!duplicated(rownames(emt6.sub.sce)), ]
emt6.sub.sce
emt6.sub.sce <- selectFeatures(emt6.sub.sce, suppress_plot = FALSE)

fibro.sce <- as.SingleCellExperiment(fibro)
rowData(fibro.sce)$feature_symbol <- rownames(fibro.sce)
#fibro.sce <- fibro.sce[!duplicated(rownames(fibro.sce)), ]
fibro.sce
fibro.sce <- selectFeatures(fibro.sce, suppress_plot = FALSE)

emt6.sub.sce <- indexCluster(emt6.sub.sce, cluster_col="celltypes")
fibro.sce <- indexCluster(fibro.sce, cluster_col="singlr_labels_Immgen") # copykat_pred
#fibro.sce <- indexCluster(fibro.sce, cluster_col="copykat_pred")


#We will project the ct26_fibro dataset to emt6 dataset:

emt6_to_ct26_fibro <- scmapCluster(
  projection = emt6.sub.sce,
  index_list = list(
    muraro = metadata(fibro.sce)$scmap_cluster_index
  )
)

plot(getSankey(colData(emt6.sub.sce)$celltypes,  emt6_to_ct26_fibro$scmap_cluster_labs[,1], plot_height=400))

#and muraro onto segerstolpe

fibro_to_emt6 <- scmapCluster(
  projection = fibro.sce,
  index_list = list(
    seger = metadata(emt6.sub.sce)$scmap_cluster_index
  )
)

plot(getSankey(colData(fibro.sce)$celltypes,  fibro_to_emt6$scmap_cluster_labs[,1], plot_height=400))


genes <- c("Hspb1","Lgals1","Cald1","Rbp1","Cd63","Ptn")


multi_dittoPlot(new, group.by = "new_celltype", vlnplot.lineweight = 0.2, jitter.size = 0.3)

#####



  

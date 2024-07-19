### Subsetting CAFs ###
#45426 sc Cells

library(Seurat)
library(dbplyr)

#data <- Read10X(data.dir = "C:/Users/haan1/Downloads/vdj_v1_hs_nsclc_5gex_filtered_gene_bc_matrices/")
#caf <- CreateSeuratObject(counts = data, project = "Ha An", min.cells = 3, min.features = 200)
#caf[["percent.mt"]] <- PercentageFeatureSet(caf, pattern = "^MT-")
#caf <- subset(caf, subset = nFeature_RNA>200 & nFeature_RNA<2500 & percent.mt <5)

caf_cont2 <- NormalizeData(caf_cont2, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
caf_cont2 <- FindVariableFeatures(caf_cont2, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
caf_cont2 <- ScaleData(caf_cont2, verbose = FALSE)
caf_cont2 <- RunPCA(caf_cont2, verbose = FALSE, features = VariableFeatures(object = caf_cont2))
caf_cont2 <- RunTSNE(object = caf_cont2, dims.use = 1:30)
caf_cont2 <- RunUMAP(caf_cont2, dims = 1:30)

caf_cont2 <- FindNeighbors(caf_cont2, dims = 1:30, force.recalc = T)
caf_cont2 <- FindClusters(caf_cont2, force.recalc = T,resolution = 0.8 , algorithm = 1)

caf_cont2$seurat_clusters <- as.factor(as.numeric(as.character(caf_cont2$seurat_clusters)) + 1)
Idents(caf_cont2) <- "seurat_clusters"

 DimPlot(caf_cont2, reduction = "umap", label = T)

markers_genes_caf_cont2_0.8 <- FindAllMarkers(caf_cont2, logfc.threshold = 0.2, test.use = "wilcox", only.pos = TRUE,assay = "RNA")
write.table(markers_genes_caf_cont2_0.8, file="marker_genes_caf_cont_0.8.txt", sep="\t",append=F, row.names = F)
cbind(table(markers_genes_caf_cont2_0.8$cluster))

## removing Contamination

Idents(caf_cont)<-"seurat_clusters"
#caf_cont <- subset(caf,ident=c(2,3,4,5,6,7,8,9,10,11,13,15,16,17,18,19,20,21,22,23,25,26,27,28,30,31,32,33,34,36,37,38,39,40,41,42,43,44,45,46))

caf_cont <- subset(caf,ident=c(1,2,3,4,5,6,7,8,9,10,12,14,15,18,19,20,21,23,24,27,28,29,30,33,34,36,38,39,40,41,42,45,46,47,49,50,51,53,55,56,57,58,59,60,61,62,63,64,65,66,67,68))
##
caf_cont1 <- subset(caf_cont,ident=c(1,2,3,4,6,7,9,10,11,12,13,14,15,16,17,18,20,21,23,24,26,27,28,29,30,33,34,35,36,37,38,39,40,41,42))

caf_cont2 <- subset(caf_cont1,ident=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,17,18,19,20,21,23,24,25))



genes = c("Pdgfra","Gpc3","Antxr1","Sdc1","Cxcl12","Dlk1","Scara5","Col1a1","Col1a2","Fap","Pdpn","Lamp5","Tgfb1","Cd9","Sema3c","Cd248","Cspg4","Epas1","Pdgfrb","Rgs5","Cspg4","Mcam")

png("heatmap_cafs.png", units="in", width=10, height=10, res=300)
dittoHeatmap(caf_cont2, genes, annot.by = c("seurat_clusters"), heatmap.colors = colorRampPalette(c("blue", "red"))(50))
dev.off()


png("heatmap_fibro_markers.png", units="in",width=10, height=15,res=300)
DoHeatmap(object = caf_cont2,genes) + scale_fill_viridis(option="magma")
 dev.off()

DotPlot(object = caf_cont2, features = genes)+ scale_colour_gradient2(low="steelblue", mid="lightgrey", high="darkgoldenrod1")

DotPlot(object = caf_cont2, features = genes, assay="RNA")+guides(color = guide_colorbar(title = 'Scaled Average Expression')) + theme(axis.text.x = element_text(angle=90))+ scale_colour_gradient2(low = "#FF00FF", mid = "#000000", high = "#FFFF00")


dittoDotPlot(caf_cont2, vars = genes, group.by = "seurat_clusters",min.color = "#2eafff",max.col="#ff4d00")

png("meyloid_markers_mycaf.png", units="in",width=8, height=10,res=300)
multi_dittoDimPlot(obj, mycaf.genes,reduction.use ="tsne",axes.labels.show=FALSE,legend.show = TRUE,min.color="lightgrey", max.color="red")
dev.off()
 
 mycafs <- c("Acta2","Tagln","Myl9","Tpm1","Tpm2","Mmp11","Postn","Hopx","Fap","Pdpn","Col12a1","Pdgfrb")
png("mycaf_markers.png", units="in",width=8, height=10,res=300) 
multi_dittoDimPlot(caf_cont2, mycafs,reduction.use ="umap",axes.labels.show=FALSE,legend.show = TRUE,min.color="lightgrey", max.color="red")
dev.off()

 png("mycaf_markers_dotPlot.png", units="in",width=10, height=10,res=300) 
DotPlot(object = caf_cont2, features = mycafs, assay="RNA")+guides(color = guide_colorbar(title = 'Scaled Average Expression')) + theme(axis.text.x = element_text(angle=90))+ scale_colour_gradient2(low = "#FF00FF", mid = "#000000", high = "#FFFF00")
dev.off()

png("mycaf_markers_ViolinPlot.png", units="in",width=10, height=10,res=300) 
multi_dittoPlot(caf_cont2, mycafs, group.by = "seurat_clusters", vlnplot.lineweight = 0.2, jitter.size = 0)
dev.off()

icafs <- c("Il6","Cxcl2","Cxcl12","Ccl2","Cfd","Lmna","Dpt","Clu","Emp1","Has1","Pdgfra") # not present = Il8,Cxcl1

png("icaf_markers.png", units="in",width=8, height=10,res=300) 
multi_dittoDimPlot(caf_cont2, icafs,reduction.use ="umap",axes.labels.show=FALSE,legend.show = TRUE,min.color="lightgrey", max.color="red")
dev.off()

 png("icaf_markers_dotPlot.png", units="in",width=10, height=10,res=300) 
DotPlot(object = caf_cont2, features = icafs, assay="RNA")+guides(color = guide_colorbar(title = 'Scaled Average Expression')) + theme(axis.text.x = element_text(angle=90))+ scale_colour_gradient2(low = "#FF00FF", mid = "#000000", high = "#FFFF00")
dev.off()

png("caf_markers_ViolinPlot.png", units="in",width=10, height=8,res=300) 
multi_dittoPlot(caf_cont2, icafs, group.by = "seurat_clusters", vlnplot.lineweight = 0.2, jitter.size = 0)
dev.off()

apcafs <- c("Cd74","Saa3")# not present H2-aa,H2-ab1,Slp1


png("apcaf_markers.png", units="in",width=8, height=10,res=300) 
multi_dittoDimPlot(caf_cont2, apcafs,reduction.use ="umap",axes.labels.show=FALSE,legend.show = TRUE,min.color="lightgrey", max.color="red")
dev.off()

 png("apcaf_markers_dotPlot.png", units="in",width=10, height=10,res=300) 
DotPlot(object = caf_cont2, features = apcafs, assay="RNA")+guides(color = guide_colorbar(title = 'Scaled Average Expression')) + theme(axis.text.x = element_text(angle=90))+ scale_colour_gradient2(low = "#FF00FF", mid = "#000000", high = "#FFFF00")
dev.off()

png("apcaf_markers_ViolinPlot.png", units="in",width=10, height=3,res=300) 
multi_dittoPlot(caf_cont2, apcafs, group.by = "seurat_clusters", vlnplot.lineweight = 0.2, jitter.size = 0)
dev.off()


fibro <- c("Col1a1","Col1a2","")

####
library(RColorBrewer)
library(colorspace)

# Define multiple palettes from RColorBrewer without exceeding their limits
palette1 <- brewer.pal(12, "Set3")
palette2 <- brewer.pal(12, "Paired")
palette3 <- brewer.pal(8, "Dark2")
palette4 <- brewer.pal(9, "Set1")
palette5 <- brewer.pal(8, "Set2")
palette6 <- brewer.pal(8, "Accent")
palette7 <- brewer.pal(7, "Pastel1")

# Combine these palettes
combined_palettes <- c(palette1, palette2, palette3, palette4, palette5, palette6, palette7)

# Check the length of combined palettes
length(combined_palettes)  # Should be less than 69

# Generate additional colors using colorspace
if (length(combined_palettes) < 69) {
  additional_colors <- rainbow_hcl(69 - length(combined_palettes))
  pal <- c(combined_palettes, additional_colors)
} else {
  pal <- combined_palettes[1:69]  # In case the combined palettes already have enough colors
}

# Check the palette
print(pal)
length(pal)  # Should be exactly 69

# Use the palette in your plot function


pdf("scdata_plot.pdf", width = 10, height = 8)
plot_scdata(caf, pal_setup = pal)
dev.off()

####



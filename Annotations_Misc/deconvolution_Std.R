  

### to check top genes expressed in each topic
topGenes(deconGexp, n = 10)



## Filtering out cell-types in pixels that contribute less than 0.05 of the pixel proportion.
deconProp <- results$theta
deconGexp <- results$beta
## visualize deconvolved cell-type proportions


pos <- GetTissueCoordinates(pdac.integrated, image="PDACH_1")
colnames(pos) <- c('x','y')
annot <- pdac.integrated$seurat_clusters

png("Check.png", units="in", width=10, height=10, res=300)
vizAllTopics(deconProp, pos,
             groups = annot, 
             group_cols = rainbow(length(levels(annot))),
             r=3)
dev.off()

rownames(deconProp)<-gsub(pattern = ".", replacement = "-",x = rownames(deconProp), fixed = TRUE)

# We can also use vizTopic() for faster plotting. Let’s visualize each deconvolved cell-type separately:
ps <- lapply(colnames(deconProp), function(celltype) {
  
  vizTopic(theta = deconProp, pos = pos, topic = celltype, plotTitle = paste0("X", celltype),
         size = 2, stroke = 1, alpha = 0.5,
         low = "white",
         high = "red") +
    
    ## remove the pixel "Groups", which is the color aesthetic for the pixel borders
    ggplot2::guides(colour = "none")
  
})
gridExtra::grid.arrange(
  grobs = ps,
  layout_matrix = rbind(c(1, 2, 3, 4),
                        c(5, 6, 7, 8),
                        c(9, 10, 11, 12),
                        c(13, 14, 15, 16))
)





######  ANNOTATIONS #########

# proxy theta for the annotated layers
mobProxyTheta <- model.matrix(~ 0 + annot)
rownames(mobProxyTheta) <- names(annot)
# fix names
colnames(mobProxyTheta) <- unlist(lapply(colnames(mobProxyTheta), function(x) {
  unlist(strsplit(x, "annot"))[2]
}))

mobProxyGexp <- counts %*% mobProxyTheta

# Strategy 1 Corealation 

## row and column names need to be characters
rownames(corMtx_beta) <- paste0("decon_", seq(nrow(corMtx_beta)))


png("st.pearson_cr_trans_prof_beta.png", units="in", width=10, height=10, res=300)
correlationPlot(mat = corMtx_beta,
                colLabs = "Deconvolved cell-types", # aka x-axis, and rows of matrix
                rowLabs = "Ground truth cell-types", # aka y-axis, and columns of matrix
                title = "Transcriptional correlation", annotation = TRUE) +
  
  ## this function returns a `ggplot2` object, so can add additional aesthetics
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0))

dev.off()


corMtx_theta <- getCorrMtx(m1 = as.matrix(deconProp), # the deconvolved cell-type `theta` (pixels x celltypes)
                          m2 = as.matrix(mobProxyTheta), # the reference `theta` (pixels x celltypes)
                          type = "t") # "b" = comparing beta matrices, "t" for thetas
## cell-type correlations based on 260 shared pixels between m1 and m2.
## row and column names need to be characters
rownames(corMtx_theta) <- paste0("decon_", seq(nrow(corMtx_theta)))

png("st.pearson_cr_trans_prof_Theta.png", units="in", width=10, height=10, res=300)
correlationPlot(mat = corMtx_theta,
                colLabs = "Deconvolved cell-types", # aka x-axis, and rows of matrix
                rowLabs = "Ground truth cell-types", # aka y-axis, and columns of matrix
                title = "Proportional correlation", annotation = TRUE) +
  
  ## this function returns a `ggplot2` object, so can add additional aesthetics
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0))

dev.off()

## order the cell-types rows based on best match (highest correlation) with each community
## cannot have more rows than columns for this pairing, so transpose
pairs <- lsatPairs(t(corMtx_theta))
m <- t(corMtx_theta)[pairs$rowix, pairs$colsix]

png("st.pearson_cr_trans_prof_highest.png", units="in", width=10, height=10, res=300)
correlationPlot(mat = t(m), # transpose back
                colLabs = "Deconvolved cell-types", # aka x-axis, and rows of matrix
                rowLabs = "Ground truth cell-types", # aka y-axis, and columns of matrix
                title = "Transcriptional correlation", annotation = TRUE) +
  
  ## this function returns a `ggplot2` object, so can add additional aesthetics
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0))

dev.off()


# Strategy 2 GSEA analysis

mobProxyLayerMarkers <- list()

## make the tissue layers the rows and genes the columns
gexp <- t(as.matrix(mobProxyGexp))

for (i in seq(length(rownames(gexp)))){
  celltype <- i
  ## log2FC relative to other cell-types
  ## highly expressed in cell-type of interest
  highgexp <- names(which(gexp[celltype,] > 10))
  ## high log2(fold-change) compared to other deconvolved cell-types and limit to top 200
  log2fc <- sort(log2(gexp[celltype,highgexp]/colMeans(gexp[-celltype,highgexp])), decreasing=TRUE)[1:200]
  
  ## for gene set of the ground truth cell-type, get the genes
  ## with log2FC > 1 (so FC > 2 over the mean exp of the other cell-types)
  markers <- names(log2fc[log2fc > 1])
  mobProxyLayerMarkers[[ rownames(gexp)[celltype] ]] <- markers
}

celltype_annotations <- annotateCellTypesGSEA(beta = results$beta, gset = mobProxyLayerMarkers, qval = 0.05)



#


pdac.integrated <- RunPCA(pdac.integrated, verbose = FALSE)
pdac.integrated <- FindNeighbors(pdac.integrated, dims = 1:30)
pdac.integrated <- FindClusters(pdac.integrated, verbose = FALSE, resolution = 0.8, algorithm = 1) # @ 0.8 res
pdac.integrated <- RunUMAP(pdac.integrated, dims = 1:30)
pdac.integrated <- RunTSNE(pdac.integrated, dims = 1:30)


pdac.integrated$seurat_clusters <- as.factor(as.numeric(as.character(pdac.integrated$seurat_clusters)) + 1)
Idents(pdac.integrated) <- "seurat_clusters"



markers_genes_st <- FindAllMarkers(pdac.integrated, logfc.threshold = 0.2, test.use = "wilcox", only.pos = TRUE)
write.table(markers_genes_st, file="marker_genes_spatial_0.5.txt", sep="\t",append=F, row.names = F)
cbind(table(markers_genes_st$cluster))

source("save_markers_to_excel_pages_by_cluster.R")
save_markers_to_excel_pages_by_cluster("marker_genes_spatial_0.8.txt","Markers_0.8.xlsx")




old.cluster.ids <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)
new.cluster.ids <- c("Stem cells", "Epithelial cells","Myeloid","Acinar","Fibroblasts","B cells","Myeloid","Fibroblasts","T cells","Epithelial cells","B cells","B cells","Endothelial cells", "Epithealial cells", "T cells","Myeloid","Stem cells","Acinar","Unknown")
pdac.integrated$seurat_clusters <- plyr::mapvalues(x = pdac.integrated$seurat_clusters, from = old.cluster.ids, to = new.cluster.ids)
kable(table(pdac.integrated$seurat_clusters))

old.cluster.ids <- c("Stem cells", "Epithelial cells","Myeloid","Acinar","Fibroblasts","B cells","Myeloid","Fibroblasts","T cells","Epithelial cells","B cells","B cells","Endothelial cells", "Epithealial cells", "T cells","Myeloid","Stem cells","Acinar","Unknown")
new.cluster.ids <- c("Stem cells", "Epithelial cells","Myeloid","Acinar","Fibroblasts","B cells","Myeloid","Fibroblasts","T cells","Epithelial cells","B cells","B cells","Endothelial cells", "Epithelial cells", "T cells","Myeloid","Stem cells","Acinar","Unknown")
pdac.integrated$seurat_clusters <- plyr::mapvalues(x = pdac.integrated$seurat_clusters, from = old.cluster.ids, to = new.cluster.ids)
kable(table(pdac.integrated$seurat_clusters))


png("clusters_Spatial.png", units="in", width=10, height=10, res=300)
dittoDimPlot(pdac.integrated, reduction = "umap", "seurat_clusters",do.label= TRUE)
dev.off()






######## step 3: output plot of how the resolution changes the number of clusters you get
n.clusters = vector(mode="numeric", length=length(set.res))
names(n.clusters) = set.res
  for(i in 1:length(n.clusters)){
      n.clusters[i] = length(table(as.vector(srobj.tmp@meta.data[,paste("integrated_snn_res.", names(n.clusters)[i], sep="")])))
  }


### Run MCP

library(dplyr)
library(ggplot2)
library(tidyr)
library(immunedeconv)
library(tibble)


#get normalised Counts

pdac.nor <- pdac[["integrated"]]@scale.data

# for sample 

pd.samp <- AggregateExpression(pdac,group.by="orig.ident",return.seurat = TRUE)
res_quantiseq.samp <- deconvolute(pd.samp[["integrated"]]@scale.data, "quantiseq", tumor = TRUE)
res_mcp_counter.samp <- deconvolute(pd.samp[["integrated"]]@scale.data, "mcp_counter")

png("pdac.sample.quatiseq.png", units="in", width=10, height=10, res=300)
res_quantiseq.samp %>%
  gather(sample, fraction, -cell_type) %>%
  # plot as stacked bar chart
  ggplot(aes(x = sample, y = fraction, fill = cell_type)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_brewer(palette = "Paired") +
  scale_x_discrete(limits = rev(levels(res_quantiseq)))
dev.off()

png("pdac.sample.mcp.png", units="in", width=10, height=15, res=300)
res_mcp_counter.samp %>%
  gather(sample, score, -cell_type) %>%
  ggplot(aes(x = sample, y = score, color = cell_type)) +
  geom_point(size = 4) +
  facet_wrap(~cell_type, scales = "free_x", ncol = 3) +
  scale_color_brewer(palette = "Paired", guide = FALSE) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()


mcp.est.samp=MCPcounter.estimate(pd.samp[["integrated"]]@scale.data,featuresType="HUGO_symbols")

png("pdac.sample.mcp2.png", units="in", width=10, height=10, res=300)
heatmap(as.matrix(mcp.est.samp),col=colorRampPalette(c("blue","white","red"))(100)) 
dev.off()



# for tumor_type2

pd.tt2 <- AggregateExpression(pdac,group.by="tumor_type2",return.seurat = TRUE)
res_quantiseq.tt2 <- deconvolute(pd.tt2[["integrated"]]@scale.data, "quantiseq", tumor = TRUE)
res_mcp_counter.tt2 <- deconvolute(pd.tt2[["integrated"]]@scale.data, "mcp_counter")

png("pdac.tt2.quatiseq.png", units="in", width=10, height=10, res=300)
res_quantiseq.tt2 %>%
  gather(sample, fraction, -cell_type) %>%
  # plot as stacked bar chart
  ggplot(aes(x = sample, y = fraction, fill = cell_type)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_brewer(palette = "Paired") +
  scale_x_discrete(limits = rev(levels(res_quantiseq)))
dev.off()

png("pdac.tt2.mcp.png", units="in", width=10, height=15, res=300)
res_mcp_counter.tt2 %>%
  gather(sample, score, -cell_type) %>%
  ggplot(aes(x = sample, y = score, color = cell_type)) +
  geom_point(size = 4) +
  facet_wrap(~cell_type, scales = "free_x", ncol = 3) +
  scale_color_brewer(palette = "Paired", guide = FALSE) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()


mcp.est.tt2=MCPcounter.estimate(pd.tt2[["integrated"]]@scale.data,featuresType="HUGO_symbols")

png("pdac.tt2.mcp2.png", units="in", width=10, height=10, res=300)
heatmap(as.matrix(mcp.est.tt2),col=colorRampPalette(c("blue","white","red"))(100)) 
dev.off()



# for tumor_type

pd.tt <- AggregateExpression(pdac,group.by="tumor_type",return.seurat = TRUE)
res_quantiseq.tt <- deconvolute(pd.tt[["integrated"]]@scale.data, "quantiseq", tumor = TRUE)
res_mcp_counter.tt <- deconvolute(pd.tt[["integrated"]]@scale.data, "mcp_counter")

png("pdac.tt.quatiseq.png", units="in", width=10, height=10, res=300)
res_quantiseq.tt %>%
  gather(sample, fraction, -cell_type) %>%
  # plot as stacked bar chart
  ggplot(aes(x = sample, y = fraction, fill = cell_type)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_brewer(palette = "Paired") +
  scale_x_discrete(limits = rev(levels(res_quantiseq)))
dev.off()

png("pdac.tt.mcp.png", units="in", width=10, height=15, res=300)
res_mcp_counter.tt %>%
  gather(sample, score, -cell_type) %>%
  ggplot(aes(x = sample, y = score, color = cell_type)) +
  geom_point(size = 4) +
  facet_wrap(~cell_type, scales = "free_x", ncol = 3) +
  scale_color_brewer(palette = "Paired", guide = FALSE) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()


mcp.est.tt=MCPcounter.estimate(pd.tt[["integrated"]]@scale.data,featuresType="HUGO_symbols")

png("pdac.tt.mcp2.png", units="in", width=10, height=10, res=300)
heatmap(as.matrix(mcp.est.tt),col=colorRampPalette(c("blue","white","red"))(100)) 
dev.off()







pd.tmr.typ2 <- AggregateExpression(pdac,group.by="tumor_type2",return.seurat = TRUE)




res_quantiseq <- deconvolute(pd.tmr.typ2[["integrated"]]@scale.data, "quantiseq", tumor = TRUE)
res_mcp_counter <- deconvolute(pd.tmr.typ2[["integrated"]]@scale.data, "mcp_counter")

Examples:

     mcp.est=MCPcounter.estimate(pdac.nor,featuresType="HUGO_symbols")
     heatmap(as.matrix(ExampleEstimates),col=colorRampPalette(c("blue","white","red"))(100)) 
     

mcp.est %>%
  gather(sample, score, -cell_type) %>%
  ggplot(aes(x = sample, y = score, color = cell_type)) +
  geom_point(size = 4) +
  facet_wrap(~cell_type, scales = "free_x", ncol = 3) +
  scale_color_brewer(palette = "Paired", guide = FALSE) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    
##make first row as colname
library(tibble)
mcp.est <- tibble::rownames_to_column(mcp.est, "cell_type")


png("example.png", units="in", width=10, height=10, res=300)
mcp.est %>%
  gather(sample, fraction, -cell_type) %>%
  # plot as stacked bar chart
  ggplot(aes(x = sample, y = fraction, fill = cell_type)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_brewer(palette = "Paired") +
  scale_x_discrete(limits = rev(levels(res_quantiseq)))
dev.off()


#######
rownames(deconProp)<-gsub(pattern = ".", replacement = "-",x = rownames(deconProp), fixed = TRUE)

Compare to transcriptional clustering

pcs.info <- stats::prcomp(t((pdac@assays$integrated@scale.data)+1), center=TRUE)
nPcs <- 30 ## let's take the top 5 PCs
pcs <- pcs.info$x[,1:nPcs]


emb <- Rtsne::Rtsne(pcs,
             is_distance=FALSE,
             perplexity=30,
             num_threads=8,
             verbose=FALSE)$Y
rownames(emb) <- rownames(pcs)
colnames(emb) <- c("x", "y")

png("std_cluster.png", units="in", width=10, height=10, res=300)
plt 
dev.off()


ps <- lapply(colnames(deconProp), function(celltype) {
  
  vizTopic(theta = deconProp, pos = emb, topic = celltype, plotTitle = paste0("X", celltype),
         size = 1, stroke = 0.5, alpha = 0.5,
         low = "white",
         high = "red") +
    
    ## remove the pixel "Groups", which is the color aesthetic for the pixel borders
    ggplot2::guides(colour = "none")
  
})
png("prop_on_cluster.png", units="in", width=10, height=10, res=300)
gridExtra::grid.arrange(
  grobs = ps,
  layout_matrix = rbind(c(1, 2, 3, 4),
                        c(5, 6, 7, 8),
                        c(9, 10, 11, 12),
                        c(13, 14, 15, 16),
                        c(17, 18,19,20))
)
dev.off()



png("cor_plot_cluster.png", units="in", width=10, height=10, res=300)
STdeconvolve::correlationPlot(mat = m,
                              colLabs = "Transcriptional clusters",
                              rowLabs = "STdeconvolve") +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 90)
)
dev.off()



# Annotation agai

corMtx_beta <- getCorrMtx(m1 = as.matrix(deconGexp), # the deconvolved cell-type `beta` (celltypes x genes)
                          m2 = t(as.matrix(com_proxyTheta)), # the reference `beta` (celltypes x genes)
                          type = "t") # "b" = comparing beta matrices, "t" for thetas


png("cor_plot_2.png", units="in", width=10, height=10, res=300)
correlationPlot(mat = corMtx_beta,
                colLabs = "Deconvolved cell-types", # aka x-axis, and rows of matrix
                rowLabs = "Ground truth cell-types", # aka y-axis, and columns of matrix
                title = "Proportional correlation", annotation = TRUE) +
  
  ## this function returns a `ggplot2` object, so can add additional aesthetics
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0))
dev.off()

#####

#MCPcounter
#CRC Liver Mets 
Idents(pdac) <- "tumor_type"

kable(table(pdac$tumor_type))
crc.l.m <- subset(pdac, idents="CRC Liver Mets")
crc.l.m.ag <- AggregateExpression(crc.l.m,group.by="orig.ident",return.seurat = TRUE)

crc.l.m.est=MCPcounter.estimate(as.matrix(crc.l.m.ag@assays$integrated@counts),featuresType="HUGO_symbols")

png("crc_liver_mets.mcp.png", units="in", width=10, height=10, res=300)
heatmap(as.matrix(crc.l.m.est),col=colorRampPalette(c("blue","white","red"))(100)) 
dev.off()

png("crc_liver_mets.mcp.png", units="in", width=10, height=10, res=300)
pheatmap(crc.l.m.est,cutree_rows = 4,cutree_cols = 3)
dev.off()

#Lymph Node
ly.n <- subset(pdac, idents="Lymph Node")
ly.n.ag <- AggregateExpression(ly.n,group.by="orig.ident",return.seurat = TRUE)

ly.n.est=MCPcounter.estimate(as.matrix(ly.n.ag@assays$integrated@counts),featuresType="HUGO_symbols")

png("lymph_nodes.mcp.png", units="in", width=20, height=20, res=300)
heatmap(as.matrix(ly.n.est),col=colorRampPalette(c("blue","white","red"))(100)) 
dev.off()

png("lymph_nodes.mcp.pheatmap.png", units="in", width=10, height=10, res=300)
pheatmap(ly.n.est,cutree_rows = 4,cutree_cols = 3)
dev.off()




#tumorType2

#tu.2 <- subset(pdac, idents="tumor_type2")
tu.2.ag <- AggregateExpression(pdac,group.by="tumor_type2",return.seurat = TRUE)

tu.2.est=MCPcounter.estimate(as.matrix(tu.2.ag@assays$integrated@counts),featuresType="HUGO_symbols")

png("tumor_type2.mcp.png", units="in", width=10, height=10, res=300)
heatmap(as.matrix(tu.2.est),col=colorRampPalette(c("blue","white","red"))(100))
dev.off()

png("tumorType2.mcp.pheatmap.png", units="in", width=10, height=10, res=300)
pheatmap(tu.2.est,cutree_rows = 4,cutree_cols = 3)
dev.off()


orig.ident <- subset(pdac, idents="orig.ident")
orig.ident.ag <- AggregateExpression(pdac,group.by="orig.ident",return.seurat = TRUE)

orig.ident.est=MCPcounter.estimate(as.matrix(orig.ident.ag@assays$integrated@counts),featuresType="HUGO_symbols")


png("samp.mcp.pheatmap.png", units="in", width=10, height=10, res=300)
pheatmap(as.matrix(heat),scale="row", annotation_col = anno,color=colorRampPalette(c("navy", "white", "red"))(50))
dev.off()


col = list(Samples = colorRampPalette(brewer.pal(12,"Paired"))
ha <- HeatmapAnnotation( Samples = ann$Samples, TissueType = ann$TissueType, col = list(Samples = c("Liver Mets" = "yellow", "lymph node" = "green", "Normal Liver"="#1e81b0
","Normap Pancreas"="#873e23","PDAC"="#760D0D")) )

png("samp.mcp.pheatmap1.png", units="in", width=20, height=10, res=300)
Heatmap(t(heat), name = "MCP Estimates", top_annotation = ha, col=hmcol)
dev.off()


######
# PDAC all
library("MCPcounter")
pdac.mcp <- subset(pdac, idents=c("Normal Pancreas","PDAC","PDAC Liver Mets"))
pdac.m.ag <- AggregateExpression(pdac.mcp,group.by="orig.ident",return.seurat = TRUE)
pdac.m.est=MCPcounter.estimate(as.matrix(pdac.m.ag@assays$integrated@counts),featuresType="HUGO_symbols")
write.csv(pdac.m.est,"pdac.m.est.csv")

#Set annotation
samp.anno<- read.csv("samp_anno.csv",header=TRUE,sep=',')
colours <- list('Tissue' = c('Normal Pancreas' = 'limegreen', 'PDAC' = 'red', 'PDAC Liver Mets' = 'blue'))
colAnn <- HeatmapAnnotation(df = samp.anno,which = 'col', col = colours, annotation_width = unit(c(1, 4), 'cm'),gap = unit(1, 'mm'))

hmap <- Heatmap(pdac.m.est, name = "MCP Counter results", show_row_names = TRUE,show_column_names = TRUE,cluster_rows = TRUE,cluster_columns = FALSE,show_column_dend = TRUE,show_row_dend = TRUE,row_dend_reorder = TRUE,column_dend_reorder = TRUE,clustering_method_rows = "ward.D2",clustering_method_columns = "ward.D2",width = unit(100, "mm"),top_annotation=colAnn)

png("pdac.all.mcp.png", units="in", width=10, height=10, res=300)
draw(hmap, heatmap_legend_side="left", annotation_legend_side="right")
dev.off()


# PDAC only
library("MCPcounter")
pdac.only.mcp <- subset(pdac, idents=c("Normal Pancreas","PDAC"))
pdac.only.m.ag <- AggregateExpression(pdac.only.mcp,group.by="orig.ident",return.seurat = TRUE)
pdac.only.m.est=MCPcounter.estimate(as.matrix(pdac.only.m.ag@assays$integrated@counts),featuresType="HUGO_symbols")
write.csv(pdac.only.m.est,"pdac.only.m.est.csv")

#Set annotation
samp.anno<- read.csv("samp_anno.csv",header=TRUE,sep=',')
colours <- list('Tissue' = c('Normal Pancreas' = 'limegreen', 'PDAC' = 'red', 'PDAC Liver Mets' = 'blue'))
colAnn <- HeatmapAnnotation(df = samp.anno,which = 'col', col = colours, annotation_width = unit(c(1, 4), 'cm'),gap = unit(1, 'mm'))

hmap <- Heatmap(pdac.only.m.est, name = "MCP Counter results", show_row_names = TRUE,show_column_names = TRUE,cluster_rows = TRUE,cluster_columns = FALSE,show_column_dend = TRUE,show_row_dend = TRUE,row_dend_reorder = TRUE,column_dend_reorder = TRUE,clustering_method_rows = "ward.D2",clustering_method_columns = "ward.D2",width = unit(100, "mm"),top_annotation=colAnn)

png("pdac.only.mcp.png", units="in", width=10, height=10, res=300)
draw(hmap, heatmap_legend_side="left", annotation_legend_side="right")
dev.off()

# PDAC all samples
library("MCPcounter")
pdac.all.mcp <- subset(pdac, idents=c("Normal Pancreas","PDAC"))
pdac.all.m.ag <- AggregateExpression(pdac,group.by="orig.ident",return.seurat = TRUE)
pdac.all.m.est=MCPcounter.estimate(as.matrix(pdac.all.m.ag@assays$integrated@counts),featuresType="HUGO_symbols")
write.csv(pdac.all.m.est,"pdac.only.all.est.csv")

#Set annotation
samp.anno<- read.csv("samp_anno.csv",header=TRUE,sep=',')
colours <- list('Tissue' = c('Normal Pancreas' = 'limegreen', 'PDAC' = 'red', 'Liver Mets' = 'blue','Lymph Node'='yellow','CRC Liver'='#ff66ff','Normal Liver' = 'cyan',))
colAnn <- HeatmapAnnotation(df = samp.anno,which = 'col', annotation_width = unit(c(1, 4), 'cm'),gap = unit(1, 'mm'))

hmap <- Heatmap(pdac.all.m.est, name = "MCP Counter results", show_row_names = TRUE,show_column_names = TRUE,cluster_rows = TRUE,cluster_columns = FALSE,show_column_dend = TRUE,show_row_dend = TRUE,row_dend_reorder = TRUE,column_dend_reorder = TRUE,clustering_method_rows = "ward.D2",clustering_method_columns = "ward.D2",width = unit(100, "mm"),top_annotation=colAnn)

png("pdac.allSamples.mcp.png", units="in", width=15, height=10, res=300)
draw(hmap, heatmap_legend_side="left", annotation_legend_side="right")
dev.off()


|Var1            | Freq|
|:---------------|----:|
|CRC Liver       |    4|
|Liver Mets      |   11|
|Lymph node      |    5|
|Normal Liver    |    2|
|Normal Liver    |    1|
|Normal Pancreas |    3|
|PDAC            |   10|
> 





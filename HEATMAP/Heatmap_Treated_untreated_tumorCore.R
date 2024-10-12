library(Seurat)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(grDevices)

# List of genes to plot
genes_to_plot <- c("ALDH1A3", "AXL", "CD24", "CD44", "CLDN2", "DDR1", "EPCAM", "EPHB2", "GAS6", "ITGA3", "LGR5", "MET", "OLFM4", "PROM1", "RNF43", "S100A4", "SOX9", "TGFBI", "ALDH1A1", "AXIN2", "CXCR4", "SERPINA3")

# Subset the Seurat object
tumor_core_subset <- subset(tumor_core, features = genes_to_plot)

# Extract normalized expression data
expression_data <- GetAssayData(tumor_core_subset, slot = "data", assay = "SCT")

# Calculate average expression for each gene in each NAC treatment group
nac_groups <- c("FOLFIRINOX", "Gem_Abx", "Gemcitabine", "Mix", "Naive")
avg_expression <- sapply(nac_groups, function(group) {
  rowMeans(expression_data[, tumor_core_subset$nac_treatment == group, drop = FALSE])
})

# Calculate Z-scores
z_scores <- t(scale(t(avg_expression)))

# Define colors
heatmap_colors <- colorRampPalette(c("#0000FF", "white", "#FF0000"))(100)

# Create annotation for NAC treatment groups and treated status
annotation_col <- data.frame(
  NAC_Treatment = factor(nac_groups, levels = nac_groups),
  Treated = ifelse(nac_groups == "Naive", "No", "Yes")
)
rownames(annotation_col) <- colnames(z_scores)

# Define colors for annotation
ann_colors <- list(
  NAC_Treatment = c("FOLFIRINOX" = "#E41A1C", "Gem_Abx" = "#377EB8", 
                    "Gemcitabine" = "#4DAF4A", "Mix" = "#984EA3", "Naive" = "#FF7F00"),
  Treated = c("No" = "#377EB8", "Yes" = "#E41A1C")
)

# Open PDF device
pdf("gene_expression_heatmap_nac_treated.pdf", width = 12, height = 10)

# Create the heatmap
pheatmap(z_scores,
         color = heatmap_colors,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_colnames = TRUE,
         show_rownames = TRUE,
         fontsize_row = 10,
         fontsize_col = 12,
         angle_col = 45,
         border_color = NA,
         cellwidth = 40,
         cellheight = 12,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         main = "Gene Expression Z-scores Across NAC Treatment Groups")

# Close the PDF device
dev.off()


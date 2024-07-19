
library(GeneNMF)
library(Seurat)
library(ggplot2)
library(UCell)
library(patchwork)
library(Matrix)
library(RcppML)
library(viridis)
library(Seurat)
library(viridis)
library(ggplot2)
library(patchwork)

print("reading data")

st <- readRDS("/Users/akhaliq/Desktop/Gene_NMf/subset.rds")
st1 <- SCTransform(st, assay = "Spatial", verbose = TRUE, method = "poisson")

print("reading data")

ndim <- 15
seu <- FindVariableFeatures(st1, nfeatures = 1000)


print("reading data")
# Check the variable features identified

variable_genes <- VariableFeatures(seu)
print(variable_genes)

seu <- runNMF(seu, k = ndim, assay = "SCT", hvg = variable_genes)
seu <- RunUMAP(seu, reduction = "NMF", dims=1:ndim, reduction.name = "NMF_UMAP", reduction.key = "nmfUMAP_")



DimPlot(seu, reduction = "NMF_UMAP", group.by = "orig.ident", label=T) + theme(aspect.ratio = 1,
axis.text = element_blank(),axis.title = element_blank(),axis.ticks = element_blank()) + ggtitle("NMF UMAP") + NoLegend()

seu.list <- SplitObject(seu, split.by = "orig.ident")

geneNMF.programs <- multiNMF(seu.list, assay="SCT", slot="data", k=4:5, L1=c(0,0),
                    do_centering=TRUE, nfeatures = 2000)

geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs,
                                        nprograms=7,
                                        max.genes=200,
                                        hclust.method="ward.D2",
                                        min.confidence=0.1)

pdf("Plot_meta.pdf", width = 8, height = 6)
ph <- plotMetaPrograms(geneNMF.metaprograms, jaccard.cutoff = c(0,0.8))
dev.off()


######### Hirarchial Clustering for MP and Spearman Coorelation

library(Seurat)
library(pheatmap)
library(RColorBrewer)
library(dendextend)
library(grid)
library(gridExtra)

# Extract MP scores
mp_scores <- seu@meta.data[, c("MP1", "MP2", "MP3", "MP4", "MP5", "MP6", "MP7")]

# Calculate Spearman correlations between MPs
cor_matrix <- cor(mp_scores, method = "spearman")

# Perform hierarchical clustering
hc <- hclust(as.dist(1 - cor_matrix), method = "complete")

# Create annotation for CompositionCluster_CC based on the provided data
cc_data <- data.frame(
  MP = c("MP1", "MP2", "MP3", "MP4", "MP5", "MP6", "MP7"),
  CC07 = c(22006, 1391, 526, 313, 21, 3816, 108),
  CC10 = c(20849, 1650, 0, 1069, 0, 6, 0)
)

cc_annotation <- data.frame(CC = ifelse(cc_data$CC07 > cc_data$CC10, "CC07", "CC10"))
rownames(cc_annotation) <- cc_data$MP

# Ensure CC is a factor with both levels
cc_annotation$CC <- factor(cc_annotation$CC, levels = c("CC07", "CC10"))

# Create color palette for annotation and correlation
annotation_colors <- list(CC = c("CC07" = "#4DAF4A", "CC10" = "#E41A1C"))
correlation_colors <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#F7F7F7", "#F4A582", "#D6604D", "#B2182B", "#67001F"))(100)

# Create and save the heatmap as PDF for better quality
pdf("MP_correlation_heatmap_with_clustering.pdf", width = 10, height = 10)

# Create the main heatmap
main_heatmap <- pheatmap(cor_matrix,
                         color = correlation_colors,
                         cluster_rows = hc,
                         cluster_cols = hc,
                         annotation_col = cc_annotation,
                         annotation_row = cc_annotation,
                         annotation_colors = annotation_colors,
                         cellwidth = 30,
                         cellheight = 30,
                         fontsize = 12,
                         fontsize_number = 10,
                         number_color = "black",
                         display_numbers = TRUE,
                         number_format = "%.2f",
                         border_color = "black",
                         legend = TRUE,
                         annotation_legend = TRUE,
                         annotation_names_col = FALSE,
                         annotation_names_row = FALSE,
                         main = "Spearman Correlation of Metaprograms (MPs)\nwith CompositionCluster_CC Annotation and Hierarchical Clustering",
                         fontsize_row = 12,
                         fontsize_col = 12,
                         treeheight_row = 50,
                         treeheight_col = 50)

# Print the heatmap
print(main_heatmap)

dev.off()

print("Enhanced heatmap with clustering has been saved as 'MP_correlation_heatmap_with_clustering.pdf'")

#########



#### Spearman Coorelation to check the Treated MPs in CC7/10


library(Seurat)
library(pheatmap)
library(dplyr)
library(tibble)
library(RColorBrewer)
library(gridExtra)

# Assuming 'seu' is your Seurat object
# Extract the MP scores and treatment information
mp_data <- seu@meta.data[, c("MP1", "MP2", "MP3", "MP4", "MP5", "MP6", "MP7", "treated")]

# Calculate Spearman correlation
cor_matrix <- cor(mp_data[, 1:7], method = "spearman")

# Calculate proportion of treated cells for each MP
mp_treatment <- mp_data %>%
  group_by(treated) %>%
  summarise(across(starts_with("MP"), mean)) %>%
  filter(treated == "Yes") %>%
  select(-treated) %>%
  as.matrix()

# Transpose mp_treatment to make it a column
mp_treatment <- t(mp_treatment)

# Create annotation data frame
annotation_df <- data.frame(
  Treated_Proportion = mp_treatment[,1]
)
rownames(annotation_df) <- rownames(cor_matrix)

# Create color palettes
correlation_colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
annotation_colors <- list(
  Treated_Proportion = colorRampPalette(c("#3A539B", "#FFFFFF", "#CB4335"))(100)
)

# Create and save the heatmap as PDF
pdf("MP_correlation_heatmap_with_treatment.pdf", width = 10, height = 9)

# Create the heatmap
heatmap <- pheatmap(cor_matrix,
                    color = correlation_colors,
                    cluster_rows = TRUE,
                    cluster_cols = TRUE,
                    cellwidth = 40,
                    cellheight = 40,
                    fontsize = 12,
                    fontsize_number = 10,
                    number_color = "black",
                    display_numbers = TRUE,
                    number_format = "%.2f",
                    border_color = "white",
                    annotation_row = annotation_df,
                    annotation_colors = annotation_colors,
                    annotation_names_row = TRUE,
                    annotation_legend = TRUE,
                    legend = TRUE,
                    main = "Spearman Correlation of Metaprograms (MPs)",
                    fontsize_row = 12,
                    fontsize_col = 12,
                    angle_col = 45,
                    treeheight_row = 20,
                    treeheight_col = 20)

# Add a title for the annotation bar
annotation_title <- textGrob("Proportion\nof Treated\nCells", gp = gpar(fontsize = 10, fontface = "bold"))
heatmap_grob <- grid.arrange(annotation_title, heatmap$gtable, ncol = 2, widths = c(0.1, 0.9))

# Print the arranged heatmap
print(heatmap_grob)

dev.off()

print("Enhanced heatmap of MP Spearman correlations with treatment annotation has been saved as 'MP_correlation_heatmap_with_treatment.pdf'")
###########

###
######### Heatmap, BArplots, Donut


library(Seurat)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(scales)

# Assuming 'seu' is your Seurat object
mp_data <- seu@meta.data[, c("MP1", "MP2", "MP3", "MP4", "MP5", "MP6", "MP7", "treated", "assigned_MP")]

# Calculate Spearman correlation
cor_matrix <- cor(mp_data[, 1:7], method = "spearman")

# Calculate proportion of treated cells for each MP
mp_treatment <- mp_data %>%
  group_by(assigned_MP) %>%
  summarise(Treated = mean(treated == "Yes") * 100,
            Untreated = mean(treated == "No") * 100) %>%
  arrange(match(assigned_MP, paste0("MP", 1:7)))

# Create annotation data frame for heatmap
annotation_df <- data.frame(
  Treated_Proportion = mp_treatment$Treated / 100
)
rownames(annotation_df) <- mp_treatment$assigned_MP

# Create color palettes
correlation_colors <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#F7F7F7", "#F4A582", "#D6604D", "#B2182B", "#67001F"))(100)
annotation_colors <- list(
  Treated_Proportion = colorRampPalette(c("#FEF0D9", "#FD8D3C", "#B30000"))(100)
)

# Create the heatmap
heatmap <- pheatmap(cor_matrix,
                    color = correlation_colors,
                    cluster_rows = TRUE,
                    cluster_cols = TRUE,
                    cellwidth = 40,
                    cellheight = 40,
                    fontsize = 14,
                    fontsize_number = 12,
                    number_color = "black",
                    display_numbers = TRUE,
                    number_format = "%.2f",
                    border_color = "white",
                    annotation_row = annotation_df,
                    annotation_colors = annotation_colors,
                    annotation_names_row = TRUE,
                    annotation_legend = TRUE,
                    legend = TRUE,
                    main = "Spearman Correlation of Metaprograms (MPs)",
                    silent = TRUE)

# Extract the dendrogram order
mp_order <- rownames(cor_matrix)[heatmap$tree_row$order]

# Reorder mp_treatment based on heatmap clustering
mp_treatment <- mp_treatment %>%
  mutate(assigned_MP = factor(assigned_MP, levels = mp_order)) %>%
  arrange(assigned_MP)

# Create the enhanced bar plot
bar_plot <- ggplot(mp_treatment, aes(x = assigned_MP)) +
  geom_bar(aes(y = 100, fill = "Untreated"), stat = "identity", width = 0.7) +
  geom_bar(aes(y = Treated, fill = "Treated"), stat = "identity", width = 0.7) +
  geom_text(aes(y = Treated/2, label = sprintf("%.1f%%", Treated)), 
            color = "white", fontface = "bold", size = 4) +
  geom_text(aes(y = Treated + (100-Treated)/2, label = sprintf("%.1f%%", Untreated)), 
            color = "white", fontface = "bold", size = 4) +
  scale_fill_manual(values = c("Untreated" = "#4DAF4A", "Treated" = "#E41A1C")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "top",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    plot.margin = unit(c(0.5,0.5,0.5,0), "cm"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  ) +
  labs(x = "Metaprogram", y = "Percentage", fill = "Treatment Status", 
       title = "Treatment Distribution Across Metaprograms") +
  coord_flip()

# Create enhanced donut plots
create_donut_plot <- function(mp, treated_pct, untreated_pct) {
  ggplot(data.frame(
    category = c("Treated", "Untreated"),
    value = c(treated_pct, untreated_pct)
  ), aes(x = 2, y = value, fill = category)) +
    geom_bar(stat = "identity", width = 0.7) +
    coord_polar(theta = "y", start = 0) +
    scale_fill_manual(values = c("Treated" = "#E41A1C", "Untreated" = "#4DAF4A")) +
    theme_void() +
    theme(legend.position = "none") +
    geom_text(aes(label = sprintf("%.1f%%", value)), 
              position = position_stack(vjust = 0.5), color = "white", fontface = "bold", size = 4) +
    ggtitle(mp) +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
    xlim(0.5, 2.5)  # This creates the inner hollow
}

donut_plots <- lapply(1:nrow(mp_treatment), function(i) {
  create_donut_plot(mp_treatment$assigned_MP[i], mp_treatment$Treated[i], mp_treatment$Untreated[i])
})

# Create a smaller legend for the Treated_Proportion color scale
treated_prop_legend <- ggplot(data.frame(x = 1:100, y = 1, z = 1:100), aes(x, y, fill = z)) +
  geom_tile() +
  scale_fill_gradientn(colors = annotation_colors$Treated_Proportion,
                       labels = c("0%", "50%", "100%"),
                       breaks = c(1, 50, 100),
                       limits = c(1, 100)) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.key.width = unit(2, "cm"),
        legend.key.height = unit(0.5, "cm"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12, face = "bold")) +
  labs(fill = "Treated Proportion")
# Arrange the plots
pdf("MP_Analysis_Revised_Final_treated_Donut_BArplot_Sperman.pdf", width = 20, height = 15)
grid.arrange(
  grob = heatmap$gtable,
  arrangeGrob(
    bar_plot,
    arrangeGrob(grobs = donut_plots, ncol = 4),
    ncol = 1,
    heights = c(2, 1)
  ),
  treated_prop_legend,
  ncol = 2,
  widths = c(3, 2),
  heights = c(10, 0.5),
  top = textGrob("Metaprogram Correlations and Treatment Proportions", 
                 gp = gpar(fontsize = 24, fontface = "bold")),
  bottom = textGrob("Figure 1: (A) Heatmap: Spearman correlations between Metaprograms (MPs). Side bar shows proportion of treated cells. 
                    (B) Bar plot: Percentage distribution of treated and untreated cells for each MP. 
                    (C) Donut charts: Compact view of treatment proportions for individual MPs.
                    Color scale at bottom represents the treated proportion, ranging from 0% (light orange, Non Treated) to 100% (dark red, Treated).", 
                    gp = gpar(fontsize = 12), just = "left", x = 0.02)
)
dev.off()

print("Revised plot has been saved as 'MP_Analysis_Revised_Final_treated_Donut_BArplot_Sperman.pdf'")


##########
# PATHWAY ANALYSIS 
##########



library(msigdbr)
library(fgsea)

as.data.frame(msigdbr_collections()) # check what DBs you want to use

C5 --> GO:BP --> go BP
H --> Hallmark
C2 --> CP:REACTOME --> Reactome
C2 --> CP:KEGG --> KEGG


top_p <- lapply(geneNMF.metaprograms$metaprograms.genes, function(program) {
  runGSEA(program, universe=rownames(seu), category = "C2", subcategory = "CP:KEGG")
})


# Load required library
library(writexl)

# Assuming top_p is already loaded in your R environment
# If not, load it first (you may need to adjust the file path)
# load("path/to/your/top_p.RData")

# Rename the elements of top_p to ensure correct sheet naming
names(top_p) <- paste0("MP", 1:length(top_p))

# Create an Excel file with separate sheets for each MP
write_xlsx(top_p, path = "MP_KEGG.xlsx")

print("Excel file 'MP_data.xlsx' has been created with separate sheets for each MP.")

# Optional: Print out the structure of the first MP to verify
print(str(top_p[[1]]))


## plots for HAllmark 

# Install and load required packages
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("forcats", quietly = TRUE)) install.packages("forcats")
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
if (!requireNamespace("scales", quietly = TRUE)) install.packages("scales")
if (!requireNamespace("ggtext", quietly = TRUE)) install.packages("ggtext")

library(ggplot2)
library(dplyr)
library(forcats)
library(patchwork)
library(scales)
library(ggtext)

# Function to create a more beautiful bar plot for each MP
create_beautiful_bar_plot <- function(df, mp_name) {
  df <- df %>%
    arrange(padj) %>%
    slice_head(n = 15) %>%  # Top 15 pathways
    mutate(
      pathway = fct_reorder(pathway, padj),
      log_padj = -log10(padj),
      overlap_ratio = overlap / size,
      label = sprintf("%s", pathway),
      info = sprintf("Overlap: %d/%d (%.1f%%)", overlap, size, overlap_ratio*100)
    )
  
  ggplot(df, aes(x = log_padj, y = label)) +
    geom_bar(stat = "identity", aes(fill = overlap_ratio), width = 0.7) +
    geom_text(aes(label = sprintf("%.2f", log_padj)), hjust = -0.1, size = 3, fontface = "bold") +
    geom_text(aes(label = info, x = 0), hjust = 1.1, size = 2.5, color = "darkgrey") +
    scale_fill_gradientn(colours = c("#FDE725FF", "#B4DE2CFF", "#6DCD59FF", "#35B779FF", "#1F9E89FF", "#26828EFF", "#31688EFF", "#3E4989FF", "#482878FF", "#440154FF"),
                         name = "Overlap\nRatio", labels = percent) +
    labs(
      title = mp_name,
      subtitle = "Top 15 Enriched Pathways",
      x = "-log10(Adjusted p-value)",
      y = NULL,
      caption = "Bar length represents significance. Color intensity shows overlap ratio."
    ) +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = "ghostwhite", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      axis.text.y = element_markdown(size = 7, color = "black"),
      axis.text.x = element_text(size = 8, color = "black"),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold", color = "navy"),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "darkblue"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_line(color = "lightgrey", linetype = "dashed"),
      legend.position = "right",
      legend.key.size = unit(0.5, "cm"),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 8),
      plot.caption = element_text(size = 8, color = "darkgrey", hjust = 1),
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt")
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.2)))
}

# Create plots for each MP
plot_list <- lapply(names(top_p), function(mp_name) {
  create_beautiful_bar_plot(top_p[[mp_name]], mp_name)
})

# Combine plots using patchwork
combined_plot <- wrap_plots(plot_list, ncol = 2) +
  plot_annotation(
    title = 'Comprehensive Pathway Enrichment Analysis Across MPs',
    subtitle = 'Visualizing Significance, Overlap, and Set Size',
    theme = theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5, color = "navy"),
      plot.subtitle = element_text(size = 16, hjust = 0.5, color = "darkblue"),
      plot.background = element_rect(fill = "ghostwhite", color = NA)
    )
  ) &
  theme(plot.background = element_rect(fill = "ghostwhite", color = NA))

# Save the combined plot as a PDF
ggsave("pathway_enrichment_all_MPs_KEGG.pdf", combined_plot, width = 18, height = 10 * ceiling(length(plot_list)/2), limitsize = FALSE)

print("Enhanced bar plots for all MPs have been combined and saved as 'pathway_enrichment_all_MPs_Hallmark.pdf'.")

####



##########
Coorelation B/w Regav Programs and Our mp_scores

##########


# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# Read the data
tumor_core <- read.csv("metaprograms_genes_TumorCore.csv")
regav <- read.csv("output_data_Regav.csv")

# Function to create a binary matrix for genes in metaprograms
create_binary_matrix <- function(data, gene_col, mp_col) {
  data %>%
    distinct({{gene_col}}, {{mp_col}}) %>%
    mutate(present = 1) %>%
    pivot_wider(names_from = {{mp_col}}, values_from = present, values_fill = 0)
}

# Create binary matrices
tumor_matrix <- create_binary_matrix(tumor_core, Genes, Meta_programs)
regav_matrix <- create_binary_matrix(regav, Genes, Regvs_programs)

# Get common genes
common_genes <- intersect(tumor_matrix$Genes, regav_matrix$Genes)

# Filter matrices to include only common genes
tumor_filtered <- tumor_matrix %>% filter(Genes %in% common_genes) %>% select(-Genes)
regav_filtered <- regav_matrix %>% filter(Genes %in% common_genes) %>% select(-Genes)

# Calculate correlation
correlation_matrix <- cor(tumor_filtered, regav_filtered)

# Write correlation matrix to CSV
write.csv(correlation_matrix, "metaprogram_correlation_matrix.csv")

# Enhanced Heatmap
correlation_df <- as.data.frame(correlation_matrix) %>%
  rownames_to_column("Tumor_MP") %>%
  pivot_longer(-Tumor_MP, names_to = "Regav_MP", values_to = "Correlation")

# Create custom color palette
color_palette <- colorRampPalette(c("#4575B4", "white", "#D73027"))(100)

ggplot(correlation_df, aes(x = Regav_MP, y = Tumor_MP, fill = Correlation)) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_gradientn(colors = color_palette, limits = c(-1, 1), 
                       breaks = seq(-1, 1, by = 0.5), labels = scales::number_format(accuracy = 0.1)) +
  geom_text(aes(label = sprintf("%.2f", Correlation)), 
            size = 3, color = ifelse(abs(correlation_df$Correlation) > 0.5, "white", "black")) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  labs(
    title = "Correlation Heatmap",
    subtitle = "Tumor_core vs Regav Metaprograms",
    fill = "Correlation"
  )

ggsave("correlation_heatmap_enhanced.png", width = 16, height = 12, dpi = 300, bg = "white")

###############

Network Plots b/w regav and ours 

#########

# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(igraph)
library(ggraph)
library(tidygraph)

# (Keep the data loading and correlation calculation part from the previous script)

# Enhanced Heatmap
# (Keep the heatmap code from the previous response)

# Enhanced Network plot
network_df <- correlation_df %>%
  filter(abs(Correlation) > 0.3) %>%
  select(Tumor_MP, Regav_MP, Correlation)

g <- as_tbl_graph(network_df, directed = FALSE) %>%
  activate(nodes) %>%
  mutate(type = ifelse(name %in% unique(network_df$Tumor_MP), "Tumor", "Regav"))

set.seed(42)  # for reproducibility
layout <- create_layout(g, layout = "fr")

ggraph(layout) +
  geom_edge_link(aes(edge_alpha = abs(Correlation), color = Correlation), 
                 edge_width = 1, show.legend = TRUE) +
  geom_node_point(aes(color = type), size = 10, alpha = 0.8) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3, fontface = "bold") +
  scale_edge_colour_gradientn(colours = colorRampPalette(c("#4575B4", "white", "#D73027"))(100),
                              limits = c(-1, 1), 
                              breaks = seq(-1, 1, by = 0.5), 
                              labels = scales::number_format(accuracy = 0.1)) +
  scale_color_manual(values = c("Tumor" = "#FF9933", "Regav" = "#33CCCC")) +
  scale_edge_alpha(range = c(0.3, 1)) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.background = element_rect(fill = "white", color = NA),
    legend.box.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 10)),
    plot.subtitle = element_text(hjust = 0.5, size = 12, margin = margin(b = 20)),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  labs(
    title = "Metaprogram Correlation Network",
    subtitle = "Strong correlations (|r| > 0.3) between Tumor_core and Regav metaprograms",
    color = "Dataset",
    edge_colour = "Correlation",
    edge_alpha = "Correlation Strength"
  )

ggsave("correlation_network_enhanced.png", width = 16, height = 12, dpi = 300, bg = "white")

## interpretation of network pltos

This network plot visualizes the correlations between metaprograms from two datasets: Tumor_core and Regav. Here's how to interpret it:

Nodes: Each circle represents a metaprogram. The color indicates which dataset it belongs to:

Orange nodes are from the Tumor_core dataset (MP1, MP2, MP4, MP5, MP7)
Blue nodes are from the Regav dataset (Squamoid, Neuroendocrine-like, Acinar-like, Classical-like, Ribosomal, signaling, TNF-NFkB, progenitor)


Edges (Lines): The lines connecting the nodes represent correlations between metaprograms.

The presence of a line indicates a strong correlation (|r| > 0.3)
The color of the line represents the direction and strength of the correlation:

Red lines indicate positive correlations
Blue lines would indicate negative correlations (not visible in this plot)


The thickness of the line represents the strength of the correlation (thicker = stronger)


Layout: The spatial arrangement of nodes is based on the strength of their connections. Closely related metaprograms are positioned nearer to each other.
Key observations:

MP1 seems to be a central node, showing strong correlations with several Regav metaprograms (Classical-like, Ribosomal, signaling)
MP2 has a strong positive correlation with TNF-NFkB
MP3 shows a positive correlation with Acinar-like
MP7 is correlated with the progenitor metaprogram
MP4 and MP5 show weaker correlations with Squamoid and Neuroendocrine-like respectively


Interpretation:

Metaprograms connected by lines likely share similar gene expression patterns or biological functions across the two datasets
The strongest correlations (thickest red lines) suggest the most significant similarities between metaprograms
Isolated nodes or those with few connections may represent unique or dataset-specific metaprograms



This visualization helps identify relationships between metaprograms across different datasets, potentially revealing shared biological processes or regulatory mechanisms.

#####


1_Plot_meta.pdf
10_correlation_corrplot.png
9_correlation_heatmap_enhanced.png
8_pathway_enrichment_all_MPs_KEGG.pdf
7_pathway_enrichment_all_MPs_Go_BP.pdf
9_pathway_enrichment_all_MPs_REACTOME.pdf
6_pathway_enrichment_all_MPs_Hallmark.pdf
5_MP_Analysis_Revised_Final_treated_Donut_BArplot_Sperman.pdf
4_MP_correlation_heatmap_with_clustering.pdf
11_hypergeometric_test_plot.png
3_split_umap_by_metaprogram.pdf
2_feature_plot2.pdf
####

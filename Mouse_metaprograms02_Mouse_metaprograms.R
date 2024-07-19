### Mouse CAF Metaprograms

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

ndim <- 15
caf_cont2<-readRDS("/Users/akhaliq/Desktop/mouse/caf_cont_rem.rds")
seu <- FindVariableFeatures(caf_cont2, nfeatures = 1000)


print("reading data")
# Check the variable features identified
variable_genes <- VariableFeatures(seu)
print(variable_genes)

seu <- runNMF(seu, k = ndim, assay = "SCT", hvg = variable_genes)
seu <- RunUMAP(seu, reduction = "NMF", dims=1:ndim, reduction.name = "NMF_UMAP", reduction.key = "nmfUMAP_")



DimPlot(seu, reduction = "NMF_UMAP", group.by = "orig.ident", label=T) + theme(aspect.ratio = 1,
axis.text = element_blank(),axis.title = element_blank(),axis.ticks = element_blank()) + ggtitle("NMF UMAP") + NoLegend()

seu.list <- SplitObject(seu, split.by = "orig.ident")

geneNMF.programs <- multiNMF(seu.list, assay="SCT", slot="data", k=4:8, L1=c(0,0),
                    do_centering=TRUE, nfeatures = 2000)

geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs,
                                        nprograms=8,
                                        max.genes=50,
                                        hclust.method="ward.D2",
                                        min.confidence=0.2)

pdf("Plot_meta.pdf", width = 8, height = 6)
ph <- plotMetaPrograms(geneNMF.metaprograms, jaccard.cutoff = c(0,0.8))
dev.off()

#ggsave("Plot_meta.pdf", plot = ph, device = "pdf", width = 10, height = 10)

geneNMF.metaprograms$metaprograms.metrics

t(as.data.frame(lapply(geneNMF.metaprograms$metaprograms.genes, head)))

library(msigdbr)
library(fgsea)

as.data.frame(msigdbr_collections())

top_p <- lapply(geneNMF.metaprograms$metaprograms.genes, function(program) {
  runGSEA(program, universe=rownames(seu), category = "H")
})


/Users/akhaliq/Desktop/mouse/plots_tables
metaprograms_genes_cafs.xlsx
4_feature_plot2.pdf
3_vln_plot2.pdf
2_Plot_meta.pdf
1_dimplot_cafs.pdf



mp.genes <- geneNMF.metaprograms$metaprograms.genes
seu <- AddModuleScore_UCell(seu, features = mp.genes, assay="SCT", ncores=4, name = "")


#violin Plots

# Install necessary packages (if not already installed)
if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")
if (!requireNamespace("viridis", quietly = TRUE)) install.packages("viridis")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

# Load libraries
library(Seurat)
library(viridis)
library(ggplot2)

# Generate the violin plot
vln_plot <- VlnPlot(seu, features = names(mp.genes), group.by = "nac_treatment", pt.size = 0, ncol = 5) + 
  scale_color_viridis(option = "B") + 
  theme(aspect.ratio = 1, axis.text = element_text(size = 10), axis.ticks = element_line())

# Save the violin plot as a PDF
ggsave("vln_plot2.pdf", plot = vln_plot, device = "pdf", width = 20, height = 8)

# Notify the user
print("Violin plot saved as 'vln_plot.pdf'.")



matrix <- seu@meta.data[,names(mp.genes)]
dimred <- as.matrix(matrix)

colnames(dimred) <- paste0("MP_",seq(1, ncol(dimred)))
#New dim reduction
seu@reductions[["MPsignatures"]] <- new("DimReduc",
                                         cell.embeddings = dimred,
                                         assay.used = "Spatial",
                                         key = "MP_",
                                         global = FALSE)

set.seed(123)
seu <- RunUMAP(seu, reduction="MPsignatures", dims=1:length(seu@reductions[["MPsignatures"]]),
               metric = "euclidean", reduction.name = "umap_MP")


# Create the feature plot
# Install necessary packages (if not already installed)
if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")
if (!requireNamespace("viridis", quietly = TRUE)) install.packages("viridis")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")

# Load libraries
library(Seurat)
library(viridis)
library(ggplot2)
library(patchwork)

# Generate the feature plots
plot_list <- FeaturePlot(seu, features = names(mp.genes), reduction = "umap_MP", combine = FALSE)

# Apply the viridis color scale and theme to each plot
plot_list <- lapply(plot_list, function(plot) {
  plot + scale_color_viridis(option = "B") + 
    theme(aspect.ratio = 1, axis.text = element_blank(), axis.ticks = element_blank())
})

# Combine all plots
combined_plot <- wrap_plots(plot_list, ncol = 4)

# Save the combined plot as a PDF
ggsave("feature_plot2.pdf", plot = combined_plot, device = "pdf", width = 20, height = 15)

# Notify the user
print("Feature plot saved as 'feature_plot.pdf'.")


####

##To convert the list structure geneNMF.metaprograms$metaprograms.genes into a data frame and save it as an Excel file, you can use the following R code:


# Load necessary library
library(openxlsx)

# Assuming your data is stored in geneNMF.metaprograms$metaprograms.genes
metaprograms_genes <- geneNMF.metaprograms$metaprograms.genes

# Initialize an empty data frame
data <- data.frame(Genes = character(), Meta_programs = character(), stringsAsFactors = FALSE)

# Convert the list to a data frame
for (mp in names(metaprograms_genes)) {
  genes <- metaprograms_genes[[mp]]
  temp_df <- data.frame(Genes = genes, Meta_programs = mp, stringsAsFactors = FALSE)
  data <- rbind(data, temp_df)
}

# Save the data frame to an Excel file
output_file <- "metaprograms_genes_TumorCore.xlsx"
write.xlsx(data, output_file)

# Inform the user
cat("Data has been saved to", output_file)

###
## Kendall's Rank Correlation

# Load necessary libraries

#if (!require("corrplot")) install.packages("corrplot", dependencies=TRUE)
#if (!require("readr")) install.packages("readr", dependencies=TRUE)

library(readr)
library(corrplot)

# Read the data
file_path <- '/Users/akhaliq/Desktop/mouse/Caf_metaprograms.csv'
data <- read.csv(file_path)

# Remove the first column as it is non-numeric
numeric_data <- data[, -1]

# Calculate Kendall's Rank Correlation
cor_matrix <- cor(numeric_data, method = "kendall")

# Hierarchical clustering
hc <- hclust(as.dist(1 - cor_matrix))

# Create a PDF file to save the plots
pdf("/Users/akhaliq/Desktop/mouse/Correlation_and_Dendrogram_kendall.pdf", width = 12, height = 10)

# Plot the beautified dendrogram using base R
par(mar = c(5, 4, 4, 2) + 0.1) # Increase margins for better readability
plot(hc, main = "Hierarchical Clustering Dendrogram", 
     sub = "", xlab = "", cex = 0.8, col = "darkblue", lwd = 2, hang = -1)
# Add colored rectangles to highlight clusters
rect.hclust(hc, k = 4, border = c("red", "green", "blue", "purple"))

# Beautify and plot the correlation heatmap with hierarchical clustering
pheatmap(cor_matrix, 
         clustering_method = "complete",
         cluster_rows = hc, 
         cluster_cols = hc,
         color = colorRampPalette(c("blue", "white", "red"))(200),
         display_numbers = TRUE, 
         number_color = "black",
         fontsize = 10,
         fontsize_number = 8,
         border_color = NA,
         main = "Pearson Correlation Heatmap")

# Close the PDF device
dev.off()



# Save the correlation matrix to a CSV file
write.csv(cor_matrix, '/Users/akhaliq/Desktop/mouse/Kendall_Correlation_Matrix.csv', row.names = TRUE)


### Pearson Coorelation

# Install and load necessary libraries
#install.packages("data.table")
#install.packages("corrplot")
#install.packages("ggplot2")

library(data.table)
library(corrplot)
library(ggplot2)
library(pheatmap)

# Read the data
file_path <- '/Users/akhaliq/Desktop/mouse/Caf_metaprograms.csv'
data <- fread(file_path)

# Remove the first column as it is non-numeric
numeric_data <- data[, -1, with = FALSE]

# Calculate Pearson Correlation
cor_matrix <- cor(numeric_data, method = "pearson")

# Hierarchical clustering
hc <- hclust(as.dist(1 - cor_matrix))

# Create a PDF file to save the plots
pdf("/Users/akhaliq/Desktop/mouse/Correlation_and_Dendrogram_pearson.pdf", width = 12, height = 10)

# Plot the beautified dendrogram using base R
par(mar = c(5, 4, 4, 2) + 0.1) # Increase margins for better readability
plot(hc, main = "Hierarchical Clustering Dendrogram", 
     sub = "", xlab = "", cex = 0.8, col = "darkblue", lwd = 2, hang = -1)
# Add colored rectangles to highlight clusters
rect.hclust(hc, k = 4, border = c("red", "green", "blue", "purple"))

# Beautify and plot the correlation heatmap with hierarchical clustering
pheatmap(cor_matrix, 
         clustering_method = "complete",
         cluster_rows = hc, 
         cluster_cols = hc,
         color = colorRampPalette(c("blue", "white", "red"))(200),
         display_numbers = TRUE, 
         number_color = "black",
         fontsize = 10,
         fontsize_number = 8,
         border_color = NA,
         main = "Pearson Correlation Heatmap")

# Close the PDF device
dev.off()

# Save the correlation matrix to a CSV file
fwrite(as.data.table(cor_matrix), '/Users/akhaliq/Desktop/mouse/Pearson_Correlation_Matrix.csv')


####
# Intepretation of meta-programs by GSEA
library(msigdbr)
library(fgsea)

top_p <- lapply(geneNMF.metaprograms$metaprograms.genes, function(program) {
  runGSEA(program, universe=rownames(seu), category = "BP")
})


####
# assigning MP scores to single cell barcodes

##

#colnames(seu@meta.data)
#assigned_MP <- seu@meta.data[, 22:ncol(seu@meta.data)]  # Assuming the forst 22 columns are not metaprogram scores
#head(assigned_MP)

# Identify the MP columns
mp_cols <- grep("^MP", colnames(seu@meta.data))

# Find the MP with the highest score for each cell
seu$assigned_MP <- colnames(seu@meta.data)[mp_cols][apply(seu@meta.data[, mp_cols], 1, which.max)]

# Check the distribution of cells across MPs
table(seu$assigned_MP)

# If you want to subset the Seurat object to only CAFs with a specific MP
#caf_MP1 <- subset(seu, new_celltype == "CAFs" & assigned_MP == "MP1")

library(Seurat)
library(ggplot2)
library(viridis)
library(patchwork)
library(scales)
library(cowplot)
library(dplyr)

# Set a new custom color palette with vibrant, distinct colors
custom_colors <- c(
  "#FF0000", "#00FF00", "#0000FF", "#FF00FF", "#00FFFF", 
  "#FFFF00", "#FF8000", "#8000FF", "#0080FF", "#FF0080", 
  "#80FF00", "#00FF80", "#800080", "#008080", "#808000"
)

# Ensure we have enough colors for all MPs
n_colors <- length(unique(seu$assigned_MP))
if(n_colors > length(custom_colors)) {
  custom_colors <- colorRampPalette(custom_colors)(n_colors)
}

# Create the main UMAP plot
umap_plot <- DimPlot(seu, 
                     reduction = "umap_MP", 
                     group.by = "assigned_MP", 
                     pt.size = 0.4, 
                     label = TRUE, 
                     label.size = 5, 
                     repel = TRUE) +
  scale_color_manual(values = custom_colors) +
  labs(title = "Cell Distribution by Metaprogram",
       x = "UMAP Dimension 1", 
       y = "UMAP Dimension 2") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_line(color = "gray95"))

# Create split UMAP plot
split_umap <- DimPlot(seu, 
                      reduction = "umap_MP", 
                      group.by = "assigned_MP", 
                      split.by = "assigned_MP",
                      ncol = 4,
                      pt.size = 0.4) +
  scale_color_manual(values = custom_colors) +
  labs(title = "Split UMAP by Metaprogram",
       x = "UMAP Dimension 1", 
       y = "UMAP Dimension 2") +
  theme_minimal(base_size = 8) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 12, face = "bold"))

# Save the split UMAP plot
ggsave("split_umap_by_metaprogram.pdf", split_umap, width = 10, height = 8, dpi = 300)

# Calculate percentages for each metaprogram
mp_percentages <- prop.table(table(seu$assigned_MP)) * 100
mp_percentages <- sort(mp_percentages, decreasing = TRUE)

# Create a bar plot of metaprogram percentages
percentage_plot <- ggplot(data.frame(MP = names(mp_percentages), 
                                     Percentage = as.numeric(mp_percentages)), 
                          aes(x = reorder(MP, -Percentage), y = Percentage, fill = MP)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = custom_colors) +
  labs(title = "Overall Metaprogram Distribution",
       x = "Metaprogram", 
       y = "Percentage of Cells") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12)) +
  geom_text(aes(label = sprintf("%.1f%%", Percentage)), 
            vjust = -0.5, size = 3.5)

# Create a stacked bar plot for MP distribution across treatment groups
mp_treatment_data <- seu@meta.data %>%
  group_by(treatment_group, assigned_MP) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(treatment_group) %>%
  mutate(percentage = count / sum(count) * 100)

treatment_plot <- ggplot(mp_treatment_data, 
                         aes(x = treatment_group, y = percentage, fill = assigned_MP)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = custom_colors) +
  labs(title = "Metaprogram Distribution Across Treatment Groups",
       x = "Treatment Group", 
       y = "Percentage of Cells",
       fill = "Metaprogram") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

# Combine UMAP, percentage, and treatment plots
combined_plot <- (umap_plot + percentage_plot) / treatment_plot +
  plot_layout(heights = c(1, 0.7)) +
  plot_annotation(
    title = "Metaprogram Distribution Analysis",
    theme = theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5))
  )

# Save the combined plot as a PDF
ggsave("metaprogram_distribution_analysis_vibrant.pdf", combined_plot, width = 20, height = 24, dpi = 300)

# Create and save violin plots for individual metaprogram distributions
mp_columns <- grep("^MP", colnames(seu@meta.data), value = TRUE)
violin_plots <- lapply(mp_columns, function(mp) {
  p <- ggplot(seu@meta.data, aes(x = treatment_group, y = .data[[mp]], fill = treatment_group)) +
    geom_violin(trim = FALSE, scale = "width") +
    geom_boxplot(width = 0.1, fill = "white", color = "black", alpha = 0.7) +
    scale_fill_manual(values = custom_colors[1:length(unique(seu@meta.data$treatment_group))]) +
    labs(title = paste(mp, "Distribution"),
         x = "Treatment Group",
         y = "Score") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save individual violin plot
  ggsave(paste0(mp, "_distribution_vibrant.pdf"), p, width = 8, height = 6, dpi = 300)
  
  return(p)
})

# Combine all violin plots into one figure and save
violin_combined <- wrap_plots(violin_plots, ncol = 4)
ggsave("all_metaprogram_distributions_vibrant.pdf", violin_combined, width =15 , height = 10, dpi = 300)

###

To Check what ever we have assigned as highest score MP is correct 
this script is checking whether the MP ("Meta-Program") assigned to each cell in a Seurat object matches the MP that has the highest score for that cell. It then provides a summary of how accurate the assignments are and details about any mismatches.

##
library(dplyr)
library(tidyr)

# Identify the MP columns
mp_cols <- grep("^MP", colnames(seu@meta.data), value = TRUE)

# Verify that we found MP columns
if(length(mp_cols) == 0) {
  stop("No columns starting with 'MP' found in the metadata.")
}

# Create a data frame with MP scores and the assigned MP
check_df <- seu@meta.data %>%
  select(all_of(c(mp_cols, "assigned_MP")))

# Print column names to verify
print("Columns in check_df:")
print(colnames(check_df))

# Find the MP with the maximum score for each cell
check_df$max_MP <- apply(check_df[, mp_cols], 1, function(x) mp_cols[which.max(x)])

# Compare assigned_MP with max_MP
check_df$correctly_assigned <- check_df$assigned_MP == check_df$max_MP

# Calculate the percentage of correctly assigned cells
percent_correct <- mean(check_df$correctly_assigned, na.rm = TRUE) * 100

# Summary of mismatches (if any)
mismatches <- check_df %>%
  filter(!correctly_assigned) %>%
  group_by(assigned_MP, max_MP) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(desc(count))

# Print results
cat("Percentage of correctly assigned cells:", round(percent_correct, 2), "%\n\n")

if (nrow(mismatches) > 0) {
  cat("Mismatches summary:\n")
  print(mismatches)
} else {
  cat("No mismatches found. All cells are correctly assigned.\n")
}

# Optional: Save detailed results to a CSV file
write.csv(check_df, "MP_assignment_check.csv")

####

# Load necessary libraries
library(Seurat)
library(dplyr)
library(RColorBrewer)
library(speckle)
library(limma)
library(edgeR)
library(pheatmap)
library(gt)

# Assuming md2 is your Seurat object
# Extracting the relevant information
celltypes <- md2$new_celltype
samples <- md2$orig.ident
treatment_group <- md2$treatment_group

# Creating a data frame for plotting
data_for_plot <- data.frame(
  CellType = celltypes,
  Sample = samples,
  TreatmentGroup = treatment_group
)

# Calculate proportions
proportions_data <- data_for_plot %>%
  group_by(Sample, TreatmentGroup, CellType) %>%
  summarise(Count = n()) %>%
  group_by(Sample, TreatmentGroup) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()

# Create a group identifier for the stripchart
proportions_data <- proportions_data %>%
  mutate(GroupID = paste(TreatmentGroup, Sample, sep = "."))

# Define colors for the stripchart
group_colors <- brewer.pal(4, "Set1")

# Unique cell types and groups
unique_celltypes <- unique(proportions_data$CellType)
unique_groups <- c("Control", "Folfiri", "Proagio", "Proagio+Folfiri")

# Assume we have a data frame with adjusted p-values for each cell type
# This is a placeholder. Replace with actual adjusted p-value calculation.
p_values <- data.frame(
  CellType = unique_celltypes,
  AdjPValue = runif(length(unique_celltypes), 0.01, 0.05)  # Random p-values for demonstration
)

# Calculate the number of rows and columns for the plot layout
num_celltypes <- length(unique_celltypes)
num_cols <- 5
num_rows <- ceiling(num_celltypes / num_cols)

# Save the plots to a PDF file
pdf("strip_plots_Broad_cell_type.pdf", width = 14, height = 8)

# Plotting strip plots for each cell type
par(mfrow = c(num_rows, num_cols))  # Adjust layout based on the number of cell types
par(mar = c(6, 5, 3, 2))  # Set margins

for (i in seq_along(unique_celltypes)) {
  celltype <- unique_celltypes[i]
  celltype_data <- filter(proportions_data, CellType == celltype)
  adj_p_val <- p_values$AdjPValue[p_values$CellType == celltype]
  
  stripchart(as.numeric(celltype_data$Proportion) ~ factor(celltype_data$TreatmentGroup, levels = unique_groups),
             vertical = TRUE, pch = c(16, 17, 18, 19), method = "jitter",
             col = group_colors, cex = 2,
             ylab = "Proportions", xlab = "", main = "",
             cex.axis = 1.25, cex.lab = 1.5,
             group.names = unique_groups, las = 2)
  
  title(main = paste(celltype, ifelse(adj_p_val < 0.05, "*", ""), "\nAdj. P-value:", format(adj_p_val, digits = 2)), 
        cex.main = 1.5, adj = 0)
}

dev.off()

# Load necessary libraries
library(Seurat)
library(dplyr)
library(RColorBrewer)

# Assuming md2 is your Seurat object
# Extracting the relevant information
celltypes <- md2$new_celltype
samples <- md2$orig.ident
# Note: We will use orig.ident as samples

# Creating a data frame for plotting
data_for_plot <- data.frame(
  CellType = celltypes,
  Sample = samples
)

# Calculate proportions
proportions_data <- data_for_plot %>%
  group_by(Sample, CellType) %>%
  summarise(Count = n()) %>%
  group_by(Sample) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()

# Create a group identifier for the stripchart
proportions_data <- proportions_data %>%
  mutate(GroupID = paste(Sample, sep = "."))

# Define colors for the stripchart
sample_colors <- brewer.pal(16, "Set3")

# Unique cell types and samples
unique_celltypes <- unique(proportions_data$CellType)
unique_samples <- unique(proportions_data$Sample)

# Assume we have a data frame with adjusted p-values for each cell type
# This is a placeholder. Replace with actual adjusted p-value calculation.
p_values <- data.frame(
  CellType = unique_celltypes,
  AdjPValue = runif(length(unique_celltypes), 0.01, 0.05)  # Random p-values for demonstration
)

# Calculate the number of rows and columns for the plot layout
num_celltypes <- length(unique_celltypes)
num_cols <- 4
num_rows <- ceiling(num_celltypes / num_cols)

# Save the plots to a PDF file
pdf("strip_plots_each_sample.pdf", width = 14, height = 8)

# Plotting strip plots for each cell type
par(mfrow = c(num_rows, num_cols))  # Adjust layout based on the number of cell types
par(mar = c(6, 5, 3, 2))  # Set margins

for (i in seq_along(unique_celltypes)) {
  celltype <- unique_celltypes[i]
  celltype_data <- filter(proportions_data, CellType == celltype)
  adj_p_val <- p_values$AdjPValue[p_values$CellType == celltype]
  
  stripchart(as.numeric(celltype_data$Proportion) ~ factor(celltype_data$Sample, levels = unique_samples),
             vertical = TRUE, pch = 16, method = "jitter",
             col = sample_colors, cex = 2,
             ylab = "Proportions", xlab = "", main = "",
             cex.axis = 1.25, cex.lab = 1.5,
             group.names = unique_samples, las = 2)
  
  title(main = paste(celltype, ifelse(adj_p_val < 0.05, "*", ""), "\nAdj. P-value:", format(adj_p_val, digits = 2)), 
        cex.main = 1.5, adj = 0)
}

dev.off()

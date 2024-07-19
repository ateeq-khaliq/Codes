## Aluvial Plot


install.packages("Seurat")
install.packages("ggplot2")
install.packages("ggalluvial")
install.packages("reshape2")
install.packages("ggrepel")
library(Seurat)
library(ggplot2)
library(ggalluvial)
library(reshape2)
library(ggrepel)


# Extract metadata
metadata <- md2@meta.data

# Create a table of cell counts
cell_counts <- table(metadata$new_celltype, metadata$treatment_group)

# Convert the counts to proportions
cell_proportions <- prop.table(cell_counts, 2) * 100

# Convert to a data frame
cell_proportions_df <- as.data.frame(cell_proportions)
colnames(cell_proportions_df) <- c("Cell_Type", "Treatment_Group", "Proportion")

# Reshape the data for plotting
cell_proportions_long <- melt(cell_proportions_df, id.vars = c("Cell_Type", "Treatment_Group"), value.name = "Proportion")


# Define a custom color palette
custom_colors <- c("B Cells" = "#E41A1C", "CAFs" = "#377EB8", "DC" = "#4DAF4A", "Endothelial Cells" = "#984EA3",
                   "Epithelial cells" = "#FF7F00", "Myeloids" = "#FFFF33", "T cells" = "#A65628", "Tumor Cells" = "#F781BF")

# Create the plot
p <- ggplot(cell_proportions_long,
            aes(axis1 = Cell_Type, axis2 = Treatment_Group, y = Proportion)) +
  geom_alluvium(aes(fill = Cell_Type), width = 0.25, alpha = 0.8, color = "gray") +
  geom_stratum(width = 0.35, color = "black", size = 0.5) +
  geom_text_repel(stat = "stratum", aes(label = after_stat(stratum)), size = 4, color = "black", fontface = "bold", direction = "y") +
  scale_fill_manual(values = custom_colors) +
  scale_x_discrete(limits = c("Cell_Type", "Treatment_Group"), expand = c(.2, .2)) +
  labs(title = "Proportions of Cells Across Different Treatment Groups", 
       x = "Cell Type", 
       y = "Proportion (%)",
       caption = "sc-RNAseq Proagio Mouse Analysis: Masood Lab") +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    axis.title.x = element_text(face = "bold", size = 16),
    axis.title.y = element_text(face = "bold", size = 16),
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 12),
    plot.caption = element_text(hjust = 1, size = 12, face = "italic", margin = margin(t = 10))
  )

# Save the plot to a PDF
ggsave("cell_proportions_alluvial_plot.pdf", plot = p, device = "pdf", width = 14, height = 10)

library(Seurat)
library(compositions)
library(tidyverse)
library(clustree)
library(patchwork)
library(uwot)
library(scran)
library(cluster)
library(ggrastr)
library(cowplot)

so <- readRDS("/Users/akhaliq/Desktop/meenakshi/ISCHIA_k10_only_spatial_Assay.rds")
rctd_tier2 <- read.csv("/Users/akhaliq/Downloads/RCTD_norm_weights 1.csv", header=T, row.names=1, sep=",")

metadata <- so@meta.data %>%
  select(orig.ident, CompositionCluster_CC) %>%
  rownames_to_column("row_id")

# Ensure the data is in the right format
rownames(rctd_tier2) <- make.unique(rownames(rctd_tier2))
 
# Log transformation of rctd_tier1
log_comps <- log10(rctd_tier2)

# Prepare data for summary statistics
cluster_summary_pat <- rctd_tier2 %>%
  as.data.frame() %>%
  rownames_to_column("row_id") %>%
  left_join(metadata, by = "row_id") %>%  # Join with meta_data using row_id as the key
  pivot_longer(-c(row_id, orig.ident, patient_id, neoadjuvant_chemo, CompositionCluster_CC), values_to = "ct_prop", names_to = "cell_type") %>%
  group_by(orig.ident, patient_id, neoadjuvant_chemo, CompositionCluster_CC, cell_type) %>%
  summarize(median_ct_prop = median(ct_prop, na.rm = TRUE))


  # Aggregate data for median ct prop
cluster_summary <- cluster_summary_pat %>%
  ungroup() %>%
  group_by(CompositionCluster_CC, cell_type) %>%
  summarize(patient_median_ct_prop = median(median_ct_prop, na.rm = TRUE))


  # Prepare matrix for hierarchical clustering
cluster_summary_mat <- cluster_summary %>%
  pivot_wider(values_from = patient_median_ct_prop, names_from = cell_type, values_fill = list(patient_median_ct_prop = 0)) %>%
  column_to_rownames("CompositionCluster_CC") %>%
  as.matrix()


  # Perform hierarchical clustering
cluster_order <- hclust(dist(cluster_summary_mat))$labels[hclust(dist(cluster_summary_mat))$order] # use this if you want to order the clusters based on the hierarchical clustering
ct_order <- hclust(dist(t(cluster_summary_mat)))$labels[hclust(dist(t(cluster_summary_mat)))$order]
 

 # Wilcoxon test for characteristic cell types
run_wilcox_up <- function(prop_data) {
  prop_data_group <- prop_data[["CompositionCluster_CC"]] %>% unique() %>% set_names()
  map(prop_data_group, function(g) {
    test_data <- prop_data %>%
      mutate(test_group = ifelse(CompositionCluster_CC == g, "target", "rest")) %>%
      mutate(test_group = factor(test_group, levels = c("target", "rest")))
    wilcox.test(median_ct_prop ~ test_group, data = test_data, alternative = "greater") %>%
      broom::tidy()
  }) %>% enframe("CompositionCluster_CC") %>% unnest()
}

wilcoxon_res <- cluster_summary_pat %>%
  ungroup() %>%
  group_by(cell_type) %>%
  nest() %>%
  mutate(wres = map(data, run_wilcox_up)) %>%
  dplyr::select(wres) %>%
  unnest() %>%
  ungroup() %>%
  mutate(p_corr = p.adjust(p.value)) %>%
  mutate(significant = ifelse(p_corr <= 0.15, "*", ""))


  file_path_cluster_summ <- here::here("summary_of_clusters.txt")
  file_path_wilcox_res <- here::here( "wilcoxon_res_cells_clusters.txt")

###
library(ggplot2)
library(dplyr)
library(viridis)
library(ggtext)
library(scales)

mean_ct_prop_plt <- cluster_summary %>%
  left_join(wilcoxon_res, by = c("CompositionCluster_CC", "cell_type")) %>%
  mutate(cell_type = factor(cell_type, levels = ct_order),
         CompositionCluster_CC = factor(CompositionCluster_CC, levels = cluster_order)) %>%
  ungroup() %>%
  group_by(cell_type) %>%
  mutate(scaled_pat_median = (patient_median_ct_prop - mean(patient_median_ct_prop)) / sd(patient_median_ct_prop)) %>%
  ungroup()

max_abs_value <- max(abs(mean_ct_prop_plt$scaled_pat_median), na.rm = TRUE)

ggplot_object <- ggplot(mean_ct_prop_plt, aes(x = cell_type, y = CompositionCluster_CC, fill = scaled_pat_median)) +
  geom_tile(color = "white", linewidth = 0.2) +
  geom_text(aes(label = significant), color = "black", size = 5, fontface = "bold") +  # Increased size for prominence
  theme_minimal(base_size = 14, base_family = "Arial") +  # Increased base font size
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12, face = "bold"),  # Increased size and made bold
    axis.text.y = element_text(size = 12, face = "bold"),  # Increased size and made bold
    axis.title = element_text(size = 14, face = "bold"),  # Increased size
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_markdown(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"),
    axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),  # Add more space above x-axis title
    axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0))   # Add more space to the right of y-axis title
  ) +
  scale_fill_gradient2(
    low = "#3B4992",
    mid = "white",
    high = "#EE0000",
    midpoint = 0,
    name = "Scaled Median",
    limits = c(-max_abs_value, max_abs_value),
    breaks = pretty_breaks(n = 5)
  ) +
  labs(
    x = "Cell Type",
    y = "Composition Cluster",
    title = "Heatmap of Scaled Median by<br>Cell Type and Composition Cluster",
    subtitle = "Significance indicated by asterisks"
  ) +
  coord_equal()

# Print the plot


# Function to safely save the plot
safe_save <- function(filename, plot, width, height, ...) {
  tryCatch({
    ggsave(filename, plot, width = width, height = height, ...)
    cat("Successfully saved:", filename, "\n")
  }, error = function(e) {
    cat("Failed to save", filename, "- Error:", e$message, "\n")
  })
}

# Attempt to save using different methods
safe_save("Meenakshi_heatmap_1.pdf", ggplot_object, width = 12, height = 10, units = "in", dpi = 300, device = cairo_pdf)
safe_save("Meenakshi_heatmap_2.png", ggplot_object, width = 12, height = 10, units = "in", dpi = 300)
safe_save("Meenakshi_heatmap_3.eps", ggplot_object, width = 12, height = 10, units = "in", dpi = 300, device = cairo_ps)

# Check which files were created
file_list <- list.files(pattern = "Meenakshi_heatmap_\\d\\.(pdf|png|eps)")
cat("Created files:", paste(file_list, collapse = ", "), "\n")


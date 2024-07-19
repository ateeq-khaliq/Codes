# All in one plot

# renaming the spelling mistake 
# Load the dplyr package
library(dplyr)

# Assuming your tibble is named result_long
result_long <- result_long %>%
  mutate(treatment_group = case_when(
    treatment_group == "Folferi" ~ "Folfiri",
    treatment_group == "Progio_folferi" ~ "Proagio+Folfiri",
    TRUE ~ treatment_group
  ))

# View the corrected tibble
print(result_long)

# Load required packages
library(ggplot2)
library(ggsignif)
library(tidyr)
library(dplyr)
library(ggpubr)


# Generate all possible pairwise combinations of treatment groups
treatment_combinations <- combn(treatment_groups, 2, simplify = FALSE)

# Define custom colors
custom_colors <- c("Control" = "#66c2a5", "Folfiri" = "#fc8d62", "Proagio" = "#8da0cb", "Proagio+Folfiri" = "#e78ac3")

# Create the box plot with significance levels
p <- ggplot(result_long, aes(x = cell_type, y = percent, fill = treatment_group)) +
  geom_boxplot(color = "black", outlier.color = "black", outlier.shape = 21, outlier.size = 2) +
  stat_summary(fun = median, geom = "line", position = position_dodge(0.75), color = "black", size = 0.5) +
  scale_fill_manual(values = custom_colors, name = "Treatment Group") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white"),
    axis.line = element_line(size = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 16, hjust = 0.5, color = "black", face = "bold")
  ) +
  labs(x = "Cell Type", y = "Percentage", title = "Comparison of Treatment Groups by Broad Cell Types in Mouse Single cell Data")

# Add significance bars and asterisks for all pairwise combinations of treatment groups within each cell type
for (cell_type in unique(result_long$cell_type)) {
  for (i in seq_along(treatment_combinations)) {
    p <- p + stat_compare_means(
      comparisons = list(treatment_combinations[[i]]),
      data = subset(result_long, cell_type == cell_type),
      aes(group = treatment_group),
      method = "wilcox.test",
      label = "p.signif",
      label.y = max(subset(result_long, cell_type == cell_type)$percent) + 0.05,
      size = 3,
      tip.length = 0.01,
      step.increase = 0.05
    ) +
      stat_compare_means(
        comparisons = list(treatment_combinations[[i]]),
        data = subset(result_long, cell_type == cell_type),
        aes(group = treatment_group),
        method = "wilcox.test",
        label = "p.format",
        label.y = max(subset(result_long, cell_type == cell_type)$percent) + 0.1,
        size = 3,
        tip.length = 0,
        step.increase = 0,
        hide.ns = TRUE
      )
  }
}

# Print the plot

ggsave("comparison_treatment_groups_all.pdf", plot = p, width = 15, height = 8)

#######

# individual
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)

# Assuming your data is stored in a data frame called 'result'
# Reshape the data frame from wide to long format
result_long <- result %>%
  pivot_longer(cols = -c(samples, treatment_group),
               names_to = "cell_type",
               values_to = "percent")

# Define a vibrant color palette
palette <- c("#E74C3C", "#3498DB", "#2ECC71", "#9B59B6")

# Create the plot with significance comparisons
p <- ggplot(result_long, aes(x = treatment_group, y = percent, fill = treatment_group)) +
  geom_boxplot(color = "black", outlier.color = "black") +
  stat_summary(fun = median, geom = "line", position = position_dodge(0.75), color = "black", size = 0.5) +
  facet_wrap(~cell_type, scales = "free") +
  scale_fill_manual(values = palette,
                    name = "Treatment Group") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "gray95"),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 0.6, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, color = "black", face = "bold"),
    axis.text.y = element_text(size = 14, color = "black", face = "bold"),
    axis.title = element_text(size = 16, color = "black", face = "bold"),
    strip.text = element_text(size = 16, color = "black", face = "bold"),
    legend.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 20, hjust = 0.5, color = "black", face = "bold"),
    plot.subtitle = element_text(size = 16, hjust = 0.5, color = "gray30"),
    plot.caption = element_text(size = 12, hjust = 1, color = "gray30", face = "italic"),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  ) +
  labs(x = "Treatment Group",
       y = "Percentage",
       title = "Cell Type Distribution by Treatment Group",
       subtitle = "Comparison of cell type proportions across different treatments",
       caption = "Data source: scRNA-seq Mouse Data (CT26) - Masood Lab") +
  stat_compare_means(
    aes(group = treatment_group),
    method = "wilcox.test",
    label = "p.signif",
    size = 6,  # Increase the size of the significance labels
    vjust = 0.5,
    step.increase = 0.05,
    comparisons = list(c("Control", "Folferi"),
                       c("Control", "Proagio"),
                       c("Control", "Progio_folferi"),
                       c("Folferi", "Proagio"),
                       c("Folferi", "Progio_folferi"),
                       c("Proagio", "Progio_folferi")),
    label.x.npc = "center",  # Center align the significance labels
    label.y.npc = "top",  # Position the significance labels at the top
    label.fontface = "bold",  # Make the significance labels bold
    symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
    hide.ns = FALSE  # Show non-significant comparisons as well
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))  # Add space at the top for significance labels

# Save the plot as a PDF
ggsave("boxplot_results_Individual.pdf", plot = p, width = 50, height = 60, units = "cm")

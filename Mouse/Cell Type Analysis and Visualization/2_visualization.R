# Part 2: Visualization
# Save this as "2_visualization.R"

# Load required packages
library(ggplot2)
library(dplyr)
library(viridis)
library(patchwork)

create_boxplot_with_comparisons <- function(data, comprehensive_stats) {
  # Calculate max y values for each group for bracket positioning
  max_y_values <- data %>%
    group_by(assigned_MP) %>%
    summarise(max_y = max(proportion))

  # Create boxplot
  boxplot <- ggplot(data, aes(x = treatment_group, y = proportion)) +
    # Add boxplots
    geom_boxplot(
      aes(fill = treatment_group),
      outlier.shape = NA,
      width = 0.7
    ) +
    # Add individual points
    geom_point(
      position = position_jitter(width = 0.15),
      size = 1,
      alpha = 0.5,
      color = "black"
    ) +
    # Add significance brackets for significant comparisons
    geom_segment(
      data = comprehensive_stats %>% 
        filter(min_p_adj < 0.05) %>%
        mutate(
          x1 = as.numeric(factor(group1)),
          x2 = as.numeric(factor(group2)),
          y_pos = mapply(function(mp, i) {
            max_val <- max_y_values$max_y[max_y_values$assigned_MP == mp]
            return(max_val * (1.1 + (i-1)*0.1))
          }, assigned_MP, seq_along(assigned_MP))
        ),
      aes(
        x = x1,
        xend = x1,
        y = y_pos - y_pos*0.02,
        yend = y_pos
      ),
      size = 0.5
    ) +
    geom_segment(
      data = comprehensive_stats %>% 
        filter(min_p_adj < 0.05) %>%
        mutate(
          x1 = as.numeric(factor(group1)),
          x2 = as.numeric(factor(group2)),
          y_pos = mapply(function(mp, i) {
            max_val <- max_y_values$max_y[max_y_values$assigned_MP == mp]
            return(max_val * (1.1 + (i-1)*0.1))
          }, assigned_MP, seq_along(assigned_MP))
        ),
      aes(
        x = x1,
        xend = x2,
        y = y_pos,
        yend = y_pos
      ),
      size = 0.5
    ) +
    geom_segment(
      data = comprehensive_stats %>% 
        filter(min_p_adj < 0.05) %>%
        mutate(
          x1 = as.numeric(factor(group1)),
          x2 = as.numeric(factor(group2)),
          y_pos = mapply(function(mp, i) {
            max_val <- max_y_values$max_y[max_y_values$assigned_MP == mp]
            return(max_val * (1.1 + (i-1)*0.1))
          }, assigned_MP, seq_along(assigned_MP))
        ),
      aes(
        x = x2,
        xend = x2,
        y = y_pos - y_pos*0.02,
        yend = y_pos
      ),
      size = 0.5
    ) +
    # Add significance stars
    geom_text(
      data = comprehensive_stats %>% 
        filter(min_p_adj < 0.05) %>%
        mutate(
          x1 = as.numeric(factor(group1)),
          x2 = as.numeric(factor(group2)),
          y_pos = mapply(function(mp, i) {
            max_val <- max_y_values$max_y[max_y_values$assigned_MP == mp]
            return(max_val * (1.1 + (i-1)*0.1))
          }, assigned_MP, seq_along(assigned_MP))
        ),
      aes(
        x = (x1 + x2)/2,
        y = y_pos,
        label = significance
      ),
      vjust = -0.5,
      size = 4
    ) +
    # Use viridis colors
    scale_fill_viridis_d(option = "D", alpha = 0.7) +
    # Facet by cell type
    facet_wrap(~assigned_MP, scales = "free_y", ncol = 4) +
    # Theme adjustments
    theme_minimal() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "gray80"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "gray95", color = NA),
      strip.text = element_text(face = "bold", size = 12),
      legend.position = "none"
    ) +
    labs(
      x = "Treatment Group",
      y = "Proportion",
      title = "Cell Type Proportions Across Treatment Groups",
      subtitle = "Significance levels: *** p<0.001, ** p<0.01, * p<0.05"
    )
  
  return(boxplot)
}

create_statistical_plots <- function(comprehensive_stats, data) {
  # Create all plots (including previous ones)
  # ... [Previous plot code remains the same] ...
  
  # Add boxplot
  p4 <- create_boxplot_with_comparisons(data, comprehensive_stats)
  
  # Return all plots
  return(list(
    heatmap = p1,
    effect_size = p2,
    fold_change = p3,
    boxplot = p4
  ))
}

save_plots <- function(plots, prefix = "statistical") {
  # Save individual plots
  ggsave(paste0(prefix, "_heatmap.pdf"), plots$heatmap, width = 12, height = 8)
  ggsave(paste0(prefix, "_effect_size.pdf"), plots$effect_size, width = 12, height = 8)
  ggsave(paste0(prefix, "_fold_change.pdf"), plots$fold_change, width = 12, height = 8)
  ggsave(paste0(prefix, "_boxplot.pdf"), plots$boxplot, width = 14, height = 10)
  
  # Create and save combined plot (optional)
  combined_plot <- (plots$heatmap / plots$effect_size / plots$fold_change / plots$boxplot) +
    plot_layout(heights = c(1, 1, 1, 1.2)) +
    plot_annotation(
      title = "Comprehensive Analysis of Cell Type Differences",
      subtitle = "Significance levels: *** p<0.001, ** p<0.01, * p<0.05, ns = not significant",
      theme = theme(
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 10)
      )
    )
  
  ggsave(paste0(prefix, "_combined.pdf"), combined_plot, width = 15, height = 25)
}

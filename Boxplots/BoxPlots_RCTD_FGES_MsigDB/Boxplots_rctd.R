# Read the CSV file
data <- read.csv("/Users/akhaliq/Desktop/Gene_NMf/Exploratory_analysis/mean_results_reshaped2_fges_new.csv")
# Read the CSV file
data <- read.csv("/Users/akhaliq/Desktop/cell_trek_all/mean_results_reshaped2_rctd_new.csv")#RCTD
data <- read.csv("/Users/akhaliq/Desktop/Gene_NMf/Exploratory_analysis/mean_results_reshaped2_fges_new.csv")#FGES
data<-read.csv("/Users/akhaliq/Desktop/cell_trek_all/mean_results_reshaped_mgsigdb.csv")#MSIGDB

# Filter out unnecessary columns and only include CC1, CC2, CC3, CC4, and CC5
data_filtered <- subset(data, cc_ischia_10 %in% c("CC1", "CC2","CC3","CC5"))

# Melt the filtered data frame to long format
data_melted <- melt(data_filtered, id.vars = c("rowname", "cc_ischia_10"))

# Create box plots for each variable with significance indicators
p <- ggplot(data_melted, aes(x = cc_ischia_10, y = value, fill = cc_ischia_10)) +
  geom_boxplot(color = "black", outlier.color = "black") +
  stat_summary(fun = median, geom = "line", position = position_dodge(0.75), color = "black", size = 0.5) +
  facet_wrap(~variable, scales = "free", ncol = 4) +
  scale_fill_manual(values = c("CC1" = "#9467bd", 
                               "CC2" = "#ff7f0e",
                               "CC3" = "#2ca02c",
                               "CC4" = "#1f77b4",
                               "CC5" = "#d62728"),
                    name = "CC Ischia Type") +
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(size = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6, color = "black", face = "bold"),
    axis.text.y = element_text(size = 6, color = "black", face = "bold"),
    axis.title = element_text(size = 8, color = "black", face = "bold"),
    strip.text = element_text(size = 6, color = "black", face = "bold"),
    legend.position = "none",
    plot.title = element_text(size = 10, hjust = 0.5, color = "black", face = "bold"),
    panel.spacing = unit(0.5, "lines")
  ) +
  labs(x = "CC Ischia Type",
       y = "Value") +
  stat_compare_means(
    comparisons = list(c("CC1", "CC2"), 
                       c("CC1", "CC3"), 
                       c("CC1", "CC5"), 
                       c("CC2", "CC3"), 
                       c("CC2", "CC5"),
                       c("CC3", "CC5")),
    method = "t.test",
    label = "p.signif",
    size = 1.5,
    vjust = 0.5,
    step.increase = 0.05
  )

# Save the plot as a PDF
ggsave("boxplot_results_msigdb.pdf", plot = p, width = 30, height = 80, units = "cm")

# Create an empty data frame to store the results
result_df <- data.frame(
  CC = character(),
  cell_type = character(),
  statistic = numeric(),
  p.value = numeric(),
  method = character(),
  alternative = character(),
  p_corr = numeric(),
  significant = character(),
  stringsAsFactors = FALSE
)

# Loop through each combination of CC and cell type
for (cc in unique(CC1_groups1$CC)) {
  for (cell_type in names(CC1_groups1)[-1]) {
    # Perform Wilcoxon rank sum test for each combination
    test_result <- wilcox.test(CC1_groups1[CC1_groups1$CC == cc, cell_type], 
                               CC1_groups1[CC1_groups1$CC != cc, cell_type])
    
    # Append the results to the result_df data frame
    result_df <- rbind(result_df, data.frame(
      CC = cc,
      cell_type = cell_type,
      statistic = test_result$statistic,
      p.value = test_result$p.value,
      method = "Wilcoxon rank sum test (Mann-Whitney U test)",
      alternative = ifelse(test_result$alternative == "two.sided", "two-sided", "greater"),
      p_corr = test_result$p.value * (length(unique(CC1_groups1$CC)) * length(names(CC1_groups1)[-1])),
      significant = ifelse(test_result$p.value < 0.05, "*", ""),
      stringsAsFactors = FALSE
    ))
  }
}

# Print the result_df data frame
print(result_df)

### input Data

> CC1_groups1
     CC   Liver_Mets Lymph node         PDAC Normal Pancreas
1   CC1 0.0512126462 0.01342573 0.9353616284    0.0000000000
2  CC10 0.0001029972 0.00000000 0.0001029972    0.9997940056
3   CC2 0.3632778610 0.23322029 0.3993217131    0.0041801404
4   CC3 0.3919855378 0.35545345 0.2498493522    0.0027116601
5   CC4 0.2003978120 0.13832256 0.6610309962    0.0002486325
6   CC5 0.0870258991 0.03363201 0.8793420891    0.0000000000
7   CC6 1.0000000000 0.00000000 0.0000000000    0.0000000000
8   CC7 0.5350838155 0.08025516 0.3815457647    0.0031152648
9   CC8 1.0000000000 0.00000000 0.0000000000    0.0000000000
10  CC9 0.0122557942 0.89964810 0.0880961048    0.0000000000
> 

Its proportions#
# if you want to use it in this format 
melted_data <- reshape2::melt(CC1_groups1, id.vars = "CC", variable.name = "Group", value.name = "Value")




# Define the data for CC1 only
CC1_data <- CC1_groups1[CC1_groups1$CC == "CC1", c("PDAC", "Liver_Mets")]

# Perform Mann-Whitney U test (Wilcoxon rank sum test) for CC1
wilcox_result_CC1 <- wilcox.test(CC1_data$PDAC, CC1_data$Liver_Mets, y = NULL, alternative = "greater", paired = TRUE )

wilcox.test(CC1_data$PDAC, CC1_data$Liver_Mets, y = NULL,
            alternative = c("two.sided", "less", "greater"),
            mu = 0, paired = FALSE, exact = NULL, correct = TRUE,
            conf.int = FALSE, conf.level = 0.95)


# Output the result
print(wilcox_result_CC1)


CC1_groups1_tibble %>%
  group_by(group) %>%
  get_summary_stats(weight, type = "median_iqr")

ggboxplot( cc1_groups, x = "Histology", y = "Spots", ylab = "CC", xlab = "Histology", add = "jitter")

CC1_groups1 %>%
  group_by(Histology) %>%
  get_summary_stats(Spots, type = "median_iqr")


cc6_new <- subset(merged,CC=="CC2")

stat.test <- cc6_new %>% 
  wilcox_test(Spots ~ Histology) %>%
  add_significance()
stat.test

kruskal.test(Spots ~ Histology, data = cc1_groups)

cc2_new <- subset(merged,CC=="CC2")

ggboxplot(cc2_new, x = "Histology", y = "Spots", 
          color = "Histology", palette = c("#00AFBB", "#E7B800", "#FC4E07","blue"),
          ylab = "spots")+
  stat_compare_means(comparisons = list(c("PDAC", "Liver_Mets"), c("PDAC", "Lymph node"),c("PDAC", "Normal Pancreas"), c("Normal Pancreas", "Lymph node"), c("Liver_Mets", "Lymph node"),c("Liver_Mets", "Normal Pancreas")),
                     label = "p.signif", method = "wilcox.test")


##### Working Example for all 10 CCs
library(ggplot2)
library(ggpubr)

# List to store plots
plots <- list()

# Loop through each CC
for (cc in unique(merged$CC)) {
  cc_data <- subset(merged, CC == cc)
  
  # Generate ggplot for each CC
  p <- ggboxplot(cc_data, x = "Histology", y = "Spots", 
                 color = "Histology", palette = c("#00AFBB", "#E7B800", "#FC4E07","blue"),
                 ylab = "spots", 
                 add = "jitter",  # Add jitter to points for better visualization
                 outlier.shape = NA) +  # Do not plot outliers separately
    stat_compare_means(comparisons = list(c("PDAC", "Liver_Mets"), c("PDAC", "Lymph node"),c("PDAC", "Normal Pancreas"), c("Normal Pancreas", "Lymph node"), c("Liver_Mets", "Lymph node"),c("Liver_Mets", "Normal Pancreas")),
                       label = "p.signif", method = "wilcox.test", step.increase = 0.1, 
                       position = position_dodge(width = 0.5), # Dodge position to prevent overlap
                       color = "black", fill = "lightgray") +  # Set color and fill for significance bars
    theme(axis.text.x = element_blank(),  # Remove x-axis labels
          axis.ticks.x = element_blank()) +  # Remove x-axis ticks
    ggtitle(paste("CC", cc))  # Add CC as title
  
  # Add plot to list
  plots[[cc]] <- p
}

# Combine plots into a single plot with a common legend
combined_plot <- ggarrange(plotlist = plots, nrow = 2, ncol = 5, common.legend = TRUE, legend = "right", align = "hv")

# Save the combined plot to a PDF
ggsave("combined_plots.pdf", combined_plot, width = 14, height = 7)


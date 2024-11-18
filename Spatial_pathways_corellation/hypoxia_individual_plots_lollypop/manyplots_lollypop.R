# Load necessary libraries
library(corrplot)
library(ggplot2)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(viridis)

# Read and prepare the data
es <- readRDS("/Users/akhaliq/Desktop/asif/featureplots/enrichment_scores_mgsib.rds")
es <- normalize_weights(es)
colnames(es) <- sub("^HALLMARK_", "", colnames(es))
es_df <- as.data.frame(es)

# Calculate correlation matrix
cor_matrix <- cor(es_df)
cor_with_hypoxia <- cor_matrix["HYPOXIA", ]
cor_df <- data.frame(
  Pathway = names(cor_with_hypoxia),
  Correlation = cor_with_hypoxia
) %>%
  arrange(desc(abs(Correlation)))

# Create PDF for all visualizations
pdf("pathway_correlation_analysis.pdf", width = 12, height = 15)

# 1. Correlation Matrix Heatmap
corrplot(cor_matrix,
         method = "color",
         type = "upper",
         order = "hclust",
         tl.col = "black",
         tl.srt = 45,
         tl.cex = 0.7,
         addCoef.col = "black",
         number.cex = 0.5,
         col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(200))
title("Pathway Correlation Matrix")

# 2. Enhanced Lollipop Chart
cor_df_filtered <- cor_df %>%
  filter(Pathway != "HYPOXIA") %>%
  mutate(
    Correlation_Type = ifelse(Correlation > 0, "Positive", "Negative"),
    Significance = case_when(
      abs(Correlation) >= 0.7 ~ "Strong",
      abs(Correlation) >= 0.4 ~ "Moderate",
      TRUE ~ "Weak"
    ),
    Pathway = factor(Pathway, levels = Pathway[order(Correlation)])
  )

lollipop_plot <- ggplot(cor_df_filtered, aes(x = Pathway, y = Correlation)) +
  geom_segment(aes(xend = Pathway, yend = 0, color = Significance), 
              size = 1) +
  geom_point(aes(color = Significance), size = 3) +
  scale_color_manual(values = c("Strong" = "#B2182B", 
                               "Moderate" = "#FF7F00", 
                               "Weak" = "#4575B4")) +
  coord_flip() +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 10),
    axis.title.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = "Pathway Correlations with Hypoxia",
    y = "Correlation Coefficient",
    color = "Correlation\nStrength"
  )

print(lollipop_plot)

# 3. Top Correlations Bar Plot
top_correlations <- cor_df_filtered %>%
  arrange(desc(abs(Correlation))) %>%
  head(15)

bar_plot <- ggplot(top_correlations, 
                   aes(x = reorder(Pathway, abs(Correlation)), 
                       y = Correlation,
                       fill = Correlation)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient2(low = "#2166AC", 
                      mid = "white", 
                      high = "#B2182B", 
                      midpoint = 0) +
  coord_flip() +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10),
    legend.position = "right"
  ) +
  labs(
    title = "Top 15 Pathways Correlated with Hypoxia",
    x = "Pathway",
    y = "Correlation Coefficient",
    fill = "Correlation"
  )

print(bar_plot)

dev.off()

# Generate summary statistics
summary_stats <- cor_df_filtered %>%
  summarise(
    Mean_Correlation = mean(Correlation),
    Median_Correlation = median(Correlation),
    Std_Dev = sd(Correlation),
    Strong_Positive = sum(Correlation > 0.7),
    Moderate_Positive = sum(Correlation > 0.4 & Correlation <= 0.7),
    Strong_Negative = sum(Correlation < -0.7),
    Moderate_Negative = sum(Correlation < -0.4 & Correlation >= -0.7)
  )

# Print summary
cat("\nCorrelation Summary:\n")
print(summary_stats)

# Print top correlations
cat("\nTop 5 positive correlations with HYPOXIA:\n")
print(cor_df_filtered %>% 
        filter(Correlation > 0) %>% 
        arrange(desc(Correlation)) %>% 
        head(5))

cat("\nTop 5 negative correlations with HYPOXIA:\n")
print(cor_df_filtered %>% 
        filter(Correlation < 0) %>% 
        arrange(Correlation) %>% 
        head(5))

# Save results to CSV
write.csv(cor_df_filtered, "hypoxia_pathway_correlations.csv", row.names = FALSE)
write.csv(summary_stats, "hypoxia_correlation_summary.csv", row.names = FALSE)

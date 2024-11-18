# Load necessary libraries
library(corrplot)
library(ggplot2)
library(gridExtra)
library(cowplot)

# Read and prepare the data
es <- readRDS("/Users/akhaliq/Desktop/asif/featureplots/enrichment_scores_mgsib.rds")
es <- normalize_weights(es)
colnames(es) <- sub("^HALLMARK_", "", colnames(es))
es_df <- as.data.frame(es)

# Create function for scatter plots with correlation analysis
create_scatterplot <- function(data, x_pathway, y_pathway) {
  correlation <- round(cor(data[[x_pathway]], data[[y_pathway]]), 2)
  
  # Determine correlation direction
  correlation_type <- if(correlation > 0) "Positive" else "Negative"
  
  ggplot(data, aes_string(x = x_pathway, y = y_pathway)) +
    geom_bin2d(bins = 30) +
    scale_fill_gradientn(colors = c("navy", "blue", "cyan", "green", "yellow", "red"), 
                        name = "Count") +
    geom_smooth(method = "lm", se = FALSE, color = "white", size = 0.5) +
    theme_minimal() +
    labs(x = x_pathway, 
         y = y_pathway,
         title = paste0("r = ", correlation, "\n", correlation_type, " correlation")) +
    theme(axis.title = element_text(size = 6),
          axis.text = element_text(size = 5),
          plot.title = element_text(size = 7, face = "bold"),
          legend.position = "none",
          plot.margin = unit(c(1, 1, 1, 1), "mm")) +
    coord_fixed(ratio = 1)
}

# Get all pathways except HYPOXIA
other_pathways <- setdiff(colnames(es_df), "HYPOXIA")

# Create plots for HYPOXIA vs all other pathways
plots <- list()
correlations <- data.frame(
  Pathway = character(),
  Correlation = numeric(),
  Type = character(),
  stringsAsFactors = FALSE
)

for (pathway in other_pathways) {
  plots[[pathway]] <- create_scatterplot(es_df, "HYPOXIA", pathway)
  
  # Store correlation information
  cor_value <- round(cor(es_df[["HYPOXIA"]], es_df[[pathway]]), 2)
  correlations <- rbind(correlations, data.frame(
    Pathway = pathway,
    Correlation = cor_value,
    Type = ifelse(cor_value > 0, "Positive", "Negative"),
    stringsAsFactors = FALSE
  ))
}

# Sort correlations by absolute value
correlations <- correlations[order(abs(correlations$Correlation), decreasing = TRUE), ]

# Calculate layout
n_plots <- length(plots)
n_cols <- ceiling(sqrt(n_plots))
n_rows <- ceiling(n_plots / n_cols)

# Create PDF
pdf("Hypoxia_correlation_analysis.pdf", width = 3 * n_cols, height = 3 * n_rows)

# Arrange all plots in a grid
grid_arranged_plots <- do.call(gridExtra::grid.arrange, c(plots, ncol = n_cols))

# Print the arranged plots
print(grid_arranged_plots)

# Close the PDF device
dev.off()

# Write correlation summary to file
write.csv(correlations, "hypoxia_correlations.csv", row.names = FALSE)

# Print summary of strongest correlations
cat("\nTop 5 strongest positive correlations with HYPOXIA:\n")
print(head(correlations[correlations$Correlation > 0, ], 5))

cat("\nTop 5 strongest negative correlations with HYPOXIA:\n")
print(head(correlations[correlations$Correlation < 0, ], 5))

# ref = https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02663-5#availability-of-data-and-materials
# Load necessary libraries
library(corrplot)
library(spacexr)
library(Seurat)


es <- readRDS("/Users/akhaliq/Desktop/asif/featureplots/enrichment_scores_mgsib.rds")
es <- normalize_weights(es)

colnames(es) <- sub("^HALLMARK_", "", colnames(es))

# Calculate the correlation matrix
cor_matrix <- cor(es)

# Open a PDF device
pdf("pathway_correlations.pdf", width = 12, height = 12)  # Adjust width and height as needed

# Create the correlation heatmap
corrplot(cor_matrix, 
         method = "color",
         type = "upper",
         order = "hclust",
         tl.col = "black",
         tl.srt = 90,       # Rotate text labels vertically
         tl.cex = 0.6,      # Reduce text size
         diag = FALSE)

# Close the PDF device
dev.off()

###

library(ggplot2)
library(gridExtra)
library(cowplot)

# Assuming 'es' is your data frame with pathway activities
# If it's not already a data frame, convert it:


# Select the pathways you want to plot
#selected_pathways <- c("HYPOXIA", "E2F_TARGETS","EPITHELIAL_MESENCHYMAL_TRANSITION","OXIDATIVE_PHOSPHORYLATION","INFLAMMATORY_RESPONSE")  # Add more pathways as needed

# ONCOGENIC Pathways
selected_pathways <- c("G2M_CHECKPOINT","EPITHELIAL_MESENCHYMAL_TRANSITION","E2F_TARGETS","PI3K_AKT_MTOR_SIGNALING","MYC_TARGETS_V1")

# Signalling Pathways
selected_pathways <- c("HEDGEHOG_SIGNALING","NOTCH_SIGNALING","WNT_BETA_CATENIN_SIGNALING")

#Metabolism & Development
selected_pathways <- c("GLYCOLYSIS","ANGIOGENESIS","ADIPOGENESIS","OXIDATIVE_PHOSPHORYLATION")

#Immune
selected_pathways <- c("IL2_STAT5_SIGNALING","TNFA_SIGNALING_VIA_NFKB","INTERFERON_ALPHA_RESPONSE","INTERFERON_GAMMA_RESPONSE")


# Assuming 'es' is your data frame with pathway activities
es_df <- as.data.frame(es)

create_scatterplot <- function(data, x_pathway, y_pathway) {
  correlation <- round(cor(data[[x_pathway]], data[[y_pathway]]), 2)
  
  ggplot(data, aes_string(x = x_pathway, y = y_pathway)) +
    geom_bin2d(bins = 30) +
    scale_fill_gradientn(colors = c("navy", "blue", "cyan", "green", "yellow", "red"), name = "Count") +
    geom_smooth(method = "lm", se = FALSE, color = "white", size = 0.5) +
    theme_minimal() +
    labs(x = x_pathway, y = y_pathway,
         title = paste("r =", correlation)) +
    theme(axis.title = element_text(size = 6),
          axis.text = element_text(size = 5),
          plot.title = element_text(size = 7, face = "bold"),
          legend.position = "none",
          plot.margin = unit(c(1, 1, 1, 1), "mm")) +
    coord_fixed(ratio = 1)
}

# Generate all pairwise scatterplots
plots <- list()
for (i in 1:(length(selected_pathways) - 1)) {
  for (j in (i + 1):length(selected_pathways)) {
    plots[[length(plots) + 1]] <- create_scatterplot(es_df, selected_pathways[i], selected_pathways[j])
  }
}

# Calculate the number of rows and columns needed
n_plots <- length(plots)
n_cols <- ceiling(sqrt(n_plots))
n_rows <- ceiling(n_plots / n_cols)

# Create a single-page PDF
pdf("Immune_pathway_scatterplots.pdf", width = 2.5 * n_cols, height = 2.5 * n_rows)

# Arrange all plots in a grid
grid_arranged_plots <- do.call(gridExtra::grid.arrange, c(plots, ncol = n_cols))

# Print the arranged plots
print(grid_arranged_plots)

# Close the PDF device
dev.off()

cat("Oncogenic pathway scatterplots have been saved to 'Oncogenic_pathway_scatterplots.pdf'\n")

######

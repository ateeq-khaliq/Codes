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

selected_pathways <- c("E2F_TARGETS","HYPOXIA")

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
pdf("E2F_Hypoxia_pathway_scatterplots.pdf", width = 2.5 * n_cols, height = 2.5 * n_rows)

# Arrange all plots in a grid
grid_arranged_plots <- do.call(gridExtra::grid.arrange, c(plots, ncol = n_cols))

# Print the arranged plots
print(grid_arranged_plots)

# Close the PDF device
dev.off()

cat("Oncogenic pathway scatterplots have been saved to 'Oncogenic_pathway_scatterplots.pdf'\n")

######
# We can combine all these pathways into a single analysis to check the correlations between them. 
library(ggplot2)
library(gridExtra)
library(reshape2)
library(corrplot)

# Combine all pathways
all_pathways <- c(
  # ONCOGENIC Pathways
  "G2M_CHECKPOINT", "EPITHELIAL_MESENCHYMAL_TRANSITION", "E2F_TARGETS", "PI3K_AKT_MTOR_SIGNALING", "MYC_TARGETS_V1",
  # Signalling Pathways
  "HEDGEHOG_SIGNALING", "NOTCH_SIGNALING", "WNT_BETA_CATENIN_SIGNALING",
  # Metabolism & Development
  "GLYCOLYSIS", "ANGIOGENESIS", "ADIPOGENESIS", "OXIDATIVE_PHOSPHORYLATION",
  # Immune
  "IL2_STAT5_SIGNALING", "TNFA_SIGNALING_VIA_NFKB", "INTERFERON_ALPHA_RESPONSE", "INTERFERON_GAMMA_RESPONSE"
)

# Assuming 'es' is your data frame with pathway activities
es_df <- as.data.frame(es)

# Subset the data to include only the selected pathways
es_subset <- es_df[, all_pathways]

# Calculate the correlation matrix
cor_matrix <- cor(es_subset)

# Create a correlation heatmap
pdf("All_pathways_correlation_heatmap.pdf", width = 12, height = 10)
corrplot(cor_matrix, method = "color", type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, addCoef.col = "black", 
         number.cex = 0.7, tl.cex = 0.7)
dev.off()

# Function to create scatterplot
create_scatterplot <- function(data, x_pathway, y_pathway) {
  correlation <- round(cor(data[[x_pathway]], data[[y_pathway]]), 2)
  
  ggplot(data, aes_string(x = x_pathway, y = y_pathway)) +
    geom_point(alpha = 0.3, size = 1, shape = 16) +
    geom_smooth(method = "lm", se = FALSE, color = "red", size = 0.5) +
    theme_minimal() +
    labs(x = x_pathway, y = y_pathway,
         title = paste("r =", correlation)) +
    theme(axis.title = element_text(size = 6),
          axis.text = element_text(size = 5),
          plot.title = element_text(size = 7, face = "bold"),
          plot.margin = unit(c(1, 1, 1, 1), "mm")) +
    coord_fixed(ratio = 1)
}

# Generate all pairwise scatterplots
plots <- list()
for (i in 1:(length(all_pathways) - 1)) {
  for (j in (i + 1):length(all_pathways)) {
    plots[[length(plots) + 1]] <- create_scatterplot(es_subset, all_pathways[i], all_pathways[j])
  }
}

# Calculate the number of rows and columns needed
n_plots <- length(plots)
n_cols <- ceiling(sqrt(n_plots))
n_rows <- ceiling(n_plots / n_cols)

# Create a multi-page PDF for scatterplots
pdf("All_pathways_scatterplots.pdf", width = 20, height = 20)
for(i in seq(1, length(plots), 12)) {
  print(do.call(gridExtra::grid.arrange, c(plots[i:min(i+11, length(plots))], ncol = 4)))
}
dev.off()

cat("All pathway correlations have been analyzed and saved to PDF files.\n")

##
# the script to calculate correlations for all pathways present in your 'es' dataset.

library(writexl)
library(dplyr)

# Assuming 'es' is your data frame with all pathway activities
es_df <- as.data.frame(es)

# Calculate the correlation matrix for all pathways
cor_matrix <- cor(es_df)

# Convert the correlation matrix to a long format
cor_long <- as.data.frame(as.table(cor_matrix))
names(cor_long) <- c("Pathway1", "Pathway2", "Correlation")

# Remove self-correlations and duplicate pairs
cor_long <- cor_long %>%
  filter(as.character(Pathway1) < as.character(Pathway2)) %>%
  mutate(AbsCorrelation = abs(Correlation)) %>%
  arrange(desc(AbsCorrelation))

# Function to categorize pathways
categorize_pathway <- function(pathway) {
  if (pathway %in% c("G2M_CHECKPOINT", "EPITHELIAL_MESENCHYMAL_TRANSITION", "E2F_TARGETS", "PI3K_AKT_MTOR_SIGNALING", "MYC_TARGETS_V1")) {
    return("ONCOGENIC")
  } else if (pathway %in% c("HEDGEHOG_SIGNALING", "NOTCH_SIGNALING", "WNT_BETA_CATENIN_SIGNALING")) {
    return("SIGNALLING")
  } else if (pathway %in% c("GLYCOLYSIS", "ANGIOGENESIS", "ADIPOGENESIS", "OXIDATIVE_PHOSPHORYLATION")) {
    return("METABOLISM & DEVELOPMENT")
  } else if (pathway %in% c("IL2_STAT5_SIGNALING", "TNFA_SIGNALING_VIA_NFKB", "INTERFERON_ALPHA_RESPONSE", "INTERFERON_GAMMA_RESPONSE")) {
    return("IMMUNE")
  } else {
    return("OTHER")
  }
}

# Function to categorize correlation strength
categorize_correlation <- function(cor_value) {
  if (cor_value >= 0.7) {
    return("Strong Positive")
  } else if (cor_value >= 0.3 && cor_value < 0.7) {
    return("Moderate Positive")
  } else if (cor_value > 0 && cor_value < 0.3) {
    return("Weak Positive")
  } else if (cor_value == 0) {
    return("No Correlation")
  } else if (cor_value > -0.3 && cor_value < 0) {
    return("Weak Negative")
  } else if (cor_value >= -0.7 && cor_value <= -0.3) {
    return("Moderate Negative")
  } else {
    return("Strong Negative")
  }
}

# Add categories for both pathways and correlation strength
cor_long <- cor_long %>%
  mutate(
    Category1 = sapply(as.character(Pathway1), categorize_pathway),
    Category2 = sapply(as.character(Pathway2), categorize_pathway),
    CorrelationStrength = sapply(Correlation, categorize_correlation)
  )

# Write to Excel
write_xlsx(cor_long, "all_pathway_correlations_with_strength.xlsx")

cat("Correlation data for all pathways has been saved to 'all_pathway_correlations_with_strength.xlsx'\n")

# Print summary statistics
cat("\nSummary of correlations:\n")
print(summary(cor_long$Correlation))

cat("\nDistribution of correlation strengths:\n")
print(table(cor_long$CorrelationStrength))

cat("\nTop 10 strongest positive correlations:\n")
print(head(cor_long[cor_long$Correlation > 0, c("Pathway1", "Pathway2", "Correlation", "CorrelationStrength")], 10))

cat("\nTop 10 strongest negative correlations:\n")
print(head(cor_long[cor_long$Correlation < 0, c("Pathway1", "Pathway2", "Correlation", "CorrelationStrength")], 10))

###


Top 10 strongest positive correlations:
> print(head(cor_long[cor_long$Correlation > 0, c("Pathway1", "Pathway2", "Correlation", "CorrelationStrength")], 10))
                    Pathway1                  Pathway2 Correlation
1  INTERFERON_ALPHA_RESPONSE INTERFERON_GAMMA_RESPONSE   0.7795299
2        IL2_STAT5_SIGNALING            UV_RESPONSE_DN   0.7409674
4    ESTROGEN_RESPONSE_EARLY    ESTROGEN_RESPONSE_LATE   0.7187136
6                 DNA_REPAIR OXIDATIVE_PHOSPHORYLATION   0.7118869
8        IL2_STAT5_SIGNALING     INFLAMMATORY_RESPONSE   0.7037582
9        ALLOGRAFT_REJECTION     INFLAMMATORY_RESPONSE   0.7024485
11               E2F_TARGETS            G2M_CHECKPOINT   0.6972485
13            MYC_TARGETS_V1 OXIDATIVE_PHOSPHORYLATION   0.6865774
15     INFLAMMATORY_RESPONSE   TNFA_SIGNALING_VIA_NFKB   0.6843531
17       ALLOGRAFT_REJECTION INTERFERON_GAMMA_RESPONSE   0.6727342
   CorrelationStrength
1      Strong Positive
2      Strong Positive
4      Strong Positive
6      Strong Positive
8      Strong Positive
9      Strong Positive
11   Moderate Positive
13   Moderate Positive
15   Moderate Positive
17   Moderate Positive
> 
> cat("\nTop 10 strongest negative correlations:\n")

Top 10 strongest negative correlations:
> print(head(cor_long[cor_long$Correlation < 0, c("Pathway1", "Pathway2", "Correlation", "CorrelationStrength")], 10))
                            Pathway1                  Pathway2 Correlation
3          OXIDATIVE_PHOSPHORYLATION            UV_RESPONSE_DN  -0.7198642
5              INFLAMMATORY_RESPONSE OXIDATIVE_PHOSPHORYLATION  -0.7160256
7  EPITHELIAL_MESENCHYMAL_TRANSITION OXIDATIVE_PHOSPHORYLATION  -0.7077599
10               IL2_STAT5_SIGNALING OXIDATIVE_PHOSPHORYLATION  -0.7015901
12                        DNA_REPAIR            UV_RESPONSE_DN  -0.6904678
14   REACTIVE_OXYGEN_SPECIES_PATHWAY            UV_RESPONSE_DN  -0.6851453
16                 KRAS_SIGNALING_UP OXIDATIVE_PHOSPHORYLATION  -0.6730319
18             INFLAMMATORY_RESPONSE         PROTEIN_SECRETION  -0.6673890
19 EPITHELIAL_MESENCHYMAL_TRANSITION            MYC_TARGETS_V1  -0.6664939
20                        DNA_REPAIR       IL2_STAT5_SIGNALING  -0.6639532
   CorrelationStrength
3      Strong Negative
5      Strong Negative
7      Strong Negative
10     Strong Negative
12   Moderate Negative
14   Moderate Negative
16   Moderate Negative
18   Moderate Negative
19   Moderate Negative
20   Moderate Negative
> 



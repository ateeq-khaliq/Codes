# Load required packages
library(tidyverse)
library(dplyr)

assay_matrix <- pdac[["rctd_tier1"]]@data
norm_weights <- as.data.frame(t(assay_matrix))

# Read in the aggregated data
#agg_data <- read.csv("/Users/akhaliq/Desktop/spatial_analysis/pdac_new/ischia/k10/aggregated_props_patients1.csv", header = TRUE)
df <- merge(norm_weights,pdac$cc_ischia_10, by=0)
colnames(df)[colnames(df)=='y'] <- "Compositional_clusters"

aggregated_df <- df %>%
  group_by(Compositional_clusters) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE))
aggregated_df = as.data.frame(aggregated_df)

# Calculate row-wise percentages
#agg_data_percent <- apply(aggregated_df, 1, function(x) 100 * x / sum(x))

# Exclude the "Compositional_clusters" column before calculating percentages
numeric_data <- aggregated_df %>%
  select(-Compositional_clusters)

# Calculate row-wise percentages
agg_data_percent <- apply(numeric_data, 1, function(x) 100 * x / sum(x))
colnames(agg_data_percent) <- aggregated_df$Compositional_clusters
head(agg_data_percent)


# Convert percentages back to a data frame
agg_data_percent_df <- as.data.frame(t(agg_data_percent))
agg_data_percent_df$Compositional_clusters <- row.names(agg_data_percent_df)

# Color mapping

color_mapping <- c(
  "B cells" = "#800080",                   # Purple color.
  "C1Q-TAM" = "#ff7f0e",                   # A shade of orange.
  "CD4+ cells" = "#ff00ff",                # Magenta color.
  "CD8-NK cells" = "#d62728",              # A shade of red.
  "DCs" = "#000000",                       # Black color.
  "Endothelial cells" = "#d4f739",         # A shade of yellow-green.
  "FCN1-TAM" = "#00ff00",                  # Pure green color.
  "Hepatocytes" = "#09bfe8",               # A shade of blue.
  "iCAF" = "#0000ff",                      # Pure blue color.
  "myCAF" = "#2ca02c",                     # A shade of green.
  "Normal Epithelial cells" = "#AAB7B8",   # A gray color.
  "Proliferative T cells" = "#8b0000",     # Dark red color.
  "PVL" = "#32cd32",                       # Lime green color.
  "SPP1-TAM" = "#FFD700",                  # Gold color.
  "Tumor Epithelial cells" = "#ff0000"     # Red color.
)




# Melt the data frame to long format
melted_data <- agg_data_percent_df %>%
  gather(cell_type, percentage, -Compositional_clusters)

# Create the plot using ggplot
plot <- ggplot(melted_data, aes(x = Compositional_clusters, y = percentage, fill = cell_type)) +
  geom_bar(stat = "identity", position = "fill", color = "white") +
  labs(title = "Distribution of Cell Types by CC Type",
       x = "CC Type",
       y = "Percentage",
       fill = "Cell Type") +
  scale_fill_manual(values = color_mapping) +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.y = element_line(color = "gray90"),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"),
        legend.box.background = element_rect(color = "black"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 16, hjust = 0.5, face = "bold")
  )

# Save the plot as a PDF file
ggsave("aggregated_proportions_K10_tempus.pdf", plot, width = 10, height = 6)

## for percentaged to be plotted on the bar plots

# Create the plot using ggplot
plot <- ggplot(melted_data, aes(x = Compositional_clusters, y = percentage, fill = cell_type)) +
  geom_bar(stat = "identity", position = "fill", color = "white") +
  geom_text(aes(label = scales::percent(percentage / 100)),
            position = position_fill(vjust = 0.5),
            size = 2,
            color = "black",
            hjust = -0.2,
            vjust = 0.2) +
  labs(title = "Distribution of Cell Types by CC Type",
       x = "CC Type",
       y = "Percentage",
       fill = "Cell Type") +
  scale_fill_manual(values = color_mapping) +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
        panel.grid.major.y = element_line(color = "gray90"),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"),
        legend.box.background = element_rect(color = "black"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 16, hjust = 0.5, face = "bold")
  )

# Save the plot as a PDF file
ggsave("aggregated_proportions_K10_percentplot.pdf", plot, width = 10, height = 6)



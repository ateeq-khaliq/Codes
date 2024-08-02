### Run this python script to get aggregates on the proportions python /Users/akhaliq/Desktop/spatial_analysis/rctd/codes/aggregating_proportions.py

# Load required packages
library(tidyverse)

# Load required packages
library(tidyverse)

# Read in the aggregated data
agg_data <- read.csv("/Users/akhaliq/Desktop/spatial_analysis/rctd/aggregated_props.csv", header = TRUE)

# Assign orig.ident as row names
row.names(agg_data) <- agg_data$orig.ident
agg_data$orig.ident <- NULL

# Calculate row-wise percentages
agg_data_percent <- apply(agg_data, 1, function(x) 100*x/sum(x))

# Convert percentages back to a data frame
agg_data_percent_df <- as.data.frame(t(agg_data_percent))
agg_data_percent_df$group <- row.names(agg_data_percent_df)
agg_data_percent_df$orig.ident <- row.names(agg_data_percent_df)
agg_data_percent_df <- gather(agg_data_percent_df, key = "cell_type", value = "percentage", -group, -orig.ident)

my_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", 
               "#FFFF33", "#A65628", "#F781BF", "#999999", "#66CCFF")

# Use the defined color palette in the plot
ggplot(agg_data_percent_df, aes(x = orig.ident, y = percentage, fill = cell_type)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "Percentage Stacked Bar Plot",
       x = "Sample ID",
       y = "Percentage",
       fill = "Cell Type") +
  theme_minimal() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major.x = element_line(color = "gray80"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

#####

# To plot the percentage on the barplots use the following.

# Define color palette
my_palette <- c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF", "#00B0F6", "#FFC400", "#FF7F00", "#A1CAF1", "#4DAF4A")

# Plot percentage stacked bar chart
ggplot(agg_data_percent_df, aes(x = orig.ident, y = percentage, fill = cell_type)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(aes(label = paste0(round(percentage), "%")), position = position_fill(vjust = 0.5)) +
  scale_fill_manual(values = my_colors) +
  labs(title = "Percentage Stacked Bar Plot",
       x = "Sample ID",
       y = "Percentage",
       fill = "Cell Type") +
  theme_minimal() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1))


#####






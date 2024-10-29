# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

# Extract MP1 to MP7 scores and treated status from Seurat metadata
mp_data <- tc@meta.data %>%
  select(MP1, MP2, MP3, MP4, MP5, MP6, MP7, treated)

# Calculate z-score for each MP and add to the data
mp_data_z <- as.data.frame(scale(mp_data[, 1:7]))
mp_data_z$treated <- mp_data$treated

# Prepare data in long format for plotting
mp_long <- mp_data_z %>%
  pivot_longer(cols = starts_with("MP"), 
               names_to = "Metaprogram", values_to = "Zscore")

# Calculate p-values from t-tests
ttest_results <- mp_data %>%
  summarise(across(starts_with("MP"), 
                  ~t.test(.x ~ treated)$p.value))

# Prepare p-values with simplified formatting
p_values <- ttest_results %>%
  pivot_longer(everything(), 
              names_to = "Metaprogram", 
              values_to = "p_value") %>%
  mutate(
    p_label = case_when(
      p_value < 0.001 ~ "P < 0.001",
      TRUE ~ sprintf("P = %.2f", p_value)
    )
  )

# Define colors
box_colors <- c("#E69F00", "#2C7FB8")  # Orange and Blue

# Custom theme for publication-style plot
publication_theme <- theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_blank(),
    plot.title = element_text(size = 14, face = "bold", hjust = 0),
    legend.position = "top",
    legend.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  )

# Create the plot
pdf("MP_boxplots_publication.pdf", width = 12, height = 6)

ggplot(mp_long, aes(x = treated, y = Zscore, fill = treated)) +
  # Add boxplot
  geom_boxplot(outlier.shape = 16, outlier.size = 1, alpha = 0.7, width = 0.5,
               color = "black", size = 0.5) +
  # Customize appearance
  scale_fill_manual(values = box_colors,
                   labels = c("No", "Yes")) +
  facet_wrap(~Metaprogram, scales = "free_y", ncol = 4) +
  # Add p-values with parentheses
  geom_text(
    data = p_values,
    aes(x = 1.5, y = max(mp_long$Zscore) + 0.5, 
        label = p_label),
    size = 3, color = "black", inherit.aes = FALSE
  ) +
  # Labels with corrected y-axis
  labs(x = "Treatment status",
       y = "Z-score",
       title = "Metaprogram Expression Analysis") +
  # Apply custom theme
  publication_theme +
  # Set y-axis limits to ensure consistent spacing for p-values
  coord_cartesian(clip = "off")

dev.off()

# Print a message indicating PDF has been saved
print("PDF 'MP_boxplots_publication.pdf' has been saved in your working directory.")

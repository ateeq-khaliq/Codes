#Donut Plots for individual samples

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(gridExtra)
library(grid)

# Your data (assuming meta is already loaded in your environment)

# Define a more vibrant color palette from RColorBrewer
colors <- brewer.pal(8, "Dark2")

# Prepare the data for plotting
meta_summary <- meta %>%
  count(orig.ident, new_celltype) %>%
  group_by(orig.ident) %>%
  mutate(percent = n / sum(n)) %>%
  ungroup()

# Create a function to generate individual donut plots
create_donut_plot <- function(data, sample_name, colors) {
  data_filtered <- data %>% filter(orig.ident == sample_name)
  
  ggplot(data_filtered, aes(x = 2, y = percent, fill = new_celltype)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar(theta = "y") +
    xlim(0.5, 2.5) +  # Adjust to make the donut hole larger
    scale_fill_manual(values = colors) +
    theme_void() +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 8),
          plot.title = element_text(hjust = 0.5, size = 10, face = "bold")) +
    labs(title = sample_name)
}

# List of unique samples
samples <- unique(meta_summary$orig.ident)

# Generate a list of individual plots
plots <- lapply(samples, function(sample) {
  create_donut_plot(meta_summary, sample, colors)
})

# Save the plots to a PDF
pdf("Donut_Plots_Samples_by_Treatment.pdf", width = 20, height = 20)  # Open a PDF device with a larger size

# Arrange and print the plots
grid.arrange(
  grobs = plots,
  ncol = 4,
  top = textGrob("Distribution of Cell Types in Different Samples by Treatment Group", gp = gpar(fontsize = 16, fontface = "bold"))
)

dev.off()  # Close the PDF device
####


#Donut Plots for individual samples and percentages on the plot 

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(gridExtra)
library(grid)

# Your data (assuming meta is already loaded in your environment)

# Define a more vibrant color palette from RColorBrewer
colors <- brewer.pal(8, "Dark2")

# Prepare the data for plotting
meta_summary <- meta %>%
  count(orig.ident, new_celltype) %>%
  group_by(orig.ident) %>%
  mutate(percent = n / sum(n)) %>%
  ungroup()

# Create a function to generate individual donut plots with percentage labels
create_donut_plot <- function(data, sample_name, colors) {
  data_filtered <- data %>% filter(orig.ident == sample_name)
  
  ggplot(data_filtered, aes(x = 2, y = percent, fill = new_celltype)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar(theta = "y") +
    xlim(0.5, 2.5) +  # Adjust to make the donut hole larger
    scale_fill_manual(values = colors) +
    theme_void() +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 8),
          plot.title = element_text(hjust = 0.5, size = 10, face = "bold")) +
    labs(title = sample_name) +
    geom_text(aes(label = paste0(round(percent * 100, 1), "%")), 
              position = position_stack(vjust = 0.5), size = 3)
}

# List of unique samples
samples <- unique(meta_summary$orig.ident)

# Generate a list of individual plots
plots <- lapply(samples, function(sample) {
  create_donut_plot(meta_summary, sample, colors)
})

# Save the plots to a PDF
pdf("Donut_Plots_Samples_by_Treatment.pdf", width = 20, height = 20)  # Open a PDF device with a larger size

# Arrange and print the plots
grid.arrange(
  grobs = plots,
  ncol = 4,
  top = textGrob("Distribution of Cell Types in Different Samples by Treatment Group", gp = gpar(fontsize = 16, fontface = "bold"))
)

dev.off()  # Close the PDF device
#####

####
#Donut Plots Treatment Group Wise.

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(wesanderson)

# Your data (assuming meta is already loaded in your environment)

# Define a color palette from wesanderson
colors <- wes_palette("Zissou1", 8, type = "continuous")

# Create a function to generate donut plots for each treatment group
create_donut_plot <- function(data, treatment_group_name, colors) {
  data_filtered <- data %>%
    filter(treatment_group == treatment_group_name) %>%
    count(new_celltype) %>%
    mutate(percent = n / sum(n))
  
  ggplot(data_filtered, aes(x = 2, y = percent, fill = new_celltype)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar(theta = "y") +
    xlim(0.5, 2.5) +  # Adjust to make the donut hole larger
    scale_fill_manual(values = colors) +
    theme_void() +  # Remove unnecessary axis and background
    theme(legend.position = "right",
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.background = element_rect(fill = "lightblue", color = NA),  # Background color for the plot
          panel.background = element_rect(fill = "lightblue", color = NA)) +  # Background color for the panel
    labs(title = treatment_group_name)
}

# List of treatment groups
treatment_groups <- unique(meta$treatment_group)

# Generate donut plots for each treatment group
plots <- lapply(treatment_groups, function(group) {
  create_donut_plot(meta, group, colors)
})

# Save plots to a PDF
pdf("Donut_Plots.pdf", width = 14, height = 10)  # Open a PDF device
grid.arrange(
  grobs = plots,
  ncol = 2,
  top = textGrob("Distribution of Cell Types in Different Treatment Groups", gp = gpar(fontsize = 16, fontface = "bold"))
)
dev.off()  # Close the PDF device


######

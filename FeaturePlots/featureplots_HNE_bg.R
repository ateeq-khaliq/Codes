# Load required libraries
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(viridis)  # For an alternative color palette option

###
#Data Modification

library(Seurat)
library(spacexr)
library(stringr)

pdac <- readRDS("/Users/akhaliq/Desktop/asif/04-pdac_nac_subset_CC10.rds")
msigdb <- readRDS("/Users/akhaliq/Desktop/asif/featureplots/enrichment_scores_mgsib.rds")
norm_weights <- normalize_weights(msigdb)

# Remove 'HALLMARK_' prefix from column names
colnames(norm_weights) <- str_remove(colnames(norm_weights), "^HALLMARK_")

# Display the first few rows of the modified data frame
head(norm_weights)

pdac@assays[["msigdb"]] <- CreateAssayObject(data = t(as.matrix(norm_weights)))

# Seems to be a bug in SeuratData package that the key is not set and any
# plotting function etc. will throw an error.
if (length(pdac@assays$msigdb@key) == 0) {
    pdac@assays$msigdb@key = "msigdb_"
}

DefaultAssay(pdac) <- "msigdb"
###

# Load required libraries
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(scales)  # For gradient color functions

# Define your list of image names and features
image_names <- names(pdac@images)
features <- rownames(pdac@assays$msigdb)

# Create a directory to save the PDF files
dir.create("publication_quality_plots", showWarnings = FALSE)

# Define custom theme for the plots
publication_theme <- theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 8),
    legend.position = "bottom",
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.5, "cm"),
    legend.box = "horizontal",
    plot.margin = margin(5, 5, 5, 5),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank()
  )

# Define custom color palette
custom_palette <- function(n) {
  colors <- c("#00008B", "#4169E1", "#87CEEB", "#FFFFFF", "#FFA07A", "#FF4500", "#8B0000")
  colorRampPalette(colors)(n)
}

# Function to create a single plot
create_plot <- function(feature, image_name) {
  # Get the expression data for the feature
  expr_data <- FetchData(pdac, vars = feature)
  
  # Calculate percentiles for color scaling
  percentiles <- quantile(expr_data[[feature]], probs = c(0.01, 0.1, 0.5, 0.9, 0.99), na.rm = TRUE)
  
  SpatialFeaturePlot(pdac, 
                     features = feature, 
                     stroke = 0.1,
                     pt.size.factor = 1.8,
                     images = image_name,
                     image.alpha = 0.6) +
    scale_fill_gradientn(colors = custom_palette(100),
                         limits = c(percentiles[1], percentiles[5]),
                         breaks = percentiles,
                         labels = function(x) format(x, digits = 2, scientific = TRUE),
                         oob = scales::squish,
                         name = feature,
                         guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", barwidth = 7, barheight = 0.5)) +
    labs(title = image_name,
         subtitle = feature) +
    publication_theme +
    coord_fixed()
}

# Function to create a page of plots
create_page <- function(plots, feature) {
  title <- textGrob(paste("Spatial Distribution of", feature), 
                    gp = gpar(fontsize = 16, fontface = "bold"))
  
  # Ensure we have exactly 6 grobs (fill with nullGrob if needed)
  plot_grobs <- c(plots, replicate(max(0, 6 - length(plots)), nullGrob()))
  
  grid.arrange(
    title,
    arrangeGrob(grobs = plot_grobs, ncol = 3),
    nrow = 2,
    heights = c(0.1, 0.9)
  )
}

# Loop through features
for (feature in features) {
  feature_plots <- lapply(image_names, function(img) create_plot(feature, img))
  
  # Calculate the number of pages needed (6 plots per page)
  num_pages <- ceiling(length(feature_plots) / 6)
  
  # Create a PDF file for the current feature
  pdf_name <- paste0("publication_quality_plots/", feature, "_spatial_plot_numerical.pdf")
  pdf(pdf_name, width = 16, height = 12, onefile = TRUE)  # Landscape orientation
  
  for (page in 1:num_pages) {
    start_plot <- (page - 1) * 6 + 1
    end_plot <- min(page * 6, length(feature_plots))
    page_plots <- feature_plots[start_plot:end_plot]
    
    print(create_page(page_plots, feature))
  }
  
  dev.off()
}

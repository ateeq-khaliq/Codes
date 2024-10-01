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
###

# Define your list of image names and features
image_names <- names(pdac@images)
features <- rownames(pdac@assays$msigdb)

# Create a directory to save the PDF files
dir.create("publication_quality_plots", showWarnings = FALSE)

# Define custom theme for the plots
publication_theme <- theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 20)),
    plot.subtitle = element_text(hjust = 0.5, size = 12, margin = margin(b = 10)),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.size = unit(1, "cm"),
    plot.margin = margin(20, 20, 20, 20),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank()
  )

# Define custom color palette (option 1: custom blue to red)
custom_palette <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#F7F7F7", "#F4A582", "#D6604D", "#B2182B", "#67001F"))

# Alternative color palette (option 2: viridis)
# custom_palette <- viridis::viridis

# Function to create a single plot
create_plot <- function(feature, image_name) {
  SpatialFeaturePlot(pdac, 
                     features = feature, 
                     stroke = 0.2,
                     pt.size.factor = 1.5,
                     images = image_name) +
    scale_fill_gradientn(colors = custom_palette(100), 
                         limits = c(NA, NA),
                         oob = scales::squish,
                         name = feature,
                         guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
    labs(title = image_name,
         subtitle = feature) +
    publication_theme +
    coord_fixed()  # Ensure aspect ratio is 1:1
}

# Loop through features
for (feature in features) {
  feature_plots <- lapply(image_names, function(img) create_plot(feature, img))
  
  # Calculate the number of pages needed (2 plots per page for larger, more detailed plots)
  num_pages <- ceiling(length(feature_plots) / 2)
  
  # Create a PDF file for the current feature
  pdf_name <- paste0("publication_quality_plots/", feature, "_spatial_plot.pdf")
  pdf(pdf_name, width = 12, height = 18, onefile = TRUE)
  
  for (page in 1:num_pages) {
    start_plot <- (page - 1) * 2 + 1
    end_plot <- min(page * 2, length(feature_plots))
    page_plots <- feature_plots[start_plot:end_plot]
    
    # Create a vertical layout for better use of page space
    page_grid <- grid.arrange(
      grobs = c(
        list(textGrob(paste("Spatial Distribution of", feature), 
                      gp = gpar(fontsize = 20, fontface = "bold"))),
        page_plots
      ), 
      ncol = 1,
      heights = c(1, 10, 10)  # Adjust these values to control spacing
    )
    
    print(page_grid)
  }
  
  dev.off()
}

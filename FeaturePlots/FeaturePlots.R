### pre preprations


/Users/akhaliq/Desktop/spatial_analysis/pdac_new/ischia/new_k10/hallmark/allcc/ES_allcc.rds





library(Seurat)
library(decoupleR)
library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(ggplot2)
library(pheatmap)
library(stringr)

data=pdac

p <- readRDS("/Users/akhaliq/Desktop/spatial_analysis/pdac_new/ischia/new_k10/progeny/progeny_all.rds") # for progeny
p <- readRDS("/Users/akhaliq/Desktop/spatial_analysis/pdac_new/ischia/new_k10/progeny/progeny_minuscc.rds") # for progeny minus CC8,6,9
p <- readRDS("/Users/akhaliq/Desktop/spatial_analysis/pdac_new/ischia/new_k10/hallmark/allcc/ES_allcc.rds") # for hallmark
colnames(p) <- str_remove(colnames(p), "^HALLMARK_")
p <- readRDS("/Users/akhaliq/Desktop/spatial_analysis/pdac_new/ischia/new_k10/minus6_8_9/ES.rds") # for FGES minus CC8,6,9
p <- readRDS("/Users/akhaliq/Desktop/spatial_analysis/pdac_new/ischia/K10_old/scripts/meenakshi/ES.rds") # FGES

p <- readRDS("/Users/akhaliq/Desktop/spatial_analysis/pdac_new/ischia/new_k10/genesets/adel_liding_markers/ES_histology_markers.rds") # adel li ding markers
p <- readRDS("/Users/akhaliq/Desktop/spatial_analysis/pdac_new/ischia/new_k10/genesets/all_markers/all/ES_All_TAMs.rds")# TAMs
p <- readRDS("/Users/akhaliq/Desktop/spatial_analysis/pdac_new/ischia/new_k10/genesets/all_markers/all/ES_All_T Cells.rds") # T cells
p <- readRDS("/Users/akhaliq/Desktop/spatial_analysis/pdac_new/ischia/new_k10/genesets/all_markers/all/ES_All_Basal Classical.rds") # Bassel Clasicall 
p <- readRDS("/Users/akhaliq/Desktop/spatial_analysis/pdac_new/ischia/new_k10/genesets/all_markers/all/ES_All_CAFs.rds") # CAFs
p <- readRDS("/Users/akhaliq/Desktop/spatial_analysis/pdac_new/ischia/new_k10/genesets/all_markers/all/ES_All_B_plasma_Endo.rds") # B, Plasma & Endo
p <- readRDS("/Users/akhaliq/Desktop/spatial_analysis/pdac_new/ischia/new_k10/genesets/all_markers/all/merged_all.rds") # all together

tibble_data <- p %>%
  rownames_to_column(var = "condition") %>%
  pivot_longer(cols = -condition, names_to = "source", values_to = "score") %>%
  arrange(source, condition) %>%
  select(source, condition, score)



a = pdac

a[['icms']] <- tibble_data %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

DefaultAssay(a) <- "icms"
a <- ScaleData(a)

a@assays$icms@data <- a@assays$icms@scale.data
Idents(a) <- a$cc_ischia_10


data <- a

rm(a)

Idents(data) <- data$cc_ischia_10

assay_matrix <- data[["icms"]]@data
norm_weights <- as.data.frame(t(assay_matrix))



# Install and load the required R packages if you haven't already
#install.packages("ggplot2")
#install.packages("gridExtra")
library(ggplot2)
library(gridExtra)
library(Seurat)
library(tidyverse)
library(viridis)


# Define your list of image names and features
image_names <- names(pdac@images)
image_names <- c(
  "IU_PDA_NP11", "IU_PDA_NP2", "IU_PDA_T11", "IU_PDA_T3", "IU_PDA_T8",
  "IU_PDA_HM13", "IU_PDA_HM2", "IU_PDA_HM4", "IU_PDA_HM6", "IU_PDA_HM8",
  "IU_PDA_LNM7", "IU_PDA_T1", "IU_PDA_T4", "IU_PDA_LNM10", "IU_PDA_LNM6",
  "IU_PDA_LNM8", "IU_PDA_HM12", "IU_PDA_T9", "IU_PDA_T2", "IU_PDA_HM10",
  "IU_PDA_HM9", "IU_PDA_NP10", "IU_PDA_T10", "IU_PDA_T12", "IU_PDA_T6",
  "IU_PDA_HM11", "IU_PDA_HM3", "IU_PDA_HM5", "IU_PDA_LNM12", "IU_PDA_HM2_2"
)

features <- c(
  "Androgen", "EGFR", "Estrogen", "Hypoxia", "JAK-STAT", "MAPK",
  "NFkB", "PI3K", "TGFb", "TNFa", "Trail", "VEGF", "WNT", "p53"
)

# Create a function to generate the spatial feature plot for a specific image and feature
generate_spatial_feature_plot <- function(image_name, feature) {
  # Replace this line with your actual code to generate the spatial feature plot
  # For example:
  plot <- SpatialFeaturePlot(pdac, features = feature, stroke = 1, max.cutoff = "q99", min.cutoff = "q1", images = image_name) + scale_fill_viridis(option = "A")
  return(plot)
}

# Create a directory to save the PDF files
dir.create("pdf_output", showWarnings = FALSE)

# Loop through image names and features to generate and save individual PDFs
for (image_name in image_names) {
  for (feature in features) {
    plot <- generate_spatial_feature_plot(image_name, feature)
    pdf_name <- paste0("pdf_output/", image_name, "_", feature, ".pdf")
    pdf(pdf_name)
    print(plot)
    dev.off()
  }
}


###
Without Legened
##

library(ggplot2)
library(gridExtra)
library(viridis)

# Define your list of image names and features
image_names <- c(
  "IU_PDA_NP11", "IU_PDA_NP2", "IU_PDA_T11", "IU_PDA_T3", "IU_PDA_T8",
  "IU_PDA_HM13", "IU_PDA_HM2", "IU_PDA_HM4", "IU_PDA_HM6", "IU_PDA_HM8",
  "IU_PDA_LNM7", "IU_PDA_T1", "IU_PDA_T4", "IU_PDA_LNM10", "IU_PDA_LNM6",
  "IU_PDA_LNM8", "IU_PDA_HM12", "IU_PDA_T9", "IU_PDA_T2", "IU_PDA_HM10",
  "IU_PDA_HM9", "IU_PDA_NP10", "IU_PDA_T10", "IU_PDA_T12", "IU_PDA_T6",
  "IU_PDA_HM11", "IU_PDA_HM3", "IU_PDA_HM5", "IU_PDA_LNM12", "IU_PDA_HM2_2"
)

features <- c(
  "Androgen", "EGFR", "Estrogen", "Hypoxia", "JAK-STAT", "MAPK",
  "NFkB", "PI3K", "TGFb", "TNFa", "Trail", "VEGF", "WNT", "p53"
)

features<- rownames(pdac2@assays$fges2) # hallmark

features<- c("EPCAM", "KRT8", "KRT18", "COL1A1", "CLDN5", "CD3D", "CD79A", "LYZ", "FAP", "PDPN", "PDGFRA", "COL1A2", "CXCL12", "ANTXR1", "SDC1", "SEMA3C", "CD9", "CST1", "TGFB1", "LAMP5", "SCARA5", "DLK1", "GPC3", "RGS5", "CSPG4", "PDGFRB", "CD248", "EPAS1", "CD1C", "C1QA", "C1QB", "MRC1", "S100A8", "S100A9", "APOE", "CD163", "FCGR3B", "FOXP3", "LAG3", "CTLA4", "CCR7", "SELL", "TCF7", "ANXA1", "IL7R", "LMNA", "CXCL13", "GZMK", "GZMB", "IFNG", "PRF1", "CCL4", "CXCL13", "LAG3", "HAVCR2", "CD96", "EOMES", "KLRG1")

features<- c("XBP1", "ECM1", "IL4R", "GADD45B", "CD81", "ANXA4", "BHLHE40", "COL6A1", "SNX9", "AHR", "RHOB", "BCL2L1")
features<- c("ACP5" ,"LIPA" ,"LPL" ,"CCL18" ,"C1QB" ,"APOC1")#Lamp3+Tams
features<- c("CXCL9" ,"CCL4" ,"CXCL10" ,"CD40" ,"IL1RN" ,"CCL4" ,"IL1RN" ,"CXCL10")#Inflam TAMS
features<- c("HLA-DQA1","HLA-DRB5","HLA-DRB1","CD40")#Reg tAM
features<- c("S100A4","CCL4","CTSS","AREG")# Classical TAMS

features<- rownames(a@assays$tams) # hallmark

#CC1 M-35-Fatty.Acid-Acetyl.CoA
 features<- c("M-106-Glucose-Glucose.6.phosphate", "M-5-Pyruvate-Acetyl.Coa", "M-8-Citrate-2OG", "M-13-Malate-Oxaloacetate", "M-37-Aspartate-Aspartate.OUT", "M-5-Pyruvate-Acetyl.Coa", "M-35-Fatty.Acid-Acetyl.CoA", "M-53-Leucine-Acetyl.CoA", "M-150-PRPP-UMP", "M-7-Acetyl.CoA.Oxaloacetate-Citrate", "M-34-Acetyl.CoA-Fatty.Acid", "M-153-UMP-CDP", "M-16-Serine-Pyruvate", "M-12-Fumarate-Malate", "M-4-3PD-Pyruvate", "M-2-G6P-G3P")
#CC3
features<- c("M-106-Glucose-Glucose.6.phosphate", "M-5-Pyruvate-Acetyl.Coa", "M-8-Citrate-2OG", "M-34-Acetyl.CoA-Fatty.Acid", "M-35-Fatty.Acid-Acetyl.CoA")

#CC5 
features<- c("M-48-Glutamate-Glutamine", "M-22-Glycine-Glycine.OUT", "M-41-Aspartate-B.Alanine", "M-37-Aspartate-Aspartate.OUT", "M-49-Glutamate-GABA")


  # Create a directory to save the PDF files
  dir.create("TAMS_genes", showWarnings = FALSE)
    dir.create("TAMS_genes/Classical_TAM", showWarnings = FALSE)

  # Define custom theme for the plots
  custom_theme <- theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),  # Center-align the title
      axis.text = element_blank(),  # Remove axis text
      axis.title = element_blank(),  # Remove axis title
      legend.position = "none",  # Remove legend
      panel.spacing = unit(0.5, "lines")  # Adjust spacing between panels
    )

  # Loop through image names
  for (image_name in image_names) {
    # Create a list to store plots for the current sample
    sample_plots <- list()
    
    for (feature in features) {
      # Replace this line with your actual code to generate the spatial feature plot
      # For example, assuming you have a function SpatialFeaturePlot to create the plot:
     # plot <- SpatialFeaturePlot(pdac, features = feature, stroke = 1, max.cutoff = "q99", min.cutoff = "q1", images = image_name) + scale_fill_viridis(option = "A")
      plot <- SpatialFeaturePlot(a, features = feature, stroke = 1, max.cutoff = "q99", min.cutoff = "q1", images = image_name) + scale_fill_viridis(option = "H")
      # Set the height and width of the plot
      plot <- plot + theme_void()
      plot <- plot + labs(title = feature)
      plot <- plot + custom_theme  # Apply the custom theme
      
      sample_plots <- c(sample_plots, list(plot))
    }
    
    # Calculate the number of pages needed
    num_pages <- ceiling(length(sample_plots) / 9)
    
    # Create a PDF file for the current sample
    pdf_name <- paste0("TAMS_genes/Classical_TAM/", image_name, ".pdf")
    pdf(pdf_name, width = 10, height = 10, onefile = TRUE)  # Create a single PDF file
    
    for (page in 1:num_pages) {
      # Select the plots for the current page
      start_plot <- (page - 1) * 9 + 1
      end_plot <- min(page * 9, length(sample_plots))
      page_plots <- sample_plots[start_plot:end_plot]
      
      # Create a 3x3 grid for the current page
      page_grid <- grid.arrange(grobs = page_plots, ncol = 3)
      
      # Print the page to the PDF
      print(page_grid)
    }
    
    dev.off()
  }

####
With Legend !!!
####

library(ggplot2)
library(gridExtra)

# Define your list of image names and features
image_names <- c(
  "IU_PDA_NP11", "IU_PDA_NP2", "IU_PDA_T11", "IU_PDA_T3", "IU_PDA_T8",
  "IU_PDA_HM13", "IU_PDA_HM2", "IU_PDA_HM4", "IU_PDA_HM6", "IU_PDA_HM8",
  "IU_PDA_LNM7", "IU_PDA_T1", "IU_PDA_T4", "IU_PDA_LNM10", "IU_PDA_LNM6",
  "IU_PDA_LNM8", "IU_PDA_HM12", "IU_PDA_T9", "IU_PDA_T2", "IU_PDA_HM10",
  "IU_PDA_HM9", "IU_PDA_NP10", "IU_PDA_T10", "IU_PDA_T12", "IU_PDA_T6",
  "IU_PDA_HM11", "IU_PDA_HM3", "IU_PDA_HM5", "IU_PDA_LNM12", "IU_PDA_HM2_2"
)
images_names <- names(data@images)

features <- c(
  "Androgen", "EGFR", "Estrogen", "Hypoxia", "JAK-STAT", "MAPK",
  "NFkB", "PI3K", "TGFb", "TNFa", "Trail", "VEGF", "WNT", "p53"
)

features <- colnames(norm_weights)
# Create a directory to save the PDF files
dir.create("pdf_output4", showWarnings = FALSE)

# Define custom theme for the plots with a smaller legend
custom_theme <- theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),  # Center-align the title
    axis.text = element_blank(),  # Remove axis text
    axis.title = element_blank(),  # Remove axis title
    legend.position = "bottom",  # Position legend at the bottom
    legend.box = "horizontal",  # Display legend as horizontal box
    legend.title = element_text(size = 8),  # Adjust legend title size
    legend.text = element_text(size = 8)  # Adjust legend text size
  )

# Loop through image names
for (image_name in image_names) {
  # Create a list to store plots for the current sample
  sample_plots <- list()
  
  for (feature in features) {
    # Replace this line with your actual code to generate the spatial feature plot
    # For example, assuming you have a function SpatialFeaturePlot to create the plot:
    plot <- SpatialFeaturePlot(pdac, features = feature, stroke = 1, max.cutoff = "q99", min.cutoff = "q1", images = image_name) + scale_fill_viridis(option = "H")
    
    # Set the height and width of the plot
    plot <- plot + theme_void()
    plot <- plot + labs(title = feature)
    plot <- plot + custom_theme  # Apply the custom theme
    
    sample_plots <- c(sample_plots, list(plot))
  }
  
  # Calculate the number of pages needed
  num_pages <- ceiling(length(sample_plots) / 6)  # Two rows, three columns
  
  # Create a PDF file for the current sample
  pdf_name <- paste0("pdf_output4/", image_name, ".pdf")
  pdf(pdf_name, width = 10, height = 10, onefile = TRUE)  # Create a single PDF file
  
  for (page in 1:num_pages) {
    # Select the plots for the current page
    start_plot <- (page - 1) * 6 + 1
    end_plot <- min(page * 6, length(sample_plots))
    page_plots <- sample_plots[start_plot:end_plot]
    
    # Create a 2x3 grid for the current page
    page_grid <- grid.arrange(grobs = page_plots, ncol = 3)
    
    # Print the page to the PDF
    print(page_grid)
  }
  
  dev.off()
}

####
Feature Plot output by Feature: Each feature across all Images
####
# Load required libraries

library(Seurat)  # For SpatialFeaturePlot
library(ggplot2)  # For theme_minimal and other ggplot functions
library(viridis)  # For scale_fill_viridis
library(gridExtra)  # For grid.arrange
library(grid)  # For grid.newpage and grid.text

# Define your list of image names and features
image_names <- c(
  "IU_PDA_NP11", "IU_PDA_NP2", "IU_PDA_T11", "IU_PDA_T3", "IU_PDA_T8",
  "IU_PDA_HM13", "IU_PDA_HM2", "IU_PDA_HM4", "IU_PDA_HM6", "IU_PDA_HM8",
  "IU_PDA_LNM7", "IU_PDA_T1", "IU_PDA_T4", "IU_PDA_LNM10", "IU_PDA_LNM6",
  "IU_PDA_LNM8", "IU_PDA_HM12", "IU_PDA_T9", "IU_PDA_T2", "IU_PDA_HM10",
  "IU_PDA_HM9", "IU_PDA_NP10", "IU_PDA_T10", "IU_PDA_T12", "IU_PDA_T6",
  "IU_PDA_HM11", "IU_PDA_HM3", "IU_PDA_HM5", "IU_PDA_LNM12", "IU_PDA_HM2_2"
)
features <- colnames(norm_weights)

# Create a directory to save the PDF files
dir.create("pdf_output_by_feature", showWarnings = FALSE)

# Define custom theme for the plots with a smaller legend
custom_theme <- theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 10),  # Center-align the title and reduce size
    axis.text = element_blank(),  # Remove axis text
    axis.title = element_blank(),  # Remove axis title
    legend.position = "bottom",  # Position legend at the bottom
    legend.box = "horizontal",  # Display legend as horizontal box
    legend.title = element_text(size = 8),  # Adjust legend title size
    legend.text = element_text(size = 8)  # Adjust legend text size
  )

# Loop through features
for (feature in features) {
  # Create a list to store plots for the current feature
  feature_plots <- list()
  
  for (image_name in image_names) {
    # Generate the spatial feature plot
    plot <- SpatialFeaturePlot(pdac, features = feature, stroke = 1, max.cutoff = "q99", min.cutoff = "q1", images = image_name) + 
      scale_fill_viridis(option = "H") +
      theme_void() +
      labs(title = image_name) +
      custom_theme  # Apply the custom theme
    
    feature_plots <- c(feature_plots, list(plot))
  }
  
  # Calculate the number of pages needed
  num_pages <- ceiling(length(feature_plots) / 6)  # Two rows, three columns
  
  # Create a PDF file for the current feature
  pdf_name <- paste0("pdf_output_by_feature/", feature, ".pdf")
  pdf(pdf_name, width = 15, height = 12, onefile = TRUE)  # Increased height to accommodate feature title
  
  # Add a title page with the feature name
  grid.newpage()
  grid.text(feature, x = 0.5, y = 0.5, gp = gpar(fontsize = 40, fontface = "bold"))
  
  for (page in 1:num_pages) {
    # Select the plots for the current page
    start_plot <- (page - 1) * 6 + 1
    end_plot <- min(page * 6, length(feature_plots))
    page_plots <- feature_plots[start_plot:end_plot]
    
    # Create a 2x3 grid for the current page
    page_grid <- grid.arrange(
      grobs = page_plots, 
      ncol = 3,
      top = textGrob(paste("Feature:", feature), gp = gpar(fontsize = 16, fontface = "bold"))
    )
    
    # Print the page to the PDF
    print(page_grid)
  }
  
  dev.off()
}

####

# Heatmaps

suppressPackageStartupMessages(library(escape))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dittoSeq))
suppressPackageStartupMessages(library(ggplot2))



library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(ggplot2)
library(pheatmap)

Idents(a) <- a$annot



data <- a


rm(a)


df <- t(as.matrix(data@assays$icms@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(data)) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score))
#`summarise()` has grouped output by 'cluster'. You can override using the `.groups` argument.
top_acts_mat <- df %>%
  pivot_wider(id_cols = 'cluster', names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)


my_breaks <- c(seq(-2, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 2, length.out=floor(palette_length/2)))





pdf("Heatmap1.pdf")
pheatmap(top_acts_mat, border_color = NA, color=my_color, breaks = my_breaks) 
dev.off()


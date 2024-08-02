  library(spacexr)
  library(Seurat)
  #library(STdeconvolve)

#Data Preprocessing

#loding the Spatial data
#pdac <- readRDS("/Users/akhaliq/Desktop/spatial_analysis/new_object/pdac_cl_renamed_mod.rds")
pdac <- readRDS("/Users/akhaliq/Desktop/spatial_analysis/st_new_seurat/all_new_mod/pdac_most_updated.rds")

# First we load the Single cell data
sc <- readRDS("/Users/akhaliq/Desktop/spatial_analysis/st_new_seurat/all_new_mod/sc_rctd.rds")

#since Ankur's data got decimal values i have to remove it out
#Idents(sc)<-"dataset"
#sc <- subset(sc,idents=c("peng","raghavan","sonia"))

#counts <- read.csv(file.path(refdir,"dge.csv")) # load in counts matrix
# Assuming your Seurat object is named "sc"
sc$celltype_nicheDE <- gsub("T/NK cells", "T_NK cells", sc$celltype_nicheDE)

sc$celltype_nicheDE <- gsub("CAF-S1", "Fibroblasts", sc$celltype_nicheDE)

counts <- sc@assays$RNA@counts

Idents(sc) <- "celltype_nicheDE"
cluster <- as.factor(sc$celltype_nicheDE)
names(cluster) <- colnames(sc)
nUMI <- sc$nCount_RNA

names(nUMI) <- colnames(sc)

# Create the reference Object..

reference <- Reference(counts, cluster, nUMI)

# next we load spatial transcriptomics data.

counts <- pdac@assays$Spatial@counts

# Adding Coordinates for all the 30 images into one data frame

image_names <- c("IU_PDA_HM9", "IU_PDA_HM10", "IU_PDA_HM11", "IU_PDA_HM12", "IU_PDA_HM13", "IU_PDA_HM2", "IU_PDA_HM3", "IU_PDA_HM4", "IU_PDA_HM5", "IU_PDA_HM6", "IU_PDA_HM8", "IU_PDA_LNM10", "IU_PDA_LNM12", "IU_PDA_LNM6", "IU_PDA_LNM7", "IU_PDA_LNM8", "IU_PDA_HM2_2", "IU_PDA_NP10", "IU_PDA_NP11", "IU_PDA_NP2", "IU_PDA_T1", "IU_PDA_T9", "IU_PDA_T10", "IU_PDA_T11", "IU_PDA_T12", "IU_PDA_T2", "IU_PDA_T3", "IU_PDA_T4", "IU_PDA_T6", "IU_PDA_T8")

# Create an empty list to store the coordinates for each image
coordinates_list <- list()

# Loop through each image and generate coordinates
for (image_name in image_names) {
  # Get tissue coordinates for the current image
  pos <- GetTissueCoordinates(pdac, image = image_name)
  
  # Assign column names to the coordinates
  colnames(pos) <- c('x','y')
  
  # Store the coordinates in the list
  coordinates_list[[image_name]] <- pos
}

coords <- do.call(rbind, coordinates_list)


# Remove the pdacle name from row names
row_names <- rownames(coords)
new_row_names <- gsub("^.+\\.", "", row_names)
rownames(coords) <- new_row_names

coords[is.na(colnames(coords))] <- NULL
query <- SpatialRNA(coords, counts, colSums(counts))

print(dim(query@counts))


## Examine SpatialRNA object (optional)
print(dim(query@counts)) # observe Digital Gene Expression matrix
hist(log(query@nUMI,2)) # histogram of log_2 nUMI

print(head(query@coords)) # start of coordinate data.frame
barcodes <- colnames(query@counts) # pixels to be used (a list of barcode names). 

# This list can be restricted if you want to crop the query e.g. 
# query <- restrict_query(query, barcodes) provides a basic plot of the nUMI of each pixel
# on the plot:
#plot_puck_continuous(query, barcodes, query@nUMI, ylimit = c(0,round(quantile(query@nUMI,0.9))),title ='plot of nUMI') 


RCTD <- create.RCTD(query, reference, max_cores = 8)
#RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
RCTD.full <- run.RCTD(RCTD, doublet_mode = "full")


pdac <- AddMetaData(pdac, metadata = RCTD@results$results_df)

# RCTD in Full Mode
RCTD.full <- run.RCTD(RCTD, doublet_mode = "full")

barcodes <- colnames(pdac@assays$Spatial@counts)
weights <- RCTD.full@results$weights

#We next normalize the weights using normalize_weights so that they sum to one. Each entry represents the estimated proportion of each cell type on each pixel.

norm_weights <- normalize_weights(weights)

cell_types <- c("B cells", "Endothelial cells","Hepatocytes","PVL","Tumor Epithelial cells","Fibroblasts","Monocytes","T_NK cells","Normal Epithelial cells","TAM")
print(head(norm_weights[,cell_types])) # observe weight values

# RCTD in Multi Mode

RCTD.multi <- run.RCTD(RCTD, doublet_mode = "multi")


#Plotting

png("barplots_FirstType_orig.ident.png", units="in", width=15, height=10, res=300)
dittoBarPlot(pdac,"orig.ident",group.by = "first_type")
dev.off()


png("barplots_orig.ident_FirstType.png", units="in", width=15, height=10, res=300)
dittoBarPlot(pdac,"first_type",group.by = "orig.ident")
dev.off()


png("barplots_SecondType_orig.ident.png", units="in", width=15, height=10, res=300)
dittoBarPlot(pdac,"orig.ident",group.by = "second_type")
dev.off()


png("barplots_orig.ident_SecondType.png", units="in", width=15, height=10, res=300)
dittoBarPlot(pdac,"second_type",group.by = "orig.ident")
dev.off()



###
https://indiana.sharepoint.com/:b:/r/sites/O365-MasoodLab/Shared%20Documents/General/Ateeq_khaliq_Asst.%20Professor/Primary_PDAC%20Vs.%20Mets/RCTD/spatial_Pie_plots/spatial_plots_pdac_Tier2.pdf?csf=1&web=1&e=oeQ4rv

hc <- hclust(dist(as.data.frame(df)), 'ward.D2')
topicCols  <- rainbow(nrow(as.data.frame(df)))
names(topicCols ) <- rownames(as.data.frame(df))[hc$order]
topicCols  <- topicCols [rownames(as.data.frame(df))]



pos <- GetTissueCoordinates(pdac, image = "IU_PDA_HM3")
  colnames(pos) <- c('x','y')
  
  # Generate the plot for the current image
  plt <- vizAllTopics(theta = as.data.frame(df12),
                      pos = pos12,
                      r = 3,
                      topicCols=topicCols,
                      lwd = 0.05,
                      showLegend = TRUE,
                      plotTitle = NA) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 2)) +
    ## outer border
    ggplot2::geom_rect(data = data.frame(pos),
                       ggplot2::aes(xmin = min(x) - 90, xmax = max(x) + 90,
                                    ymin = min(y) - 90, ymax = max(y) + 90),
                       fill = NA, color = "black", linetype = "solid", size = 0.2) +
    ggplot2::theme(
      plot.background = ggplot2::element_blank()
    ) +
    ## remove the pixel "groups", which is the color aesthetic for the pixel borders
    ggplot2::guides(colour = "none")


ggplot(agg_data_percent_df, aes(x = orig.ident, y = percentage, fill = cell_type)) geom_bar(stat = "identity", position = "fill") abs(title = "Percentage Stacked Bar Plot", x = "pdacle ID", y = "Percentage", fill = "Cell Type") theme_minimal() theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5))




### Adding the decon to seurat
# lets add the RTCD.full

pdac@assays[["rctd_fullfinal"]] <- CreateAssayObject(data = t(as.matrix(norm_weights)))


# Seems to be a bug in SeuratData package that the key is not set and any
# plotting function etc. will throw an error.
if (length(pdac@assays$rctd_fullfinal@key) == 0) {
    pdac@assays$rctd_fullfinal@key = "rctd_fullfinal_"
}

"B cells", "CAF-S1", "CAF-S4", "Endothelial cells", "Epithelial cells", "Hepatocytes", "Myeloid", "NK cells", "T cells" 

DefaultAssay(pdac) <- "rctd_fullfinal"
SpatialFeaturePlot(pdac, features = c("B cells", "CAF-S1", "CAF-S4", "Endothelial cells", "Epithelial cells", "Hepatocytes", "Myeloid", "NK cells", "T cells" ), pt.size.factor = 1.6, ncol = 2, crop = TRUE,images="IU_PDA_T9")


# Define the image names
imgs <- c("IU_PDA_HM10", "IU_PDA_HM11", "IU_PDA_HM12", "IU_PDA_HM13", "IU_PDA_HM14", "IU_PDA_HM3", "IU_PDA_HM4", "IU_PDA_HM5", "IU_PDA_HM6", "IU_PDA_HM7", "IU_PDA_HM9", "IU_PDA_LNM11", "IU_PDA_LNM13", "IU_PDA_LNM7", "IU_PDA_LNM8", "IU_PDA_LNM9", "IU_PDA_NH3", "IU_PDA_NP11", "IU_PDA_NP12", "IU_PDA_NP2", "IU_PDA_T1", "IU_PDA_T10", "IU_PDA_T11", "IU_PDA_T12", "IU_PDA_T13", "IU_PDA_T2", "IU_PDA_T4", "IU_PDA_T5", "IU_PDA_T7", "IU_PDA_T9")

# Define the features
features <- c("B cells", "CAF-S1", "CAF-S4", "Endothelial cells", "Epithelial cells", "Hepatocytes", "Myeloid", "NK cells", "T cells")

# Loop through each image and generate a plot with all features
for (img in imgs) {
  # Generate the plot
  plot <- SpatialFeaturePlot(pdac, features = features, pt.size.factor = 1.6, ncol = 3, crop = TRUE, images = img)

  # Add a title to the plot
  plot.new()
  title(main = img)

  # Save the plot with higher resolution
  png(paste0("plots/", img, ".png"), width = 1200, height = 800)
  print(plot)
  dev.off()
}

# Convert the Multi Mode proportions into dataframe

# Create an empty list to store the sub_weights values
sub_weights_list <- list()

# Loop over each element of RCTD.multi@results
for (i in seq_along(RCTD.multi@results)) {
  # Extract the sub_weights value and append it to the list
  sub_weights_list[[i]] <- RCTD.multi@results[[i]]$sub_weights
}

head(sub_weights_list)

###

# Define a function to convert the sub_weights list to a data frame
sub_weights_to_df <- function(sub_weights_list) {
  # Get all column names from the list
  all_col_names <- unique(unlist(lapply(sub_weights_list, names)))
  
  # Create an empty data frame with all column names
  df <- data.frame(matrix(0, ncol = length(all_col_names), nrow = length(sub_weights_list)))
  colnames(df) <- all_col_names
  
  # Fill in the data frame with the sub_weights values
  for (i in seq_along(sub_weights_list)) {
    for (j in seq_along(names(sub_weights_list[[i]]))) {
      col_name <- names(sub_weights_list[[i]])[j]
      df[i, col_name] <- sub_weights_list[[i]][[j]]
    }
  }
  
  return(df)
}
df <- sub_weights_to_df(sub_weights_list)

# Set the rownames of df to rownames(RCTD.multi@spatialRNA@coords)
rownames(df) <- rownames(RCTD.multi@spatialRNA@coords)
head(df)




### Adding the decon to seurat
# lets add the RTCD.multi

pdac@assays[["rctd_multi"]] <- CreateAssayObject(data = t(as.matrix(df)))


# Seems to be a bug in SeuratData package that the key is not set and any
# plotting function etc. will throw an error.
if (length(pdac@assays$rctd_multi@key) == 0) {
    pdac@assays$rctd_multi@key = "rctd_nicheDE_"
}

# "B cells", "CAF-S1", "CAF-S4", "Endothelial cells", "Epithelial cells", "Hepatocytes", "Myeloid", "NK cells", "T cells" 

 [1] "CAF-S1"                  "Monocytes"              
 [3] "B cells"                 "TAM"                    
 [5] "Hepatocytes"             "Endothelial cells"      
 [7] "PVL"                     "Normal Epithelial cells"
 [9] "T_NK cells"              "Tumor Epithelial cells" 

DefaultAssay(pdac) <- "rctd_multi"
SpatialFeaturePlot(pdac, features = c("CAF-S1", "B cells",  "Hepatocytes", "PVL", "T_NK cells", "Monocytes", "TAM", "Endothelial cells","Normal Epithelial cells", "Tumor Epithelial cells" ), pt.size.factor = 1.6, ncol = 2, crop = TRUE,images="IU_PDA_T9")


# Define the image names
imgs <- c("IU_PDA_T1", "IU_PDA_NP2", "IU_PDA_T2", "IU_PDA_NH2", "IU_PDA_HM2", "IU_PDA_HM3", "IU_PDA_T3", "IU_PDA_HM4", "IU_PDA_T4", "IU_PDA_HM5", "IU_PDA_HM6", "IU_PDA_LNM6", "IU_PDA_T6", "IU_PDA_LNM7", "IU_PDA_HM8", "IU_PDA_LNM8", "IU_PDA_T8", "IU_PDA_T9", "IU_PDA_HM9", "IU_PDA_HM10", "IU_PDA_LNM10", "IU_PDA_NP10", "IU_PDA_T10", "IU_PDA_HM11", "IU_PDA_T11", "IU_PDA_NP11", "IU_PDA_HM12", "IU_PDA_LNM12", "IU_PDA_T12", "IU_PDA_HM13")

# Define the features
features <- c("B cells", "CAF-S1", "CAF-S4", "Endothelial cells", "Epithelial cells", "Hepatocytes", "Myeloid", "NK cells", "T cells")

# Loop through each image and generate a plot with all features
for (img in imgs) {
  # Generate the plot
  plot <- SpatialFeaturePlot(pdac1, features = features, pt.size.factor = 1.6, ncol = 3, crop = TRUE, images = img)

  # Add a title to the plot
  plot.new()
  title(main = img)

  # Save the plot with higher resolution
  png(paste0(img, "_multi.png"), width = 1200, height = 800)
  print(plot)
  dev.off()
}

#norm_weights_multi <- normalize_weights(df)

#to create barplot for proportions
# run /Users/akhaliq/Desktop/rctd/codes/aggregating_proportions.py



# replace the old patterns with the new ones
pdac@meta.data$orig.ident <- gsub(
  pattern = paste0("(", paste(old_patterns, collapse = "|"), ")"),
  replacement = function(x) new_patterns[match(x, old_patterns)],
  x = pdac@meta.data$orig.ident
)


## plotting scatterpies


library(ggplot2)
library(SPOTlight)
library(SingleCellExperiment)
library(SpatialExperiment)
library(scater)
library(scran)
library(gridExtra)
library(grid)

pos <- GetTissueCoordinates(pdac, image = "IU_PDA_HM9")
df12 <- subset(df, rownames(df) %in% rownames(pos))

plotSpatialScatterpie(x = pos, y = df12)


# Create a PDF device to save the plots
pdf("spatial_plots_mod_col_new_tier1.pdf")
df<- as.matrix(norm_weights)
# List of image names
image_names <- c("IU_PDA_HM9", "IU_PDA_HM10", "IU_PDA_HM11", "IU_PDA_HM12", "IU_PDA_HM13", "IU_PDA_HM2",
                 "IU_PDA_HM3", "IU_PDA_HM4", "IU_PDA_HM5", "IU_PDA_HM6", "IU_PDA_HM8", "IU_PDA_LNM10",
                 "IU_PDA_LNM12", "IU_PDA_LNM6", "IU_PDA_LNM7", "IU_PDA_LNM8", "IU_PDA_NH2", "IU_PDA_NP10",
                 "IU_PDA_NP11", "IU_PDA_NP2", "IU_PDA_T1", "IU_PDA_T9", "IU_PDA_T10", "IU_PDA_T11",
                 "IU_PDA_T12", "IU_PDA_T2", "IU_PDA_T3", "IU_PDA_T4", "IU_PDA_T6", "IU_PDA_T8")


# Set a beautiful color palette
#paletteMartin <- c("#E6194B", "#3CB44B", "#4363D8", "#F58231", "#911EB4", "#42D4F4", "#F032E6", "#469990", "#9A6324", "#800000")
paletteMartin <- c("#7FC97F", "#BEAED4", "#FDC086", "#FF0000", "#386CB0", "#F0027F", "#BF5B17", "#666666", "#1B9E77", "#D95F02")
#paletteMartin <- c("#7FC97F", "#BEAED4", "#FDC086", "#FF0000", "#386CB0", "#FF1493", "#FFD700", "#666666", "#1B9E77", "#FF7F00")
#paletteMartin <- c("#006400", "#BEAED4", "#FDC086", "#FF0000", "#386CB0", "#F0027F", "#BF5B17", "#666666", "#1B9E77", "#D95F02")
#paletteMartin <- c("#FF0000", "#FF6600", "#FFCC00", "#00FF00", "#0066FF", "#CC00FF", "#FF00FF", "#00FFFF", "#FF99CC", "#9900FF")
df12<- as.matrix(norm_weights)
pal <- colorRampPalette(paletteMartin)(length(df12))
names(pal) <- colnames(df12)

# Iterate over each image
for (image_name in image_names) {
  # Get tissue coordinates for the current image
  pos <- GetTissueCoordinates(pdac, image = image_name)
  
  # Subset the data frame based on row names
  df12 <- subset(df, rownames(df) %in% rownames(pos))
  
  # Change the column name from "CAF-S1" to "CAF"
  #colnames(df12)[colnames(df12) == "CAF"] <- "Fibroblasts"
  
  # Create scatterpie plot
  scatterplot <- plotSpatialScatterpie(x = pos, y = df12, cell_types = colnames(df12), img = FALSE, scatterpie_alpha = 1, pie_scale = 0.4) +
    scale_fill_manual(values = pal, breaks = names(pal)) +
    labs(title = image_name) +
    theme_bw() +
    theme(plot.title = element_text(size = 14, face = "bold"),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  
  # Save each plot in a separate page of the PDF
  print(scatterplot)
}
theme(axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank())

# Close the PDF device
dev.off()

##### Progeny

library(progeny)
library(dplyr)
library(Seurat)
library(ggplot2)
library(tidyr)
library(readr)
library(pheatmap)
library(tibble)


CellsClusters <- data.frame(Cell = names(Idents(pdac)), CellType = as.character(Idents(pdac)), stringsAsFactors = FALSE)


## We compute the Progeny activity scores and add them to our Seurat object
## as a new assay called Progeny. 
pdac_progeny <- progeny(pdac, scale=FALSE, organism="Human", top=1000, perm=1, return_assay = TRUE,assay="Spatial")

## We can now directly apply Seurat functions in our Progeny scores. 
## For instance, we scale the pathway activity scores. 
pdac_progeny<- Seurat::ScaleData(pdac_progeny, assay = "progeny") 

## We transform Progeny scores into a data frame to better handling the results
progeny_scores_df <- 
    as.data.frame(t(GetAssayData(pdac_progeny, slot = "scale.data", 
        assay = "progeny"))) %>%
    rownames_to_column("Cell") %>%
    gather(Pathway, Activity, -Cell) 

## We match Progeny scores with the cell clusters.
progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

## We summarize the Progeny scores by cellpopulation
summarized_progeny_scores <- progeny_scores_df %>% 
    group_by(Pathway, CellType) %>%
    summarise(avg = mean(Activity), std = sd(Activity))



## We prepare the data for the plot
summarized_progeny_scores_df <- summarized_progeny_scores %>%
    dplyr::select(-std) %>%   
    spread(Pathway, avg) %>%
    data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 


paletteLength = 100
myColor = colorRampPalette(c("#008080", "white","#FFA500"))(paletteLength) #Teal, White, and Orange:
#myColor = colorRampPalette(c("#FF00FF", "white","#00FFFF"))(paletteLength) #Magenta, White, and Cyan:


progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0, 
                      length.out=ceiling(paletteLength/2) + 1),
                  seq(max(summarized_progeny_scores_df)/paletteLength, 
                      max(summarized_progeny_scores_df), 
                      length.out=floor(paletteLength/2)))
progeny_hmap = pheatmap(t(summarized_progeny_scores_df[,-1]),fontsize=14, 
                        fontsize_row = 10, 
                        color=myColor, breaks = progenyBreaks, 
                        main = "PROGENy (1000)", angle_col = 45,
                        treeheight_col = 0,  border_color = NA)

#### C-SIDE


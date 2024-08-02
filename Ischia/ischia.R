# Loading required packages
library(ISCHIA)
library(robustbase)
library(data.table)
library(ggplot2)
library(Seurat)
library(dplyr)

# calculating the K value for ischia

library(factoextra)
library(cluster)

pdac <- readRDS("/Users/akhaliq/Desktop/spatial_analysis/pdac_new/pdac_mets_rctd.rds")

## to get Assays from pdac object 

assay_matrix <- pdac[["rctd_tier1"]]@data
norm_weights <- as.data.frame(t(assay_matrix))
#write.csv(assay_df, file = "assay_rctd_full.csv")

## Elbwo plot

set.seed(123)


# Generate random data
num_points <- 200

# Function to compute total within-cluster sum of squares
wss <- function(k) {
  kmeans(norm_weights, k, nstart = 10)$tot.withinss
}

# Compute and plot wss for k = 1 to k = 20
k.values <- 1:20
wss_values <- sapply(k.values, wss)

# Plot the wss values
pdf("wss_values_plot.pdf")
plot(k.values, wss_values,
     type = "b", pch = 19, frame = FALSE, 
     xlab = "Number of clusters K",
     ylab = "Total within-cluster sum of squares")

# Find the elbow point
diff_wss <- diff(wss_values)
curvature <- diff(diff_wss)
elbow_index <- which.max(curvature) + 1

# Add an elbow line
elbow_x <- k.values[elbow_index]
elbow_y <- wss_values[elbow_index]
abline(0, elbow_y, h = elbow_y, col = "red", lty = "dashed")
abline(v = elbow_x, col = "red", lty = "dashed")
dev.off()

# Perform K-means clustering with the optimal K
k_optimal <- elbow_x
clusters <- kmeans(norm_weights, k_optimal, nstart = 10)

# Plot the K-means clusters
plot(norm_weights, col = clusters$cluster, pch = 19, frame = FALSE,
     xlab = "X", ylab = "Y", main = "K-means Clustering")

# Add cluster centroids
points(clusters$centers, col = 1:k_optimal, pch = 8, cex = 2)


***
## Calculate Gap Statistics with Silhoutte width validation,  for cluster validation and selection
set.seed(123)

# Generate random data
num_points <- 200


# Function to compute total within-cluster sum of squares
wss <- function(k) {
  kmeans(norm_weights, k, nstart = 10)$tot.withinss
}

# Compute and plot wss for k = 1 to k = 20
k.values <- 1:20
wss_values <- sapply(k.values, wss)

# Plot the wss values
plot(k.values, wss_values,
     type = "b", pch = 19, frame = FALSE, 
     xlab = "Number of Clusters (K)",
     ylab = "Total Within-Cluster Sum of Squares",
     main = "Elbow Method: Determining Optimal K")

# Find the elbow point
diff_wss <- diff(wss_values)
curvature <- diff(diff_wss)
elbow_index <- which.max(curvature) + 1

# Add an elbow line
elbow_x <- k.values[elbow_index]
elbow_y <- wss_values[elbow_index]
abline(0, elbow_y, h = elbow_y, col = "red", lty = "dashed")
abline(v = elbow_x, col = "red", lty = "dashed")

# Compute and plot Gap statistic for k = 1 to k = 20
gap_stat <- function(k) {
  km.res <- kmeans(norm_weights, k, nstart = 10)
  if (k == 1) {
    return(NA)
  }
  obs_disp <- sum(km.res$withinss)
  
  # Generate reference null distribution
  reference_disp <- replicate(10, {
    random_points <- matrix(rnorm(num_points * 2), ncol = 2)
    km.null <- kmeans(random_points, k, nstart = 10)
    sum(km.null$withinss)
  })
  reference_disp <- mean(reference_disp)
  
  # Calculate Gap statistic
  gap_statistic <- log(reference_disp) - log(obs_disp)
  return(gap_statistic)
}

k.values <- 1:20
gap_stat_values <- sapply(k.values, gap_stat)

# Plot the Gap statistic
plot(k.values, gap_stat_values,
     type = "b", pch = 19, frame = FALSE, 
     xlab = "Number of Clusters (K)",
     ylab = "Gap Statistic",
     main = "Gap Statistic: Determining Optimal K")

# Find the elbow point
elbow_index_gap <- which.max(gap_stat_values)

# Add an elbow line
elbow_x_gap <- k.values[elbow_index_gap]
elbow_y_gap <- gap_stat_values[elbow_index_gap]
abline(0, elbow_y_gap, h = elbow_y_gap, col = "blue", lty = "dashed")
abline(v = elbow_x_gap, col = "blue", lty = "dashed")


###
 ##### Find the optimal number of clusters based on the highest Calinski-Harabasz Index
 set.seed(123)

num_points <- 200

# Function to compute total within-cluster sum of squares
wss <- function(k) {
  kmeans(norm_weights, k, nstart = 10)$tot.withinss
}

# Function to calculate the Calinski-Harabasz Index
calinski_harabasz_index <- function(data, labels) {
  num_clusters <- length(unique(labels))
  num_points <- nrow(data)
  
  # Calculate the centroid of each cluster
  centroids <- tapply(data, labels, FUN = colMeans)
  
  # Calculate the between-cluster dispersion
  between_disp <- 0
  for (i in 1:num_clusters) {
    cluster_points <- data[labels == i, ]
    cluster_centroid <- centroids[i, ]
    between_disp <- between_disp + sum((colMeans(cluster_points) - cluster_centroid) ^ 2) * nrow(cluster_points)
  }
  
  # Calculate the within-cluster dispersion
  within_disp <- 0
  for (i in 1:num_clusters) {
    cluster_points <- data[labels == i, ]
    cluster_centroid <- centroids[i, ]
    within_disp <- within_disp + sum((cluster_points - cluster_centroid) ^ 2)
  }
  
  # Calculate the Calinski-Harabasz Index
  ch_index <- (between_disp / (num_clusters - 1)) / (within_disp / (num_points - num_clusters))
  return(ch_index)
}

# Create a PDF file for saving the plots
pdf("cluster_evaluation_plots.pdf")

# Compute and plot wss for k = 1 to k = 20
k.values <- 1:20
wss_values <- sapply(k.values, wss)

# Plot the wss values
plot(k.values, wss_values,
     type = "b", pch = 19, frame = FALSE, 
     xlab = "Number of Clusters (K)",
     ylab = "Total Within-Cluster Sum of Squares",
     main = "Elbow Method: Determining Optimal K")

# Find the elbow point
diff_wss <- diff(wss_values)
curvature <- diff(diff_wss)
elbow_index <- which.max(curvature) + 1

# Add an elbow line
elbow_x <- k.values[elbow_index]
elbow_y <- wss_values[elbow_index]
abline(0, elbow_y, h = elbow_y, col = "red", lty = "dashed")
abline(v = elbow_x, col = "red", lty = "dashed")

# Compute and plot Gap statistic for k = 1 to k = 20
gap_stat <- function(k) {
  km.res <- kmeans(norm_weights, k, nstart = 10)
  if (k == 1) {
    return(NA)
  }
  obs_disp <- sum(km.res$withinss)
  
  # Generate reference null distribution
  reference_disp <- replicate(10, {
    random_points <- matrix(rnorm(num_points * 2), ncol = 2)
    km.null <- kmeans(random_points, k, nstart = 10)
    sum(km.null$withinss)
  })
  reference_disp <- mean(reference_disp)
  
  # Calculate Gap statistic
  gap_statistic <- log(reference_disp) - log(obs_disp)
  return(gap_statistic)
}

dev.off()

gap_stat_values <- sapply(k.values, gap_stat)

pdf("2_GapStats_plot.pdf")
# Plot the Gap statistic
plot(k.values, gap_stat_values,
     type = "b", pch = 19, frame = FALSE, 
     xlab = "Number of Clusters (K)",
     ylab = "Gap Statistic",
     main = "Gap Statistic: Determining Optimal K")


# Find the elbow point
elbow_index_gap <- which.max(gap_stat_values)

# Add an elbow line
elbow_x_gap <- k.values[elbow_index_gap]
elbow_y_gap <- gap_stat_values[elbow_index_gap]
abline(0, elbow_y_gap, h = elbow_y_gap, col = "blue", lty = "dashed")
abline(v = elbow_x_gap, col = "blue", lty = "dashed")

dev.off()

# Calculate the Calinski-Harabasz Index for each value of k
ch_values <- sapply(k.values, function(k) {
  km.res <- kmeans(norm_weights, k, nstart = 10)
  calinski_harabasz_index(norm_weights, km.res$cluster)
})

# Plot the Calinski-Harabasz Index
plot(k.values, ch_values,
     type = "b", pch = 19, frame = FALSE,
     xlab = "Number of Clusters (K)",
     ylab = "Calinski-Harabasz Index",
     main = "Calinski-Harabasz Index: Determining Optimal K")

# Find the optimal number of clusters based on the highest Calinski-Harabasz Index
optimal_k_ch <- k.values[which.max(ch_values)]

# Add a point for the optimal k
points(optimal_k_ch, max(ch_values), col = "red", pch = 19)
text(optimal_k_ch, max(ch_values), paste("K =", optimal_k_ch), pos = 3, offset = 1, col = "red")

# Close the PDF device
dev.off()


          1 
Optimal K (Elbow Method): 4 
Optimal K (Calinski-Harabasz Index): 2 
Optimal K (Gap Statistic): 20 

#### Ischea

# Loading required packages
library(ISCHIA)
library(robustbase)
library(data.table)
library(ggplot2)
library(Seurat)
library(dplyr)

## to get Assays from pdac object 

#assay_matrix <- pdac[["rctd_fullfinal"]]@data
#norm_weights <- as.data.frame(t(assay_matrix))
#write.csv(assay_df, file = "assay_rctd_full.csv")
deconv.mat= norm_weights
# Deciding about the k

pdf("3_ccPlot.pdf")
Composition.cluster.k(deconv.mat, 20)
dev.off()

# Composition clustering of the deconvoluted spatial spots
pdac <- Composition.cluster(pdac,norm_weights, 12)
#not good...Idents(pdac) <- "CompositionCluster_CC"

table(pdac$CompositionCluster_CC)

pdac$cc_12<-pdac$CompositionCluster_CC


#SpatialDimPlot(pdac, group.by = c("CompositionCluster_CC")) + scale_fill_manual(values = c("cyan", "orange", "purple","green","yellow","blue", "red","black"))

## dim for all images in a loop
"#7FC97F", "#8F30A1", "#FDC086", "#BF5B17", "#17DEEE", "#F0027F", "#666666", "#FFD700", "#1B9E77", "#b9bf17", "#FF00FF", "#377EB8", "#4DAF4A", "#ff9a00", "#FF0000"

# Load required libraries
library(ggplot2)
library(showtext)
library(gridExtra)

# Set plot theme
theme_set(theme_minimal(base_family = "Arial"))

# Load bold title font
font_add_google("Montserrat", "Montserrat-Bold")

# Enable custom font
showtext_auto()

# Define image names
image_names <- c(
  "IU_PDA_T1", "IU_PDA_T2", "IU_PDA_HM2", "IU_PDA_HM2_2", "IU_PDA_NP2", "IU_PDA_T3", 
  "IU_PDA_HM3", "IU_PDA_T4", "IU_PDA_HM4", "IU_PDA_HM5", "IU_PDA_T6", "IU_PDA_HM6", 
  "IU_PDA_LNM6", "IU_PDA_LNM7", "IU_PDA_T8", "IU_PDA_HM8", "IU_PDA_LNM8", "IU_PDA_T9", 
  "IU_PDA_HM9", "IU_PDA_T10", "IU_PDA_HM10", "IU_PDA_LNM10", "IU_PDA_NP10", "IU_PDA_T11", 
  "IU_PDA_HM11", "IU_PDA_NP11", "IU_PDA_T12", "IU_PDA_HM12", "IU_PDA_LNM12", "IU_PDA_HM13"
)

# Custom color palette
paletteMartin <- c(
  "#FF7F0E", "#9467BD", "#69B3E2", "#FFD700", "#238B45",
  "#FF1493", "#7C8B8E", "#E41A1C", "#4B0082"
)

#https://sashamaps.net/docs/resources/20-colors/
#paletteMartin <- c("#ddc2f0", "#ff7f0e", "#2ca02c", "#d62728", "#000000", "#d4f739", "#00ff00", "#09bfe8","#1CFFBB","#F296EC","","")
#paletteMartin <- c("#ddc2f0", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#1f77b4", "#d6616b")
#paletteMartin <- c("#E74C3C", "#3498DB", "#2ECC71", "#9B59B6", "#F39C12", "#1ABC9C", "#E67E22", "#2980B9", "#27AE60", "#D35400", "#8E44AD", "#C0392B")
#paletteMartin <- c("#734F96", "#3498DB", "#2ECC71", "#64a649", "#F39C12", "#1ABC9C", "#FF5733", "#34A853", "#FFC300", "#0074D9", "#FF006E", "#009fef")
paletteMartin <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff')
#c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666", "#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#cab2d6", "#ff7f00", "#b15928", "#33a02c")
#paletteMartin <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666", "#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f")
# Define colors for all 8 CCs
all_ccs <- unique(pdac$CompositionCluster_CC)
color_mapping <- setNames(paletteMartin[1:length(all_ccs)], all_ccs)


# Output PDF file
output_file <- "4_spatial_plots_K12_new.pdf"

# Create PDF file
pdf(output_file, width = 10, height = 7, onefile = TRUE)

# Iterate over image names
for (image_name in image_names) {
  # Create the plot
  plot <- SpatialDimPlot(pdac, group.by = c("CompositionCluster_CC"), images = image_name) +
    scale_fill_manual(values = color_mapping) +
    theme(
      plot.title = element_text(face = "bold", size = 14, margin = margin(b = 10)),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "right",
      legend.text = element_text(size = 10),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white"),
      panel.border = element_rect(color = "black", fill = NA),
      axis.line = element_line(color = "black")
    )
  
  # Add bold title
  plot <- plot + ggtitle(bquote(bold(.(image_name))))
  
  # Save plot on a separate PDF page
  print(plot)
}

# Close PDF file
dev.off()

########

#> Scale for fill is already present.
#> Adding another scale for fill, which will replace the existing scale.


#Composition_cluster_enrichedCelltypes(pdac,"CC4", as.matrix(norm_weights))

### for All CCs

#c("CC1", "CC2", "CC3", "CC4", "CC5", "CC6", "CC7", "CC8", "CC9", "CC10", "CC11","CC12","CC13","CC14")

library(gridExtra)

# Define a function to generate and save the plot for a specific CC
save_cc_plot <- function(cc) {
  # Generate the plot for the given CC
  plot <- Composition_cluster_enrichedCelltypes(pdac, cc, as.matrix(norm_weights))
  
  # Create a PDF file for the current CC plot
  pdf_name <- paste0(cc, ".pdf")
  pdf(file = pdf_name)
  
  # Save the plot to the PDF
  print(plot)
  
  # Close the PDF file
  dev.off()
}

# Iterate over each CC and save its respective plot
ccs <- paste0("CC", 1:12)  # Replace with the appropriate range for your case
for (cc in ccs) {
  save_cc_plot(cc)
}

library(pdftools)

# Define the names of the individual PDF files
pdf_files <- paste0("CC", 1:12, ".pdf")  # Replace with the appropriate range for your case

# Define the name of the merged PDF file
merged_pdf <- "5_enrichedCelltypes_CC_12.pdf"

# Use the pdf_combine function to merge the PDFs
pdf_combine(pdf_files, output = merged_pdf)


pdac.umap <- Composition_cluster_umap(pdac, norm_weights)
#> Plotting scatterpies for 2185 pixels with 20 cell-types...this could take a while if the dataset is large.
pdac.umap$umap.cluster.gg

pdf("6_pie_chart.pdf")
pdac.umap$umap.deconv.gg
dev.off()



#####
## adding UMAP and TSNE to the seurat obj

#write.csv(emb,"tsne_stdcon.csv")
#write.csv(emb.umap,"umap_stdcon.csv")
write.csv(pdac.umap$umap.table,"umap_ischia14.csv")
#emb.tsne <- read.csv("tsne_stdcon.csv",header=T,sep=",",row.names=1)
emb.umap<- read.csv("umap_ischia14.csv",header=T,sep=",",row.names=1)


emb.umap = pdac.umap$umap.table
emb.umap$CompositionCluster_CC <- NULL
emb.umap$Slide <- NULL

emb.umap <- as.matrix(emb.umap)

#colnames(emb.tsne)<- c("tSNE1","tSNE2")
colnames(emb.umap)<- c("UMAP1","UMAP2")
head(emb.umap)



#pdac[['tsne.std']] <- CreateDimReducObject(embeddings = as.matrix(emb.tsne), key = 'tsne_', assay = 'integrated')
#pdac[['umap.std']] <- CreateDimReducObject(embeddings = as.matrix(emb.umap1), key = 'umap_', assay = 'integrated')

#pdac@assays$@key <- "rna_"

pdac[['umap.ischia8']] <- CreateDimReducObject(embeddings = emb.umap, key = 'umap.ischia8_', assay = 'rctd_tier1')


emb.umap<- read.csv("umap_ischia12.csv",header=T,sep=",",row.names=1)
emb.umap$x <- NULL
emb.umap$y <- NULL
emb.umap$Slide <- NULL
colnames(emb.umap) <- "cc_12"
pdac <- AddMetaData(pdac,emb.umap)

pdf("seurat_ischia_umap_12.pdf")
DimPlot(pdac, reduction = "umap.ischia12", label = FALSE,group.by="cc_12")
dev.off()

pdf("barplot_SampVsorig_12.pdf", height=12, width=20)
dittoBarPlot(pdac, "orig.ident", group.by = "cc_12")
dev.off()

pdf("barplot_origVsSamp_12.pdf", height=10, width=20)
dittoBarPlot(pdac, "cc_12", group.by = "orig.ident")
dev.off()


per <- dittoBarPlot(pdac, "orig.ident", group.by = "cc_12", data.out = TRUE)
write.csv(per$data,"percent_K12.csv")

CC4.celltype.cooccur <- spatial.celltype.cooccurence(spatial.object=pdac,deconv.prob.mat=deconv.mat, COI="CC4", prob.th= 0.05, Condition=unique(pdac$orig.ident))
plot.celltype.cooccurence(CC4.celltype.cooccur)

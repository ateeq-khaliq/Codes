library(ggplot2)
library(SPOTlight)
library(SingleCellExperiment)
library(SpatialExperiment)
library(scater)
library(scran)

dir.create("correlation")
setwd("/Users/akhaliq/Desktop/spatial_analysis/pdac_new/ischia/new_k10/correlation")

assay_matrix <- pdac[["rctd_fullfinal"]]@data
norm_weights <- as.data.frame(t(assay_matrix))

mat <- merge(pdac$cc_ischia_10,norm_weights, by = 0)

#mat1 <- t(mat)
# Assuming mat1 is your dataframe

#new_colnames <- as.character(mat1[2,])  # Get the values of the first row as new column names
#mat1 <- mat1[-1,]  # Remove the first row to keep only the data
#mat1 <- mat1[-1,]

# Set the new column names
#colnames(mat1) <- new_colnames

aggregated_df <- mat %>% group_by(x) %>% summarise(across(where(is.numeric), sum, na.rm = TRUE))
aggregated_df = as.data.frame(aggregated_df)

mat2 = aggregated_df
#rm(aggregated_df)

rownames(mat2) <- mat2$x
mat2$x <- NULL

mat3 <- t(mat2)


# Plotting from Spotlight
pdf("test4.pdf", width = 10, height = 10)
plotCorrelationMatrix(as.matrix(mat3))
dev.off()


# # Plotting general 1
cormat <- round(cor(mat3),2)


# Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }

  upper_tri <- get_upper_tri(cormat)

pdf("test2.pdf", width = 10, height = 10)

# Melt the correlation matrix
library(reshape2)
melted_cormat <- melt(upper_tri, na.rm = TRUE)

# Heatmap
library(ggplot2)
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
   midpoint = 0, limit = c(-1,1), space = "Lab", 
   name="Pearson\nCorrelation") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 12, hjust = 1))+
 coord_fixed()
dev.off()

# # Plotting general 2

reorder_cormat <- function(cormat){
# Use correlation between variables as distance
dd <- as.dist((1-cormat)/2)
hc <- hclust(dd)
cormat <-cormat[hc$order, hc$order]
}
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
   midpoint = 0, limit = c(-1,1), space = "Lab", 
    name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 12, hjust = 1))+
 coord_fixed()


pdf("test3.pdf", width = 10, height = 10)
ggheatmap + 
geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank(),
  legend.justification = c(1, 0),
  legend.position = c(0.6, 0.7),
  legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                title.position = "top", title.hjust = 0.5))

dev.off()

library(heatmaply)
heatmaply_cor(x = cormat, xlab = "Features",ylab = "Features", k_col = 2, k_row = 2)


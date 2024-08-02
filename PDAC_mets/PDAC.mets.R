#Created : Ateeq
#on:04-13-23
#Project: PDAC_METS Spatial and our pdac samples analysis
####
# part1
library(Seurat)

# Reading all the 68 files from tempus

# Set the directory where your data is stored
data_dir <- "/Users/akhaliq/Desktop/spatial_analysis/ILR/samples"

# Get a list of folders in the directory that end with "b"
folders <- list.files(data_dir, full.names = TRUE)

# Create an empty list to store the Seurat objects
seurat_objects <- list()

# Loop over the folders and read in the data for each one
for (folder in folders) {
  cat("Processing folder:", folder, "\n")
  #/Users/akhaliq/Desktop/spatial_analysis/ILR/samples/PDACLN_5/outs/filtered_feature_bc_matrix.h5
  # Read in the filtered_feature_bc_matrix.h5 file
  expr.url <- file.path(folder, "/outs/filtered_feature_bc_matrix.h5")
  expr.data <- Read10X_h5(filename =  expr.url)
  
  # Create a Seurat object
  seurat_obj <- CreateSeuratObject(counts = expr.data, project = basename(folder), 
                                    assay = 'Spatial', min.cells = 3, min.features = 200)
  
  # Read in the image data
  img.dir <- file.path(folder, "/outs/spatial")
  seurat_img <- Read10X_Image(image.dir = img.dir)
  DefaultAssay(object = seurat_img) <- 'Spatial'
  
  # Subset the image data to the cells in the Seurat object
  seurat_img <- seurat_img[colnames(seurat_obj)]
  seurat_obj[[basename(folder)]] <- seurat_img
  
  # Store the Seurat object in the list, with the name of the sample
  seurat_objects[[basename(folder)]] <- seurat_obj
}

# Save each Seurat object to a separate file, using the sample name as the filename
for (sample_name in names(seurat_objects)) {
  seurat_obj <- seurat_objects[[sample_name]]
  output_file <- file.path(data_dir, paste0(sample_name, ".rds"))
  saveRDS(seurat_obj, file = output_file)
}

# part 3
# reading .rds (68 TEMPUS files and 13 Pancreas our files = 81 samples)

# Set the directory where your data is stored
data_dir <- "/Users/akhaliq/Desktop/spatial_analysis/ILR/samples/"

# Get a list of .rds files in the directory
rds_files <- list.files(path = data_dir, pattern = "\\.rds$", full.names = TRUE)

# Loop over the .rds files and read them into separate R objects
for (rds_file in rds_files) {
  # Construct the name for the R object (use the file name without the ".rds" extension)
  r_obj_name <- basename(rds_file)
  r_obj_name <- substr(r_obj_name, 1, nchar(r_obj_name)-4)
  
  # Read the R object from the .rds file and assign it to a variable with the same name as the object
  assign(r_obj_name, readRDS(rds_file))
}

# part 4
#Merging as one object 

pdac_all <- merge(IU_PDA_HM10, c(IU_PDA_HM11, IU_PDA_HM12, IU_PDA_HM13, IU_PDA_HM14, IU_PDA_HM3, IU_PDA_HM4, IU_PDA_HM5, IU_PDA_HM6, IU_PDA_HM7, IU_PDA_HM9,IU_PDA_LNM11, IU_PDA_LNM13, IU_PDA_LNM7, IU_PDA_LNM8, IU_PDA_LNM9, IU_PDA_NH3, IU_PDA_NP11, IU_PDA_NP12, IU_PDA_NP2, IU_PDA_T1, IU_PDA_T10, IU_PDA_T11, IU_PDA_T12, IU_PDA_T13, IU_PDA_T2, IU_PDA_T4, IU_PDA_T5, IU_PDA_T7, IU_PDA_T9))


# part 5
# remove all Except pdac_all

rm(list=setdiff(ls(),"pdac_all"))
gc()


# part 7

# showed an error "Error in FUN(left, right) : non-numeric argument to binary operator"
# The following code should loop through all image names in the pdac_all Seurat object and apply the as.integer() function to the specified fields for each image.


image_names <- c("IU_PDA_HM10","IU_PDA_HM11","IU_PDA_HM12","IU_PDA_HM13","IU_PDA_HM14","IU_PDA_HM3","IU_PDA_HM4","IU_PDA_HM5","IU_PDA_HM6","IU_PDA_HM7","IU_PDA_HM9","IU_PDA_LNM11","IU_PDA_LNM13","IU_PDA_LNM7","IU_PDA_LNM8","IU_PDA_LNM9","IU_PDA_NH3","IU_PDA_NP11","IU_PDA_NP12","IU_PDA_NP2","IU_PDA_T1","IU_PDA_T10","IU_PDA_T11","IU_PDA_T12","IU_PDA_T13","IU_PDA_T2","IU_PDA_T4","IU_PDA_T5","IU_PDA_T7","IU_PDA_T9")

for (i in image_names) {
  pdac_all@images[[i]]@coordinates$tissue <- as.integer(pdac_all@images[[i]]@coordinates$tissue)
  pdac_all@images[[i]]@coordinates$row <- as.integer(pdac_all@images[[i]]@coordinates$row)
  pdac_all@images[[i]]@coordinates$col <- as.integer(pdac_all@images[[i]]@coordinates$col)
  pdac_all@images[[i]]@coordinates$imagerow <- as.integer(pdac_all@images[[i]]@coordinates$imagerow)
  pdac_all@images[[i]]@coordinates$imagecol <- as.integer(pdac_all@images[[i]]@coordinates$imagecol)
}

#In this code, we first define a vector image_names that contains all the names of the images we want to modify. Then, we use a for loop to iterate through each image name in pdac_all For each image, we use double brackets to access the image object within the pdac_all@images list, and then apply the code to modify its coordinates.

png("2_Before_Qc_Images.png", units="in", width=100, height=30, res=300)
SpatialFeaturePlot(pdac_all, features = c("nCount_Spatial"))
dev.off()

# part 8
#Check for mito and ribo contamination

pdac_all <- PercentageFeatureSet(pdac_all, "^mt-", col.name = "percent_mito")
pdac_all <- PercentageFeatureSet(pdac_all, "^Hb.*-", col.name = "percent_hb")

# no mito or ribo is seen

#


# Filter Bl1
pdac_all <- pdac_all[!grepl("Bc1", rownames(pdac_all)), ]

# Filter Mitocondrial
pdac_all <- pdac_all[!grepl("^mt-", rownames(pdac_all)), ]

# Filter Hemoglobin gene (optional if that is a problem on your data)
pdac_all <- pdac_all[!grepl("^Hb.*-", rownames(pdac_all)), ]

dim(pdac_all)

saveRDS(pdac_all,"pdac_all.rds")

pdac_all <- SCTransform(pdac_all, assay = "Spatial", verbose = TRUE, method = "poisson")
pdac_all <- RunPCA(pdac_all, assay = "SCT", verbose = FALSE)
pdac_all <- FindNeighbors(pdac_all, reduction = "pca", dims = 1:30)
pdac_all <- FindClusters(pdac_all, verbose = FALSE)
pdac_all <- RunUMAP(pdac_all, reduction = "pca", dims = 1:30)

saveRDS(pdac_all,"pdac_all.rds")

png("3_umap_beforeQC_Images.png", units="in", width=100, height=30, res=300)
DimPlot(pdac_all, reduction = "umap", group.by = c("ident", "orig.ident"))
dev.off()

# create a list of the original data that we loaded to start with
st.list <- SplitObject(pdac_all, split.by = "orig.ident")

# run SCT on both datasets
st.list = lapply(st.list, SCTransform, assay = "Spatial", method = "poisson")


# need to set maxSize for PrepSCTIntegration to work
options(future.globals.maxSize = 20000000000000 * 1024^2)  # set allowed size to 2K MiB


st.features = SelectIntegrationFeatures(st.list, nfeatures = 3000, verbose = FALSE)
st.list <- PrepSCTIntegration(object.list = st.list, anchor.features = st.features,
    verbose = FALSE)
int.anchors <- FindIntegrationAnchors(object.list = st.list, normalization.method = "SCT",
    verbose = FALSE, anchor.features = st.features)
pdac.integrated <- IntegrateData(anchorset = int.anchors, normalization.method = "SCT",
    verbose = FALSE)

rm(int.anchors, st.list)
gc()

pdac.integrated <- RunPCA(pdac.integrated, verbose = FALSE)
pdac.integrated <- FindNeighbors(pdac.integrated, dims = 1:30)
pdac.integrated <- FindClusters(pdac.integrated, verbose = FALSE)
pdac.integrated <- RunUMAP(pdac.integrated, dims = 1:30)
pdac.integrated <- RunTSNE(pdac.integrated, dims = 1:30)

png("4_umap_aftrerQC_Images.png", units="in", width=100, height=30, res=300)
DimPlot(pdac.integrated, reduction = "umap", group.by = c("ident", "orig.ident"))
dev.off()

pdac$seurat_clusters <- as.factor(as.numeric(as.character(pdac$seurat_clusters)) + 1)

saveRDS(pdac.integrated,"pdac_all.rds")

# There are 900 Images i need to keep only 30

pdac <- readRDS("pdac_all.rds")

SpatialDimPlot(pdac,group.by = "res_0.8",images=c("PDACP_1","PDACP_3.3"))



# to check the image Coordinate present or absent 

# Create a vector of image names
#image_names <- c("PDACNP_1.1", "PDACNP_10.1", "PDACNP_2.1", "PDACP_1.1", "PDACP_10.1", "PDACP_11.1", "PDACP_12.1", "PDACP_13.1", "PDACP_3.1", "PDACP_4.1", "PDACP_5.1", "PDACP_8.1", "PDACP_9.1", "MV001b.2")
# Get the names of the elements in pdac@images as a character vector
pdac_images_names <- names(pdac@images)

# Create a vector with the image names
image_names <- as.vector(pdac_images_names)

# Create a data frame to store image names and their presence or absence of coordinates
image_coordinates_df <- data.frame(Image_Name = character(0), Coordinates_Present = character(0))

# Iterate through the image names
for (image_name in image_names) {
  # Get the coordinates using GetTissueCoordinates function
  coordinates <- GetTissueCoordinates(pdac, image = image_name)
  
  # Check if the coordinates data frame is empty (i.e., not NULL)
  if (!is.null(coordinates) && nrow(coordinates) > 0) {
    # If not NULL and not empty, add the image name and "Coordinates Present" to image_coordinates_df
    image_coordinates_df <- rbind(image_coordinates_df, cbind(Image_Name = image_name, Coordinates_Present = "Coordinates Present"))
  } else {
    # If NULL or empty, add the image name and "Coordinates Absent" to image_coordinates_df
    image_coordinates_df <- rbind(image_coordinates_df, cbind(Image_Name = image_name, Coordinates_Present = "Coordinates Absent"))
  }
}

# Print the data frame with image names and whether coordinates are present or absent
print(head(image_coordinates_df))

# Pull the image names where coordinates are absent
image_names_absent <- image_coordinates_df$Image_Name[image_coordinates_df$Coordinates_Present == "Coordinates Absent"]
image_names_present <- image_coordinates_df$Image_Name[image_coordinates_df$Coordinates_Present == "Coordinates Present"]


length(image_names_present)
length(image_names_absent)

#######

retain_names = image_names_present

# List of image names to retain
#retain_names <- c("IU_PDA_HM10", "IU_PDA_HM11.1", "IU_PDA_HM12.2", "IU_PDA_HM13.3", "IU_PDA_HM2.4", "IU_PDA_HM3.5", "IU_PDA_HM4.6", "IU_PDA_HM5.7", "IU_PDA_HM6.8", "IU_PDA_HM8.9", "IU_PDA_HM9.10", "IU_PDA_LNM10.11", "IU_PDA_LNM12.12", "IU_PDA_LNM6.13", "IU_PDA_LNM7.14", "IU_PDA_LNM8.15", "IU_PDA_NH2.16", "IU_PDA_NP10.17", "IU_PDA_NP11.18", "IU_PDA_NP2.19", "IU_PDA_T1.20", "IU_PDA_T10.21", "IU_PDA_T11.22", "IU_PDA_T12.23", "IU_PDA_T2.24", "IU_PDA_T3.25", "IU_PDA_T4.26", "IU_PDA_T6.27", "IU_PDA_T8.28", "IU_PDA_T9.29")

# Loop through all image names in the pdac object and delete those not in the list of retained names
for (name in names(pdac@images)) {
  if (!(name %in% retain_names)) {
    eval(parse(text=paste0("pdac@images$", name, " <- NULL")))
  }
}

# since the new integrated data saves the images names with "." lets delete it and make it look good.
# Create vector of new names
new_names <- gsub("\\..*", "", names(pdac@images))

# Rename images
names(pdac@images) <- new_names

#
Idents(myeloid) <- "seurat_clusters"



###
********************************************
#Running ILR
********************************************

#modify the integrated compositions 
integrated_compositions <- readRDS("/Users/akhaliq/Desktop/spatial_analysis/ILR/integrated_compositions_new.rds")
rownames(integrated_compositions) <- make.unique(rownames(integrated_compositions))
rownames(integrated_compositions) <- gsub("\\.", "_", rownames(integrated_compositions))

# remove the string from column names
colnames(integrated_compositions) <- gsub("q05cell_abundance_w_sf_", "", colnames(integrated_compositions))

print("1-Loading the Library")
library(Seurat)
library(STdeconvolve)
library(compositions)
library(tidyverse)
library(clustree)
library(uwot)
library(scran)
library(cluster)
print("1-Done")

print("2-Loading the Rdata")
#load("/N/project/akhaliq/spatial_samples/st_analysis/monad/STdeconvolve.pdac_integ18.RData")
print("2-Done")

print("3-Starting the Actual Work NOW")
#rownames(deconProp)<-gsub(pattern = ".", replacement = "-",x = rownames(deconProp), fixed = TRUE)
#std.meta <- merge(deconProp,pdac.integrated@meta.data,by=0)
#head(std.meta)

#integrated_compositions = deconProp
#integrated_compositions <- readRDS("integrated_compositions.rds")

#integrated_compositions = read.csv("/Users/akhaliq/Desktop/spatial_analysis/new_monads/merged.csv",header=T,sep=",",row.names=1)
#colnames(integrated_compositions) <- sub("q05cell_abundance_w_sf_", "", colnames(integrated_compositions))
baseILR <- ilrBase(x = integrated_compositions, method = "basic")
head(baseILR)
cell_ilr <- as.matrix(ilr(integrated_compositions, baseILR))
colnames(cell_ilr) <- paste0("ILR_", 1:ncol(cell_ilr))
k_vect <- c(10,15,25)
k_vect <- set_names(k_vect, paste0("k_",k_vect))
head(cell_ilr)
dim(integrated_compositions)
head(integrated_compositions)
baseILR
# Make community graph
#k_vect <- c(10,15,25,50,75,100)
k_vect <- set_names(k_vect, paste0("k_",k_vect))
cluster_info <- map(k_vect, function(k) {
print(k)
print("4-Generating SNN")
snn_graph <- scran::buildSNNGraph(x = t(cell_ilr %>% as.data.frame() %>% as.matrix()), k = k)
print("Louvain clustering")
clust.louvain <- igraph::cluster_louvain(snn_graph)
clusters <- tibble(cluster = clust.louvain$membership,
spot_id = rownames(cell_ilr))
})
print("4-Done")

print("5-saving the RDS")
#saveRDS(snn_graph,"snn_graph.rds")
print("5-Done")

print(head(cluster_info))

cluster_info <- cluster_info %>%
enframe() %>%
unnest() %>%
pivot_wider(names_from = name,
values_from = cluster)
# For each k, make a 60% subsampling -----------------------------------------------------

set.seed(241099)
k_vect <- set_names(names(k_vect))

subsampling_mean_ss <- map(k_vect, function(k) {
#print(k)
cluster_info_summary <- cluster_info %>%
group_by_at(k) %>%
summarize(ncells = floor(n() * 0.3))
#print(cluster_info_summary)
cells <- cluster_info %>%
dplyr::select_at(c("spot_id", k)) %>%
group_by_at(k) %>%
nest() %>%
left_join(cluster_info_summary) %>%
mutate(data = map(data, ~ .x[[1]])) %>%
mutate(selected_cells = map2(data, ncells, function(dat,n) {
sample(dat, n)
})) %>%
pull(selected_cells) %>%
unlist()
#print(cells)
dist_mat <- dist(cell_ilr[cells, ])
#print(dist_mat)
k_vect <- purrr::set_names(cluster_info[[k]], cluster_info[["spot_id"]])[cells]
#print(k_vect)
sil <- cluster::silhouette(x = k_vect, dist = dist_mat)
#saveRDS(sil,paste(k,"_sil.rds"))
#png(paste(k,".png"))
#plot(sil)
#dev.off()
mean(sil[, 'sil_width'])
print(mean(sil[, 'sil_width']))
})

subsampling_mean_ss <- enframe(subsampling_mean_ss) %>%
unnest() %>%
dplyr::filter()
monad_resolution <- dplyr::filter(subsampling_mean_ss,
value == max(value)) %>%
pull(name)
head(monad_resolution)

# Create UMAP and plot the compositions

comp_umap <- umap(as.matrix(cell_ilr),
n_neighbors = 30,
n_epochs = 1000,
metric = "cosine") %>%
as.data.frame() %>%
mutate(row_id = rownames(cell_ilr))
comp_umap <- comp_umap[, 1:3]
comp_umap <- comp_umap %>%
left_join(cluster_info,
by = c("row_id" = "spot_id"))
#write_csv(comp_umap, file = "./results/monad_mapping/ct_monads/umap_compositional.csv")

write_csv(comp_umap, file = "umap_compositional.csv")
#pdf("./results/monad_mapping/ct_monads/ct_ILR_umap.pdf", height = 6, width = 7)

pdf("ct_ILR_umap.pdf", height = 6, width = 7)
log_comps <- log10(integrated_compositions)
plt <- comp_umap %>%
ggplot(aes(x = V1, y = V2,
color = as.character(monad_resolution))) +
ggrastr::geom_point_rast(size = 0.1) +
theme_classic() +
xlab("UMAP1") +
ylab("UMAP2") +
guides(colour = guide_legend(override.aes = list(size=4)))
plot(plt)
plt
dev.off()
#dev.off()

comp_umap %>%
ggplot(aes(x = V1, y = V2,
color = as.character(monad_resolution))) +
ggrastr::geom_point_rast(size = 0.1) +
theme_classic() +
xlab("UMAP1") +
ylab("UMAP2") +
guides(colour = guide_legend(override.aes = list(size=4)))
dev.off()
 #pdf("ct_ILR_umap.pdf", height = 6, width = 7)
 pdf("ct_ILR_umap.pdf", height = 6, width = 7)
comp_umap %>%
ggplot(aes(x = V1, y = V2,
color = as.character(monad_resolution))) +
ggrastr::geom_point_rast(size = 0.1) +
theme_classic() +
xlab("UMAP1") +
ylab("UMAP2") +
guides(colour = guide_legend(override.aes = list(size=4)))
dev.off()

cts <- set_names(colnames(integrated_compositions))
walk(cts, function(ct){
plot_df <- comp_umap %>%
mutate(ct_prop = log_comps[ , ct])
plt <- plot_df %>%
ggplot(aes(x = V1, y = V2,
color = ct_prop)) +
ggrastr::geom_point_rast(size = 0.07) +
theme_classic() +
ggtitle(ct) +
xlab("UMAP1") +
ylab("UMAP2")
plot(plt)
})

dev.off()
cluster_info <- comp_umap %>%
dplyr::select(c("row_id",monad_resolution)) %>%
dplyr::rename("monad" = monad_resolution) %>%
dplyr::mutate(ct_monad = paste0("monad_", monad))
head(cluster_info)
monad_summary_pat <- integrated_compositions %>%
as.data.frame() %>%
rownames_to_column("row_id") %>%
pivot_longer(-row_id,values_to = "ct_prop",
names_to = "cell_type") %>%
left_join(cluster_info) %>%
mutate(orig.ident = strsplit(row_id, "[..]") %>%
map_chr(., ~ .x[1])) %>%
group_by(orig.ident, ct_monad, cell_type) %>%
summarize(median_ct_prop = median(ct_prop))
monad_summary <- monad_summary_pat %>%
ungroup() %>%
group_by(ct_monad, cell_type) %>%
summarise(patient_median_ct_prop = median(median_ct_prop))
head(monad_summary)
monad_summary_mat <- monad_summary %>%
pivot_wider(values_from = patient_median_ct_prop,
names_from = cell_type, values_fill = 0) %>%
column_to_rownames("ct_monad") %>%
as.matrix()
monad_order <- hclust(dist(monad_summary_mat))
monad_order <- monad_order$labels[monad_order$order]
ct_order <- hclust(dist(t(monad_summary_mat)))
ct_order <- ct_order$labels[ct_order$order]
# Find characteristic cell types of each monad
# We have per patient the proportion of each cell-type in each monad
run_wilcox_up <- function(prop_data) {
prop_data_group <- prop_data[["ct_monad"]] %>%
unique() %>%
set_names()
map(prop_data_group, function(g) {
test_data <- prop_data %>%
mutate(test_group = ifelse(ct_monad == g,
"target", "rest")) %>%
mutate(test_group = factor(test_group,
levels = c("target", "rest")))
wilcox.test(median_ct_prop ~ test_group,
data = test_data,
alternative = "greater") %>%
broom::tidy()
}) %>% enframe("ct_monad") %>%
unnest()
}
wilcoxon_res <- monad_summary_pat %>%
ungroup() %>%
group_by(cell_type) %>%
nest() %>%
mutate(wres = map(data, run_wilcox_up)) %>%
dplyr::select(wres) %>%
unnest() %>%
ungroup() %>%
mutate(p_corr = p.adjust(p.value))
wilcoxon_res <- wilcoxon_res %>%
mutate(significant = ifelse(p_corr <= 0.15, "*", ""))
write.table(monad_summary_pat, file = "monad_summary_pat.txt",col.names = T, row.names = F, quote = F, sep = "\t")
write.table(wilcoxon_res, file = "wilcoxon_res_cells_monads.txt",
col.names = T, row.names = F, quote = F, sep = "\t")
mean_ct_prop_plt <- monad_summary %>%
left_join(wilcoxon_res, by = c("ct_monad", "cell_type")) %>%
mutate(cell_type = factor(cell_type, levels = ct_order),
ct_monad = factor(ct_monad, levels = monad_order)) %>%
ungroup() %>%
group_by(cell_type) %>%
mutate(scaled_pat_median = (patient_median_ct_prop - mean(patient_median_ct_prop))/sd(patient_median_ct_prop)) %>%
ungroup() %>%
ggplot(aes(x = cell_type, y = ct_monad, fill = scaled_pat_median)) +
geom_tile() +
geom_text(aes(label = significant)) +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
legend.position = "bottom",
plot.margin = unit(c(0, 0, 0, 0), "cm"),
axis.text.y = element_text(size=12)) +
scale_fill_gradient(high = "#9b59b6", low = "#f9cb9c")  ##9b59b6,  #f9cb9c
cluster_counts <- cluster_info %>%
dplyr::select_at(c("row_id", "ct_monad")) %>%
group_by(ct_monad) %>%
summarise(nspots = length(ct_monad)) %>%
mutate(prop_spots = nspots/sum(nspots))
write_csv(cluster_counts, file = "monad_prop_summary.csv")
barplts <- cluster_counts %>%
mutate(ct_monad = factor(ct_monad, levels = monad_order)) %>%
ggplot(aes(y = ct_monad, x = prop_spots)) +
geom_bar(stat = "identity") +
theme_classic() + ylab("") +
theme(axis.text.y = element_blank(),
plot.margin = unit(c(0, 0, 0, 0), "cm"),
axis.text.x = element_text(size=12))
monad_summary_plt <- cowplot::plot_grid(mean_ct_prop_plt, barplts, align = "hv", axis = "tb")
pdf("characteristic_ct_monads.pdf", height = 10, width = 20)
plot(monad_summary_plt)
dev.off()
#dev.off()
#monad_summary_plt
plt <- monad_summary_pat %>%
ggplot(aes(x = ct_monad, y = median_ct_prop)) +
geom_boxplot() +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
facet_wrap(.~cell_type, ncol = 3,scales = "free_y")

pdf("box_plot_monads.pdf",height = 10, width = 20)
plt
dev.off()

print("monad Completed !!!")

savehistory("monad_script.RHistory")
save.image("monad_k25.RData")


********************************************
# END of  ILR
********************************************
#change barcodes in Seurat obj


#Since the cell2Loc barcodes are different than the seurat barcodes lets change the barcodes in seurat object
write.csv(pdac@meta.data,"pdac.meta.csv")
#then delete all except the barcodes and create two colnames one with "old" and another with "new", then modify the barcodes as the requirment, save the csv and read it
barcode_csv <- read.csv("/Users/akhaliq/Desktop/spatial_analysis/ILR/pdac.meta.csv")
#Create a named vector from the CSV file, where the old barcode is the name of the vector element, and the new barcode is the value of the vector element. You can use the setNames() function to create the named vector.

barcode_vector <- setNames(barcode_csv$new, barcode_csv$old)

#Use the RenameCells() function from the Seurat package to update the barcode information in the Seurat object. Pass the named vector created in step 2 to the RenameCells() function.

pdac <- RenameCells(pdac, new.names = barcode_vector)

# check if this is fine
head(pdac)
# This should update the barcode information in the Seurat object with the new barcode information from the CSV file. Note that the new.names argument in the RenameCells() function expects a named vector where the names correspond to the old barcode names and the values correspond to the new barcode names.

saveRDS(pdac,"/Users/akhaliq/Desktop/spatial_analysis/new_object/pdac_cl_renamed_mod.rds")

###
#Subset only those barcodes present in both seurat and cell2loc.

library(data.table)
umi_list <- fread("/Users/akhaliq/Desktop/spatial_analysis/ILR/umap_compositional.csv")
umi_vec <- umi_list$row_id
pdac_sub <- subset(pdac, cells = umi_vec)


#####
## adding the UMAP and TSNE to the seurat obj

write.csv(emb,"tsne_stdcon.csv")
write.csv(emb.umap,"umap_stdcon.csv")

emb.tsne <- read.csv("tsne_stdcon.csv",header=T,sep=",",row.names=1)
emb.umap1<- read.csv("umap_stdcon.csv",header=T,sep=",",row.names=1)

colnames(emb.tsne)<- c("tSNE1","tSNE2")
colnames(emb.umap1)<- c("UMAP1","UMAP2")

pdac[['tsne.std']] <- CreateDimReducObject(embeddings = as.matrix(emb.tsne), key = 'tsne_', assay = 'integrated')
pdac[['umap.std']] <- CreateDimReducObject(embeddings = as.matrix(emb.umap1), key = 'umap_', assay = 'integrated')


library(dittoSeq)
#add clusters
std_clusters <- data.frame(tempCom)
colnames(std_clusters) <- "std.clusters"
pdac <- AddMetaData(pdac,std_clusters)



png("monads_pdac.png", units="in", width=10, height=10, res=300)
dittoDimPlot(pdac_sub, reduction = "umap.c2l", "monads")
dev.off()

png("stdcon_tsne_cluster1.png", units="in", width=10, height=10, res=300)
dittoDimPlot(pdac, reduction = "tsne.std", "std.clusters")
dev.off()


png("stdcon_umap_cluster1.png", units="in", width=10, height=10, res=300)
dittoDimPlot(pdac, reduction = "umap.std", "std.clusters")
dev.off()

save.image("stdcon.recent.RData")


*******
# bar plots for each monad

library(ggplot2)
library(gridExtra)

# Read the data from the CSV file
data <- read.csv("monad_summary.csv")

# Create a list to store the plots
plots <- list()

# Create a barplot for each ct_monad
for (monad in unique(data$ct_monad)) {
  monad_data <- subset(data, ct_monad == monad)
  plot <- ggplot(monad_data, aes(x = cell_type, y = patient_median_ct_prop, fill = cell_type)) +
    geom_bar(stat = "identity") +
    ggtitle(paste("Monad", monad)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    xlab("Cell Type") +
    ylab("Patient Median CT Proportion")
  plots[[monad]] <- plot
}

# Combine all the plots into one image
combined_plot <- do.call(grid.arrange, c(plots, ncol = 3))

# Save the image to a file
ggsave("monad_plots.png", combined_plot, width = 12, height = 8)


###
library(ggplot2)
library(gridExtra)

# Read the data
monad_summary <- read.csv("monad_summary.csv", header = TRUE, sep = "\t")

# Define a larger color palette with 25 distinct colors
my_colors <- rainbow(length(unique(monad_summary$monads)), s = 0.8, v = 0.8)

# Create a list of ggplot objects for each monad
monad_plots <- lapply(unique(monad_summary$monads), function(m) {
  data <- subset(monad_summary, monads == m)
  ggplot(data, aes(x=cell_type, y=median_ct_prop, fill=monads)) + 
    geom_bar(stat="identity") +
    ggtitle(paste("Cell type proportions for", m)) +
    scale_fill_manual(values = my_colors) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
})

# Combine all plots into a single PNG file
png("monad_summary.png", width = 1600, height = 1000)
grid.arrange(grobs = monad_plots, ncol=5)
dev.off()


# box plots

library(ggplot2)
library(dplyr)

# plot
ggplot(df, aes(x = cell_type, y = median_ct_prop, fill = monads)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Monad", y = "Median CT Proportion", fill = "Monads") +
  guides(fill = FALSE) # hide legend


# Bar Plots
png("barplots_monads_pdac_origin.png", units="in", width=10, height=10, res=300)
dittoBarPlot(pdac_sub, "Origin", group.by = "monads")
dev.off()


png("barplots_monads_pdac_histology.png", units="in", width=10, height=10, res=300)
dittoBarPlot(pdac_sub, "Histology", group.by = "monads")
dev.off()

png("umapc2l_monads_pdac.png", units="in", width=10, height=10, res=300)
DimPlot(pdac_sub,reduction="umap.c2l",group.by="monads",label=T)
dev.off()

png("umap_seuratclusters_pdac.png", units="in", width=10, height=10, res=300)
DimPlot(pdac_sub,reduction="umap",group.by="seurat_clusters",label=T)
dev.off()


#####


plot1 <- VlnPlot(p2, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(p2, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

p2 <- SCTransform(p2, assay = "Spatial", verbose = FALSE)

p2 <- RunPCA(p2, assay = "SCT", verbose = FALSE)
p2 <- FindNeighbors(p2, reduction = "pca", dims = 1:30)
p2 <- FindClusters(p2, verbose = FALSE)
p2 <- RunUMAP(p2, reduction = "pca", dims = 1:30)

p1 <- DimPlot(p2, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(p2, label = TRUE, label.size = 3)
p1 + p2


###
pdac_allres <- FindClusters(pdac, resolution = c(0,0.2,0.25,0.3,0.35,0.4,0.5,0.6,0.7,0.8) , algorithm = 1, graph.name= "integrated_snn")

pdac_allres$integrated_snn_res.0.2  pdac_allres$integrated_snn_res.0.4  pdac_allres$integrated_snn_res.0.6  pdac_allres$integrated_snn_res.0.8  
pdac_allres$integrated_snn_res.0.3  pdac_allres$integrated_snn_res.0.5  pdac_allres$integrated_snn_res.0.7  

pdac_allres$integrated_snn_res.0.2 <- as.factor(as.numeric(as.character(pdac_allres$integrated_snn_res.0.2)) + 1)
pdac_allres$integrated_snn_res.0.25 <- as.factor(as.numeric(as.character(pdac_allres$integrated_snn_res.0.25)) + 1)
pdac_allres$integrated_snn_res.0.3 <- as.factor(as.numeric(as.character(pdac_allres$integrated_snn_res.0.3)) + 1)
pdac_allres$integrated_snn_res.0.35 <- as.factor(as.numeric(as.character(pdac_allres$integrated_snn_res.0.35)) + 1)
pdac_allres$integrated_snn_res.0.4 <- as.factor(as.numeric(as.character(pdac_allres$integrated_snn_res.0.4)) + 1)
pdac_allres$integrated_snn_res.0.5 <- as.factor(as.numeric(as.character(pdac_allres$integrated_snn_res.0.5)) + 1)
pdac_allres$integrated_snn_res.0.6 <- as.factor(as.numeric(as.character(pdac_allres$integrated_snn_res.0.6)) + 1)
pdac_allres$integrated_snn_res.0.7 <- as.factor(as.numeric(as.character(pdac_allres$integrated_snn_res.0.7)) + 1)
pdac_allres$integrated_snn_res.0.8 <- as.factor(as.numeric(as.character(pdac_allres$integrated_snn_res.0.8)) + 1)

library(clustree)

clustree(pdac_allres@meta.data, prefix = "integrated_snn_res.")



png("res_dim.png", units="in", width=20, height=15, res=300)
cowplot::plot_grid(ncol = 3,
  DimPlot(pdac_allres, reduction = "umap", group.by = "integrated_snn_res.0.2",label=T,label.box=T,label.size = 6, repel = T)+ggtitle("res.0.2"),
  DimPlot(pdac_allres, reduction = "umap", group.by = "integrated_snn_res.0.3",label=T,label.box=T,label.size = 6, repel = T)+ggtitle("res.0.3"),
  DimPlot(pdac_allres, reduction = "umap", group.by = "integrated_snn_res.0.4",label=T,label.box=T,label.size = 6, repel = T)+ggtitle("res.0.4"),

  DimPlot(pdac_allres, reduction = "umap", group.by = "integrated_snn_res.0.5",label=T,label.box=T,label.size = 6, repel = T)+ggtitle("res.0.5"),
  DimPlot(pdac_allres, reduction = "umap", group.by = "integrated_snn_res.0.6",label=T,label.box=T,label.size = 6, repel = T)+ggtitle("res.0.6"),
  DimPlot(pdac_allres, reduction = "umap", group.by = "integrated_snn_res.0.7",label=T,label.box=T,label.size = 6, repel = T)+ggtitle("res.0.7"),

  DimPlot(pdac_allres, reduction = "umap", group.by = "integrated_snn_res.0.8",label=T,label.box=T,label.size = 6, repel = T)+ggtitle("res.0.8"),
  DimPlot(pdac_allres, reduction = "umap", group.by = "integrated_snn_res.0.25",label=T,label.box=T,label.size = 6, repel = T)+ggtitle("res.0.25"),
  DimPlot(pdac_allres, reduction = "umap", group.by = "integrated_snn_res.0.35",label=T,label.box=T,label.size = 6, repel = T)+ggtitle("res.0.35")
)
dev.off()



png("res_dim.tsne.png", units="in", width=20, height=15, res=300)
cowplot::plot_grid(ncol = 3,
  DimPlot(pdac_allres, reduction = "tsne", group.by = "integrated_snn_res.0.2",label=T,label.box=T,label.size = 6, repel = T)+ggtitle("res.0.2"),
  DimPlot(pdac_allres, reduction = "tsne", group.by = "integrated_snn_res.0.3",label=T,label.box=T,label.size = 6, repel = T)+ggtitle("res.0.3"),
  DimPlot(pdac_allres, reduction = "tsne", group.by = "integrated_snn_res.0.4",label=T,label.box=T,label.size = 6, repel = T)+ggtitle("res.0.4"),

  DimPlot(pdac_allres, reduction = "tsne", group.by = "integrated_snn_res.0.5",label=T,label.box=T,label.size = 6, repel = T)+ggtitle("res.0.5"),
  DimPlot(pdac_allres, reduction = "tsne", group.by = "integrated_snn_res.0.6",label=T,label.box=T,label.size = 6, repel = T)+ggtitle("res.0.6"),
  DimPlot(pdac_allres, reduction = "tsne", group.by = "integrated_snn_res.0.7",label=T,label.box=T,label.size = 6, repel = T)+ggtitle("res.0.7"),

  DimPlot(pdac_allres, reduction = "tsne", group.by = "integrated_snn_res.0.8",label=T,label.box=T,label.size = 6, repel = T)+ggtitle("res.0.8"),
  DimPlot(pdac_allres, reduction = "tsne", group.by = "integrated_snn_res.0.25",label=T,label.box=T,label.size = 6, repel = T)+ggtitle("res.0.25"),
  DimPlot(pdac_allres, reduction = "tsne", group.by = "integrated_snn_res.0.35",label=T,label.box=T,label.size = 6, repel = T)+ggtitle("res.0.35")
)
dev.off()



##Marker genes

Idents(pdac_allres) <- "integrated_snn_res.0.3"

table(pdac_allres$integrated_snn_res.0.25)
png("UMAP_tcells.0.8_clus24_rem.png", units="in", width=10, height=10, res=300)
DimPlot(pdac_allres, reduction = "umap", label=T, group.by="seurat_clusters")
dev.off()

markers_genes_pdac_0.3 <- FindAllMarkers(pdac_allres, logfc.threshold = 0.2, test.use = "wilcox", only.pos = TRUE,assay = "Spatial")
write.table(markers_genes_pdac_0.3, file="markers_genes_pdac_0.3.txt", sep="\t",append=F, row.names = F)
cbind(table(markers_genes_pdac_0.3$cluster))

 
Idents(pdac_allres) <- "integrated_snn_res.0.3"

markers_genes_pdac_0.3 <- FindAllMarkers(pdac_allres, logfc.threshold = 0.2, test.use = "wilcox", only.pos = TRUE,assay = "Spatial")
write.table(markers_genes_pdac_0.3, file="markers_genes_pdac_0.3.txt", sep="\t",append=F, row.names = F)
cbind(table(markers_genes_pdac_0.3$cluster))


Idents(pdac_allres) <- "integrated_snn_res.0.35"

markers_genes_pdac_0.35 <- FindAllMarkers(pdac_allres, logfc.threshold = 0.2, test.use = "wilcox", only.pos = TRUE,assay = "Spatial")
write.table(markers_genes_pdac_0.35, file="markers_genes_pdac_0.35.txt", sep="\t",append=F, row.names = F)
cbind(table(markers_genes_pdac_0.35$cluster))


Idents(pdac_allres) <- "integrated_snn_res.0.4"

markers_genes_pdac_0.4 <- FindAllMarkers(pdac_allres, logfc.threshold = 0.2, test.use = "wilcox", only.pos = TRUE,assay = "Spatial")
write.table(markers_genes_pdac_0.4, file="markers_genes_pdac_0.4.txt", sep="\t",append=F, row.names = F)
cbind(table(markers_genes_pdac_0.4$cluster))

 save.image("pdac_st_till_rece.RData")
 


suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(cowplot)
    library(ggplot2)
    library(pheatmap)
    library(rafalib)
})

 markers_genes_pdac_0.4 %>%
    group_by(cluster) %>%
    top_n(-25, p_val_adj) -> top25
top25


#png("2_Barplots_top_exp_genes_each_clusters_0.3.png", units="in", width=20, height=10, res=300)
mypar(3, 5, mar = c(4, 6, 3, 1))
for (i in unique(top25$cluster)) {
    barplot(sort(setNames(top25$avg_log2FC, top25$gene)[top25$cluster == i], F),
        horiz = T, las = 1, main = paste0(i, " vs. rest"), border = "red", yaxs = "i")
    abline(v = c(0, 0.25), lty = c(1, 2))
}
dev.off()

#heatmap

markers_genes_pdac_0.2 %>%
    group_by(cluster) %>%
    top_n(-5, p_val_adj) -> top5


# create a scale.data slot for the selected genes
pdac_allres <- ScaleData(pdac_allres, features = as.character(unique(top5$gene)), assay = "Spatial")

png("3_heatmaps_top_exp_genes_each_clusters.png", units="in", width=20, height=20, res=300)
DoHeatmap(pdac_allres, features = as.character(unique(top5$gene)), group.by = "integrated_snn_res.0.2",
    assay = "Spatial")
dev.off()

png("4_Dotplots_top_exp_genes_each_clusters.png", units="in", width=20, height=20, res=300)
DotPlot(pdac_allres, features = rev(as.character(unique(top5$gene))), group.by = sel.clust,
    assay = "RNA") + coord_flip()
dev.off()


png("pdac_allres_markergenes_compartment.png", units="in", width=40, height=10, res=300)
cowplot::plot_grid(ncol = 6,
FeaturePlot(pdac_allres,features = "KRT18" ,min.cutoff="q9", cols=c("lightgrey", "#02a7c0"), label=FALSE)+ labs(title = "Epithelial Cells"),
FeaturePlot(pdac_allres,features = "CD3D" ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE)+ labs(title = "T Cells"),
FeaturePlot(pdac_allres,features = "CD79A" ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE)+ labs(title = "B Cells"),
FeaturePlot(pdac_allres,features = "LYZ" ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE)+ labs(title = "Myeloid"),
FeaturePlot(pdac_allres,features = "COL1A1" ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE)+ labs(title = "Fibroblast"),
FeaturePlot(pdac_allres,features = "CLDN5" ,min.cutoff="q9", cols=c("lightgrey", "red"), label=FALSE)+ labs(title = "Endothelial"))
dev.off()

dittoBarPlot(pdac_allres, "integrated_snn_res.0.2", group.by = "Histology")
dittoBarPlot(pdac_allres, "Histology", group.by = "integrated_snn_res.0.2")


###

# Spatial Dimension Plot Generation

This script generates spatial dimension plots for a list of images using the `SpatialDimPlot` function. Each plot is saved in a single PDF file with the corresponding image name as the title.


#code adds the image name as the title of each plot using the ggtitle() function. Each plot will now have the corresponding image name as its title in the PDF file. Please ensure you have the necessary data object (pdac_allres) and packages loaded before running the code.

# List of image names
image_names <- c("IU_PDA_T1", "IU_PDA_NP2", "IU_PDA_T2", "IU_PDA_NH2", "IU_PDA_HM2", "IU_PDA_HM3", "IU_PDA_T3", "IU_PDA_HM4", "IU_PDA_T4", "IU_PDA_HM5", "IU_PDA_HM6", "IU_PDA_LNM6", "IU_PDA_T6", "IU_PDA_LNM7", "IU_PDA_HM8", "IU_PDA_LNM8", "IU_PDA_T8", "IU_PDA_T9", "IU_PDA_HM9", "IU_PDA_HM10", "IU_PDA_LNM10", "IU_PDA_NP10", "IU_PDA_T10", "IU_PDA_HM11", "IU_PDA_T11", "IU_PDA_NP11", "IU_PDA_HM12", "IU_PDA_LNM12", "IU_PDA_T12", "IU_PDA_HM13")

# Open the PDF device
pdf("spatial_dim_plots.pdf")

# Loop over each image
for (image_name in image_names) {
  command <- paste0("SpatialDimPlot(pdac_allres, group.by = 'integrated_snn_res.0.4', images = '", image_name, "', facet.highlight = F, label = TRUE, label.size = 4, repel = T)")
  plot <- eval(parse(text = command))
  
  # Set the image name as the plot title
  plot <- plot + ggtitle(image_name)
  
  # Save the plot to the PDF file
  print(plot)
}

# Close the PDF device
dev.off()

########

# Changing the name of the images in St object 


# Define the old and new image name mapping
image_mapping <- c(
"IU_PDA_NP11" = "IU_PDA_NP11_t",
"IU_PDA_NP2" = "IU_PDA_NP2_t",
"IU_PDA_T11" = "IU_PDA_T11_t",
"IU_PDA_T3" = "IU_PDA_T3_t",
"IU_PDA_T8" = "IU_PDA_T8_t",
"IU_PDA_HM13" = "IU_PDA_HM13_t",
"IU_PDA_HM2" = "IU_PDA_HM2_t",
"IU_PDA_HM4" = "IU_PDA_HM4_t",
"IU_PDA_HM6" = "IU_PDA_HM6_t",
"IU_PDA_HM8" = "IU_PDA_HM8_t",
"IU_PDA_LNM7" = "IU_PDA_LNM7_t",
"IU_PDA_T1" = "IU_PDA_T1_t",
"IU_PDA_T4" = "IU_PDA_T4_t",
"IU_PDA_LNM10" = "IU_PDA_LNM10_t",
"IU_PDA_LNM6" = "IU_PDA_LNM6_t",
"IU_PDA_LNM8" = "IU_PDA_LNM8_t",
"IU_PDA_HM12" = "IU_PDA_HM12_t",
"IU_PDA_T9" = "IU_PDA_T9_t",
"IU_PDA_T2" = "IU_PDA_T2_t",
"IU_PDA_HM10" = "IU_PDA_HM10_t",
"IU_PDA_HM9" = "IU_PDA_HM9_t",
"IU_PDA_NP10" = "IU_PDA_NP10_t",
"IU_PDA_T10" = "IU_PDA_T10_t",
"IU_PDA_T12" = "IU_PDA_T12_t",
"IU_PDA_T6" = "IU_PDA_T6_t",
"IU_PDA_HM11" = "IU_PDA_HM11_t",
"IU_PDA_HM3" = "IU_PDA_HM3_t",
"IU_PDA_HM5" = "IU_PDA_HM5_t",
"IU_PDA_LNM12" = "IU_PDA_LNM12_t",
"IU_PDA_HM2_2" = "IU_PDA_HM2_2_t")

# Rename the image names in pdac_allres
names(pdac@images) <- image_mapping[names(pdac@images)]

# now remove the _t from the image to the following

# Define the old and new image name mapping
image_mapping <- c(
  "IU_PDA_NP11_t" = "IU_PDA_HM9",
"IU_PDA_NP2_t" = "IU_PDA_HM10",
"IU_PDA_T11_t" = "IU_PDA_HM11",
"IU_PDA_T3_t" = "IU_PDA_HM12",
"IU_PDA_T8_t" = "IU_PDA_HM13",
"IU_PDA_HM13_t" = "IU_PDA_HM2",
"IU_PDA_HM2_t" = "IU_PDA_HM3",
"IU_PDA_HM4_t" = "IU_PDA_HM4",
"IU_PDA_HM6_t" = "IU_PDA_HM5",
"IU_PDA_HM8_t" = "IU_PDA_HM6",
"IU_PDA_LNM7_t" = "IU_PDA_HM8",
"IU_PDA_T1_t" = "IU_PDA_LNM10",
"IU_PDA_T4_t" = "IU_PDA_LNM12",
"IU_PDA_LNM10_t" = "IU_PDA_LNM6",
"IU_PDA_LNM6_t" = "IU_PDA_LNM7",
"IU_PDA_LNM8_t" = "IU_PDA_LNM8",
"IU_PDA_HM12_t" = "IU_PDA_HM2_2",
"IU_PDA_T9_t" = "IU_PDA_NP10",
"IU_PDA_T2_t" = "IU_PDA_NP11",
"IU_PDA_HM10_t" = "IU_PDA_NP2",
"IU_PDA_HM9_t" = "IU_PDA_T1",
"IU_PDA_NP10_t" = "IU_PDA_T9",
"IU_PDA_T10_t" = "IU_PDA_T10",
"IU_PDA_T12_t" = "IU_PDA_T11",
"IU_PDA_T6_t" = "IU_PDA_T12",
"IU_PDA_HM11_t" = "IU_PDA_T2",
"IU_PDA_HM3_t" = "IU_PDA_T3",
"IU_PDA_HM5_t" = "IU_PDA_T4",
"IU_PDA_LNM12_t" = "IU_PDA_T6",
"IU_PDA_HM2_2_t" = "IU_PDA_T8"
)

# Rename the image names in pdac_allres
names(pdac@images) <- image_mapping[names(pdac@images)]


## Stacked Bar Plots for Cell Type Percentages by orig.ident


install.packages("ggplot2")
install.packages("dplyr")


library(ggplot2)
library(dplyr)

#updated metadata from the file
metadata <- read.csv("/Users/akhaliq/Desktop/spatial_analysis/st_new_seurat/all_new_mod/meta.csv")
prop <- read.csv("df.csv")
merged_data <- merge(metadata, prop , by="barcode")


cell_types <- c("CAF", "Monocytes", "B_cells", "TAM", "Hepatocytes", 
                "Endothelial_cells", "PVL", "Normal_Epithelial_cells", 
                "T_NK_cells", "Tumor_Epithelial_cells")

percentage_data <- merged_data %>%
  group_by(orig.ident) %>%
  summarize(across(all_of(cell_types), sum, na.rm = TRUE)) %>%
  tidyr::pivot_longer(cols = -orig.ident, names_to = "cell_type", values_to = "count") %>%
  group_by(orig.ident) %>%
  mutate(percentage = count / sum(count) * 100)


percentage_data$cell_type <- factor(percentage_data$cell_type, levels = cell_types)


bar_plot <- ggplot(percentage_data, aes(x = orig.ident, y = percentage, fill = cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Set3") +
  labs(x = "orig.ident", y = "Percentage", title = "Stacked Bar Plot for Each orig.ident") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


ggsave("stacked_bar_plots.pdf", plot = bar_plot, width = 10, height = 6)


# To plot patient/origin/hsitologg

# Install the required packages if not already installed
install.packages("ggplot2")
install.packages("dplyr")
install.packages("RColorBrewer")

# Load the required libraries
library(ggplot2)
library(dplyr)
library(RColorBrewer)

# Read the updated metadata from the file
metadata <- read.csv("/Users/akhaliq/Desktop/spatial_analysis/st_new_seurat/all_new_mod/meta.csv")

# Merge the metadata with the merged_data dataframe
merged_data <- merge(merged_data, metadata, by = "orig.ident")

# Define the order of cell types for stacking
cell_types <- c("CAF", "Monocytes", "B_cells", "TAM", "Hepatocytes", 
                "Endothelial_cells", "PVL", "Normal_Epithelial_cells", 
                "T_NK_cells", "Tumor_Epithelial_cells")

# Calculate the percentage of each cell type within an 'orig.ident'
percentage_data <- merged_data %>%
  group_by(orig.ident, Origin, Histology, patient) %>%
  summarize(across(all_of(cell_types), sum, na.rm = TRUE)) %>%
  tidyr::pivot_longer(cols = -c(orig.ident, Origin, Histology, patient), 
                      names_to = "cell_type", values_to = "count") %>%
  group_by(orig.ident, Origin, Histology, patient) %>%
  mutate(percentage = count / sum(count) * 100)

# Reorder the levels of 'cell_type' based on 'cell_types'
percentage_data$cell_type <- factor(percentage_data$cell_type, levels = cell_types)

# Create a stacked bar plot for each 'orig.ident' with vibrant colors
bar_plot <- ggplot(percentage_data, aes(x = orig.ident, y = percentage, fill = cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = brewer.pal(n = length(unique(percentage_data$cell_type)), name = "Set3")) +
  labs(x = "orig.ident", y = "Percentage", title = "Stacked Bar Plot for Each orig.ident") +
  facet_grid(Origin ~ Histology ~ patient, scales = "free", labeller = label_wrap_gen(width = 10)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Save the plot as a PDF
ggsave("stacked_bar_plots.pdf", plot = bar_plot, width = 12, height = 8)

#### other way of representation

# Install the required packages if not already installed
install.packages("ggplot2")
install.packages("dplyr")
install.packages("RColorBrewer")

# Load the required libraries
library(ggplot2)
library(dplyr)
library(RColorBrewer)

# Read the updated metadata from the file
metadata <- read.csv("/Users/akhaliq/Desktop/spatial_analysis/st_new_seurat/all_new_mod/meta.csv")

# Merge the metadata with the merged_data dataframe
prop <- read.csv("df.csv")
merged_data <- merge(metadata, prop, by = "barcode")

# Define the order of cell types for stacking
cell_types <- c("CAF", "Monocytes", "B_cells", "TAM", "Hepatocytes", 
                "Endothelial_cells", "PVL", "Normal_Epithelial_cells", 
                "T_NK_cells", "Tumor_Epithelial_cells")

# Calculate the percentage of each cell type within an 'orig.ident'
percentage_data <- merged_data %>%
  group_by(orig.ident, Origin, Histology, patient) %>%
  summarize(count = sum(count)) %>%
  group_by(orig.ident, Histology) %>%
  mutate(percentage = count / sum(count) * 100)

# Reorder the levels of 'cell_type' based on 'cell_types'
percentage_data$cell_type <- factor(percentage_data$cell_type, levels = cell_types)

# Create a stacked bar plot for each 'orig.ident' grouped by 'Histology' with vibrant colors
bar_plot <- ggplot(percentage_data, aes(x = orig.ident, y = percentage, fill = cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = brewer.pal(n = length(unique(percentage_data$cell_type)), name = "Set3")) +
  labs(x = "orig.ident", y = "Percentage", title = "Stacked Bar Plot for Each orig.ident") +
  facet_grid(.~ Histology, scales = "free", labeller = label_wrap_gen(width = 10)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Save the plot as a PDF
ggsave("stacked_bar_plots.pdf", plot = bar_plot, width = 12, height = 8)

##

# pie chart
install.packages("ggplot2")
install.packages("dplyr")

library(ggplot2)
library(dplyr)

# Read metadata and property data
metadata <- read.csv("/Users/akhaliq/Desktop/spatial_analysis/st_new_seurat/all_new_mod/meta.csv")
prop <- read.csv("df.csv")

# Merge metadata and property data
merged_data <- merge(metadata, prop, by = "barcode")

# Define cell types
cell_types <- c("CAF", "Monocytes", "B_cells", "TAM", "Hepatocytes", 
                "Endothelial_cells", "PVL", "Normal_Epithelial_cells", 
                "T_NK_cells", "Tumor_Epithelial_cells")

# Calculate cell type percentages by orig.ident
percentage_data <- merged_data %>%
  group_by(orig.ident) %>%
  summarize(across(all_of(cell_types), sum, na.rm = TRUE)) %>%
  tidyr::pivot_longer(cols = -orig.ident, names_to = "cell_type", values_to = "count") %>%
  group_by(orig.ident) %>%
  mutate(percentage = count / sum(count) * 100)

# Set the order of cell types
percentage_data$cell_type <- factor(percentage_data$cell_type, levels = cell_types)

# Create pie plots for each orig.ident
pie_plots <- lapply(unique(percentage_data$orig.ident), function(orig) {
  data_subset <- subset(percentage_data, orig.ident == orig)
  
  pie_plot <- ggplot(data_subset, aes(x = "", y = percentage, fill = cell_type)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    scale_fill_brewer(palette = "Set3") +
    labs(title = orig, fill = "Cell Type") +
    theme_void()
  
  pie_plot
})

# Save pie plots in a single PDF file
pdf("pie_plots.pdf", width = 10, height = 6)
for (plot in pie_plots) {
  print(plot)
}
dev.off()

###

## to get Assays from pdac_allres object 

assay_matrix <- pdac_allres[["rctd_full"]]@data
assay_df <- as.data.frame(t(assay_matrix))
write.csv(assay_df, file = "assay_rctd_full.csv")


######## ISCHIA


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





## to get Assays from pdac_allres object 

assay_matrix <- pdac_allres[["rctd_fullfinal"]]@data
norm_weights <- as.data.frame(t(assay_matrix))
write.csv(assay_df, file = "assay_rctd_full.csv")



# Deciding about the k
Composition.cluster.k(deconv.mat, 20)

# Composition clustering of the deconvoluted spatial spots
pdac_allres <- Composition.cluster(pdac_allres,norm_weights, 14)
#not good...Idents(pdac_allres) <- "CompositionCluster_CC"

table(pdac_allres$CompositionCluster_CC)
#> 
#> CC1 CC2 CC3 CC4 CC5 CC6 CC7 CC8 
#> 302 296 274 382 108 238 426 159

SpatialDimPlot(pdac_allres, group.by = c("CompositionCluster_CC")) + scale_fill_manual(values = c("cyan", "orange", "purple","green","yellow","blue", "red","black"))

## dim for all images in a loop
"#7FC97F", "#8F30A1", "#FDC086", "#BF5B17", "#17DEEE", "#F0027F", "#666666", "#FFD700", "#1B9E77", "#b9bf17", "#FF00FF", "#377EB8", "#4DAF4A", "#ff9a00", "#FF0000"
library(ggplot2)
library(showtext)
library(gridExtra)

# Set the theme for the plot
theme_set(theme_minimal(base_family = "Arial"))

# Load the font for bold title
font_add_google("Montserrat", "Montserrat-Bold")

# Enable the loaded font
showtext_auto()

# Define the image names
image_names <- c("IU_PDA_HM9", "IU_PDA_HM10", "IU_PDA_HM11", "IU_PDA_HM12", "IU_PDA_HM13", "IU_PDA_HM2",
                 "IU_PDA_HM3", "IU_PDA_HM4", "IU_PDA_HM5", "IU_PDA_HM6", "IU_PDA_HM8", "IU_PDA_LNM10",
                 "IU_PDA_LNM12", "IU_PDA_LNM6", "IU_PDA_LNM7", "IU_PDA_LNM8", "IU_PDA_NH2", "IU_PDA_NP10",
                 "IU_PDA_NP11", "IU_PDA_NP2", "IU_PDA_T1", "IU_PDA_T9", "IU_PDA_T10", "IU_PDA_T11",
                 "IU_PDA_T12", "IU_PDA_T2", "IU_PDA_T3", "IU_PDA_T4", "IU_PDA_T6", "IU_PDA_T8")

# Define the color palette
color_palette <- c("#FF00FF", "#FF0000", "#FFA500", "#FFFF00", "#00FF00", "#00FFFF", "#0000FF", "#800080", "#FFD700", "#FF4500", "#F7F6C5", "#00FF7F", "#9400D3", "#00CED1","#FF0000")


# Set the output PDF file
output_file <- "spatial_plots_K10_new.pdf"

# Create a new PDF file
pdf(output_file, onefile = TRUE)

# Iterate over the image names
for (image_name in image_names) {
  # Create the plot
  plot <- SpatialDimPlot(pdac_allres, group.by = c("CompositionCluster_CC"), images = image_name) +
    scale_fill_manual(values = color_palette)
  
  # Add a title with image name in bold
  plot <- plot + ggtitle(bquote(bold(.(image_name))))
  
  # Save the plot on a separate page of the PDF file
  print(plot)
}

# Close the PDF file
dev.off()

########
Ask him to give aderinal to pancreatic sections,
#> Scale for fill is already present.
#> Adding another scale for fill, which will replace the existing scale.


#Composition_cluster_enrichedCelltypes(pdac_allres,"CC4", as.matrix(norm_weights))

### for All CCs

#c("CC1", "CC2", "CC3", "CC4", "CC5", "CC6", "CC7", "CC8", "CC9", "CC10", "CC11","CC12","CC13","CC14")

library(gridExtra)

# Define a function to generate and save the plot for a specific CC
save_cc_plot <- function(cc) {
  # Generate the plot for the given CC
  plot <- Composition_cluster_enrichedCelltypes(pdac_allres, cc, as.matrix(norm_weights))
  
  # Create a PDF file for the current CC plot
  pdf_name <- paste0(cc, ".pdf")
  pdf(file = pdf_name)
  
  # Save the plot to the PDF
  print(plot)
  
  # Close the PDF file
  dev.off()
}

# Iterate over each CC and save its respective plot
ccs <- paste0("CC", 1:10)  # Replace with the appropriate range for your case
for (cc in ccs) {
  save_cc_plot(cc)
}

library(pdftools)

# Define the names of the individual PDF files
pdf_files <- paste0("CC", 1:10, ".pdf")  # Replace with the appropriate range for your case

# Define the name of the merged PDF file
merged_pdf <- "enrichedCelltypes_cc_ischia_10.pdf"

# Use the pdf_combine function to merge the PDFs
pdf_combine(pdf_files, output = merged_pdf)


pdac_allres.umap <- Composition_cluster_umap(pdac_allres, norm_weights)
#> Plotting scatterpies for 2185 pixels with 20 cell-types...this could take a while if the dataset is large.
pdac_allres.umap$umap.cluster.gg

pdf("pie_chart.pdf")
pdac_allres.umap$umap.deconv.gg
dev.off()



#####
## adding UMAP and TSNE to the seurat obj

#write.csv(emb,"tsne_stdcon.csv")
#write.csv(emb.umap,"umap_stdcon.csv")
write.csv(pdac_allres.umap$umap.table,"umap_ischia14.csv")
#emb.tsne <- read.csv("tsne_stdcon.csv",header=T,sep=",",row.names=1)
emb.umap<- read.csv("umap_ischia14.csv",header=T,sep=",",row.names=1)

emb.umap$CompositionCluster_CC <- NULL
emb.umap$Slide <- NULL

emb.umap <- as.matrix(emb.umap)

#colnames(emb.tsne)<- c("tSNE1","tSNE2")
colnames(emb.umap)<- c("UMAP1","UMAP2")
head(emb.umap)



#pdac[['tsne.std']] <- CreateDimReducObject(embeddings = as.matrix(emb.tsne), key = 'tsne_', assay = 'integrated')
#pdac[['umap.std']] <- CreateDimReducObject(embeddings = as.matrix(emb.umap1), key = 'umap_', assay = 'integrated')

pdac_allres[['umap.ischia10']] <- CreateDimReducObject(embeddings = emb.umap,key = 'umap.ischia10_', assay = 'rctd_fullfinal')


emb.umap<- read.csv("umap_ischia10.csv",header=T,sep=",",row.names=1)
emb.umap$x <- NULL
emb.umap$y <- NULL
emb.umap$Slide <- NULL
colnames(emb.umap) <- "cc_ischia_10"
pdac_allres <- AddMetaData(pdac_allres,emb.umap)

pdf("seurat_ischia_umap_10.pdf")
DimPlot(pdac_allres, reduction = "umap.ischia10", label = FALSE,group.by="cc_ischia_10")
dev.off()

pdf("barplot_SampVsorig_10.pdf")
dittoBarPlot(pdac_allres, "orig.ident", group.by = "cc_ischia_10")
dev.off()

pdf("barplot_origVsSamp_10.pdf")
dittoBarPlot(pdac_allres, "cc_ischia_10", group.by = "orig.ident")
dev.off()

per <- dittoBarPlot(pdac_allres, "orig.ident", group.by = "cc_ischia_10", data.out = TRUE)
write.csv(per$data,"percent_K10.csv")

CC4.celltype.cooccur <- spatial.celltype.cooccurence(spatial.object=pdac_allres,deconv.prob.mat=deconv.mat, COI="CC4", prob.th= 0.05, Condition=unique(pdac_allres$orig.ident))
plot.celltype.cooccurence(CC4.celltype.cooccur)

###
"Here's how the test is applied in your code:

Grouping Data: The data is grouped based on the Histology variable, which represents different histology categories (e.g., Liver_Mets, Lymph node, Normal Liver, Normal Pancreas, PDAC).

Iteration over Histology Categories: The code uses the lapply function to iterate over each unique histology category. This allows you to perform the Wilcoxon rank-sum test for each category separately.

Subset Data for Current Category: Within each iteration, the merged data is subsetted to include only the data corresponding to the current histology category. This ensures that the test is performed on the appropriate subset of data for each category.

Performing the Wilcoxon Rank-Sum Test: For each cell type, the code performs a Wilcoxon rank-sum test between the data of the current histology category and the entire merged data. This test evaluates whether the distribution of normalized weights for a particular cell type differs significantly between the current category and the rest of the data.

Calculating P-Values: The resulting test statistic and p-value are obtained from the Wilcoxon rank-sum test. The p-value represents the probability of observing the observed test statistic or a more extreme value under the null hypothesis that there is no difference in the distributions of normalized weights between the two groups.

Multiple Comparison Adjustment: To account for the issue of multiple comparisons (performing multiple statistical tests), the p-values are adjusted using the Benjamini-Hochberg procedure. This procedure controls the false discovery rate (FDR) and helps identify statistically significant results while minimizing the chances of false positives.

By performing these tests for each histology category and each cell type, the code allows you to identify which cell types have significantly different normalized weights between different histology categories. The adjusted p-values provide a measure of statistical significance, indicating which comparisons are likely to be true positives. This information can help you understand the relationships between histology and cell type-specific weights in your data."


# Merge norm_weights with meta based on row names
merged_data <- merge(norm_weights, meta, by = 0, sort = FALSE)

# Perform Wilcoxon rank-sum test for each histology category
results <- lapply(unique(meta$cc_ischia_10), function(category) {
  # Subset the merged data for the current category
  category_data <- merged_data[merged_data$cc_ischia_10 == category, ]
  
  # Perform Wilcoxon rank-sum test on each cell type
  p_values <- sapply(names(norm_weights)[-1], function(Histology) {
    test <- wilcox.test(category_data[[Histology]], merged_data[[Histology]])
    test$p.value
  })
  
  # Adjust p-values for multiple comparisons using the Benjamini-Hochberg procedure
  adjusted_p_values <- p.adjust(p_values, method = "BH")
  
  # Return the adjusted p-values for the current category
  data.frame(Category = category, Cell_Type = names(norm_weights)[-1], Adjusted_P_Value = adjusted_p_values)
})

# Combine the results into a single data frame
results_df <- do.call(rbind, results)

# Print the results
print(results_df)

library(tidyverse)

# Reshape the data to long format
long_data <- merged_data %>%
  pivot_longer(cols = c("B cells", "C1Q-TAM", "CD4+ cells", "CD8-NK cells", "DCs", "Endothelial cells",
                        "FCN1-TAM", "Hepatocytes", "iCAF", "myCAF", "Normal Epithelial cells",
                        "Proliferative T cells", "PVL", "SPP1-TAM", "Tumor Epithelial cells"),
               names_to = "Cell_Type", values_to = "Normalized_Weight")

# Define colors for the box plots
cell_type_colors <- c("steelblue", "darkorange", "darkgreen", "purple", "red", "blue", "green",
                      "cyan", "magenta", "brown", "gold", "pink", "grey", "olivedrab", "dodgerblue")

# Plot all cell types on a single page with improved aesthetics
plot <- ggplot(long_data, aes(x = Histology, y = Normalized_Weight, fill = Cell_Type)) +
  geom_boxplot(outlier.shape = NA) +
  labs(x = "Histology", y = "Normalized Weight") +
  ggtitle("Wilcoxon Rank-Sum Test Results by Cell Type and Histology") +
  theme_bw() +
  scale_fill_manual(values = cell_type_colors) +
  facet_grid(Cell_Type ~ ., scales = "free_y", space = "free", switch = "y") +
  theme(strip.text.y = element_text(angle = 0, vjust = 0.5),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 20)),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        panel.spacing = unit(0.5, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.background = element_blank(),
        plot.margin = margin(1, 1, 1, 1, "cm"))

# Save the plot to PDF with high resolution and suitable dimensions for publication
ggsave(filename = "box_plots_all_cell_types_Test.pdf", plot = plot, width = 10, height = 12, units = "in", dpi = 300)
#### Heatmaop of the same

library(tidyverse)

# Reshape the data to long format
long_data <- merged_data %>%
  pivot_longer(cols = c("B cells", "C1Q-TAM", "CD4+ cells", "CD8-NK cells", "DCs", "Endothelial cells",
                        "FCN1-TAM", "Hepatocytes", "iCAF", "myCAF", "Normal Epithelial cells",
                        "Proliferative T cells", "PVL", "SPP1-TAM", "Tumor Epithelial cells"),
               names_to = "Cell_Type", values_to = "Normalized_Weight")

# Create a custom color palette
my_palette <- colorRampPalette(c("#FFECF0", "#FFA5B2", "#FF627C", "#FF1744", "#B3002D"))

# Create a heatmap
heatmap_plot <- ggplot(long_data, aes(x = Cell_Type, y = Histology, fill = Normalized_Weight)) +
  geom_tile() +
  labs(x = "Cell Type", y = "Histology", fill = "Normalized Weight") +
  ggtitle("Heatmap of Normalized Weights") +
  theme_bw() +
  scale_fill_gradientn(colors = my_palette(100)) +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10, angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10))

# Save the heatmap to PDF with high resolution and suitable dimensions for publication
ggsave(filename = "heatmap_normalized_weights_cc.pdf", plot = heatmap_plot, width = 8, height = 10, units = "in", dpi = 300)

#####
Stacked_bar monad_plots

library(tidyverse)

# Define distinct colors for each cell type
cell_type_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
                      "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#1b9e77", "#d95f02",
                      "#7570b3", "#e7298a", "#66a61e")

# Reshape the data to long format
long_data <- merged_data %>%
  pivot_longer(cols = c("B cells", "C1Q-TAM", "CD4+ cells", "CD8-NK cells", "DCs", "Endothelial cells",
                        "FCN1-TAM", "Hepatocytes", "iCAF", "myCAF", "Normal Epithelial cells",
                        "Proliferative T cells", "PVL", "SPP1-TAM", "Tumor Epithelial cells"),
               names_to = "Cell_Type", values_to = "Normalized_Weight")

# Calculate the percentage of each cell type within histology categories
percentages <- long_data %>%
  group_by(Histology, Cell_Type) %>%
  summarize(Percentage = sum(Normalized_Weight)) %>%
  group_by(Histology) %>%
  mutate(Percentage = Percentage / sum(Percentage) * 100)

# Create a stacked bar plot with distinct colors and improved aesthetics
stacked_bar_plot <- ggplot(percentages, aes(x = Histology, y = Percentage, fill = Cell_Type)) +
  geom_bar(stat = "identity") +
  labs(x = "Histology", y = "Percentage", fill = "Cell Type") +
  ggtitle("Percentage of Cell Types within CC Categories") +
  theme_bw() +
  scale_fill_manual(values = cell_type_colors) +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        axis.line = element_line(color = "black"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.position = "right",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Add x and y axis lines
stacked_bar_plot <- stacked_bar_plot + 
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.grid.major.y = element_line(color = "gray", linetype = "dashed", size = 0.5),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank())

# Save the stacked bar plot with distinct colors to PDF with high resolution and suitable dimensions for publication
ggsave(filename = "stacked_bar_plot_percentages_cc.pdf", plot = stacked_bar_plot, width = 8, height = 6, units = "in", dpi = 300)


#####################

Wilcoxon test between Ischia cc and Histology

## Enrichment Analysis and Box Plot Generation

This code performs enrichment analysis and generates box plots for a given dataset. It helps to identify significant differences in proportions between different categories of a variable.

### Requirements

The following R packages are required to run this code:

- ggplot2
- gridExtra

Please ensure that these packages are installed before running the code.

### Usage

1. Prepare the Data:

   - Make sure you have the following data:
     - `norm_weights`: A data frame containing normalized weights.
     - `meta`: A data frame containing metadata, including the columns "Histology" and "cc_ischia_10".

2. Run the Code:

   - Copy and paste the provided code into your R environment or script.
   - Make sure to replace `norm_weights` and `meta` with your actual data frames.
   - Run the code.

3. Interpret the Results:

   - The code will perform enrichment analysis and generate box plots for each histology category.
   - The box plots show the distribution of p-values for each `cc_ischia_10` category.
   - Significant differences in proportions between different categories are indicated by lower p-values.
   - The box plots are saved as a PDF file named "boxplots.pdf" in your working directory.

### Additional Notes

- The code uses the `merge()` function to combine the `norm_weights` and `meta` data frames based on the spot IDs.
- It then subsets the merged data frame to include only the relevant columns for analysis.
- Enrichment analysis is performed by comparing proportions using the Wilcoxon test.
- The code uses custom colors for the box plots. You can modify the `colors` vector to customize the colors as per your preference.
- The resulting box plots are arranged in a grid using the `gridExtra` package.
- The grid of box plots is saved as a PDF file using the `ggsave()` function.

###############################

meta <- data.frame(pdac_allres@meta.data)
# Merge the norm_weights dataframe with the meta dataframe based on the spot IDs
merged_df <- merge(norm_weights, meta, by = "row.names", all.x = TRUE)

# Subset the merged dataframe to include only the relevant columns
subset_df <- subset(merged_df, select = c("Histology", "cc_ischia_10", colnames(norm_weights)))

# Create an empty list to store the enrichment results
enrichment_results <- list()

# Iterate over each histology category
for (histology in unique(subset_df$Histology)) {
  # Subset the data for the current histology category
  histology_subset <- subset_df[subset_df$Histology == histology, ]

  # Create a matrix to store the Wilcoxon test scores for each cc_ischia_10 category
  scores <- matrix(NA, nrow = length(unique(histology_subset$cc_ischia_10)),
                   ncol = ncol(norm_weights),
                   dimnames = list(unique(histology_subset$cc_ischia_10), colnames(norm_weights)))

  # Iterate over each cc_ischia_10 category
  for (cc_category in unique(histology_subset$cc_ischia_10)) {
    # Subset the data for the current cc_ischia_10 category
    cc_subset <- histology_subset[histology_subset$cc_ischia_10 == cc_category, ]

    # Iterate over each proportion column and perform Wilcoxon test
    for (col in colnames(norm_weights)) {
      group1 <- cc_subset[, col]
      group2 <- subset_df[subset_df$Histology != histology & subset_df$cc_ischia_10 == cc_category, col]

      # Check if both groups have enough observations for the test
      if (length(group1) > 1 && length(group2) > 1) {
        # Perform Wilcoxon test and store the score (e.g., U statistic)
        score <- wilcox.test(group1, group2)$statistic
        scores[cc_category, col] <- score
      } else {
        # Set the score to NA if there are not enough observations
        scores[cc_category, col] <- NA
      }
    }
  }

  # Store the enrichment results for the current histology category
  enrichment_results[[histology]] <- scores
}



library(ggplot2)
library(gridExtra)

# Create a list to store the box plot data
boxplot_data <- list()

# Define custom colors for the box plots
colors <- c("Liver_Mets" = "#FF0000", "Lymph node" = "#00FF00", "Normal Pancreas" = "#0000FF", "PDAC" = "#FF00FF")

# Iterate over each histology category
for (histology in names(enrichment_results)) {
  # Get the enrichment results for the current histology category
  histology_results <- enrichment_results[[histology]]
  
  # Create a data frame with the p-values for each cc_ischia_10 category
  p_values_df <- data.frame(cc_ischia_10 = rownames(histology_results),
                            p_value = as.numeric(histology_results))
  
  # Remove rows with NA values
  p_values_df <- na.omit(p_values_df)
  
  # Store the data frame in the box plot data list
  boxplot_data[[histology]] <- p_values_df
}

# Create a list to store the plots
plots <- list()

# Create box plots for each histology category
for (histology in names(boxplot_data)) {
  # Create a new plot for each histology category
  boxplot_plot <- ggplot(boxplot_data[[histology]], aes(x = cc_ischia_10, y = p_value)) +
    geom_boxplot(fill = colors[histology], color = "black", alpha = 0.8, outlier.shape = NA) +
    labs(title = histology, x = "Composition Clusters", y = "Wilcox Score") +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  
  # Add the plot to the list of plots
  plots[[histology]] <- boxplot_plot
}

# Combine the plots into a single page with extra spacing
grid_plots <- grid.arrange(grobs = plots, ncol = 2, top = "Enrichment Analysis", padding = unit(2, "lines"))

# Save the plots as a PDF file
ggsave("Wilcoxon_CCVsHistology_16.pdf", grid_plots, width = 12, height = 10)

# Print the success message
cat("Box plots saved as 'boxplots_16.pdf'\n")

#########  ***** for individual

  # Merge the norm_weights dataframe with the meta dataframe based on the spot IDs
  merged_df <- merge(norm_weights, meta, by = "row.names", all.x = TRUE)

  # Subset the merged dataframe to include only the relevant columns
  subset_df <- subset(merged_df, select = c("cc_ischia_10", "Histology", colnames(norm_weights)))

  # Create an empty list to store the enrichment results
  enrichment_results <- list()

  # Iterate over each cc_ischia_10 category
  for (cc_category in unique(subset_df$cc_ischia_10)) {
    # Subset the data for the current cc_ischia_10 category
    cc_subset <- subset_df[subset_df$cc_ischia_10 == cc_category, ]
    
    # Create a matrix to store the Wilcoxon test statistics for each Histology category
    wilcoxon_stats <- matrix(NA, nrow = length(unique(cc_subset$Histology)),
                             ncol = ncol(norm_weights),
                             dimnames = list(unique(cc_subset$Histology), colnames(norm_weights)))
    
    # Iterate over each Histology category
    for (histology in unique(cc_subset$Histology)) {
      # Subset the data for the current Histology category
      histology_subset <- cc_subset[cc_subset$Histology == histology, ]
      
      # Iterate over each proportion column and perform Wilcoxon test
      for (col in colnames(norm_weights)) {
        group1 <- histology_subset[, col]
        group2 <- subset_df[subset_df$cc_ischia_10 != cc_category & subset_df$Histology == histology, col]
        
        # Check if both groups have enough observations for the test
        if (length(group1) > 1 && length(group2) > 1) {
          # Perform Wilcoxon test and store the test statistic
          wilcoxon_stat <- wilcox.test(group1, group2)$statistic
          wilcoxon_stats[histology, col] <- wilcoxon_stat
        } else {
          # Set the test statistic to NA if there are not enough observations
          wilcoxon_stats[histology, col] <- NA
        }
      }
    }
    
    # Store the Wilcoxon test statistics in the enrichment results list
    enrichment_results[[cc_category]] <- wilcoxon_stats
  }


  library(ggplot2)
  library(gridExtra)

# Create a list to store the plots
plots <- list()

# Define the histology categories and their corresponding colors
histology_colors <- c("Liver_Mets" = "skyblue", "Lymph node" = "pink","Normal Pancreas" = "orange", "PDAC" = "purple")

# Create box plots for each orig.ident category
for (orig_ident in names(enrichment_results)) {
  # Get the enrichment results for the current orig.ident category
  orig_ident_results <- enrichment_results[[orig_ident]]
  
  # Create a data frame with the p-values for each cc_ischia_10 category
  p_values_df <- data.frame(cc_ischia_10 = rownames(orig_ident_results),
                            p_value = as.numeric(orig_ident_results))
  
  # Remove rows with NA values
  p_values_df <- na.omit(p_values_df)
  
  # Create a new plot for each orig.ident category
  boxplot_plot <- ggplot(p_values_df, aes(x = cc_ischia_10, y = p_value, fill = cc_ischia_10)) +
    geom_boxplot(color = "black", alpha = 0.8, outlier.shape = NA) +
    geom_text(data = subset(p_values_df, p_value < 0.05), aes(label = sprintf("%.2f", p_value)),
              vjust = -0.5, size = 3, color = "black", show.legend = FALSE) +
    labs(title = orig_ident, x = "Origin", y = "Proportions") +
    scale_fill_manual(values = histology_colors, labels = names(histology_colors)) +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          legend.position = "none")  # Hide the legend to avoid duplicate labels
  
  # Add the plot to the list of plots
  plots[[orig_ident]] <- boxplot_plot
}

  # Determine the number of plots and calculate the grid layout dimensions
  num_plots <- length(plots)
  num_cols <- 2
  num_rows <- ceiling(num_plots / num_cols)

  # Calculate the plot width and height based on desired dimensions
  plot_width <- 6
  plot_height <- 5

  # Adjust the plot size and aspect ratio
  for (orig_ident in names(plots)) {
    plots[[orig_ident]] <- plots[[orig_ident]] +
      theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
      theme(plot.background = element_blank(),
            panel.background = element_blank())
    
    ggsave(paste0(orig_ident, ".pdf"), plots[[orig_ident]], width = plot_width, height = plot_height)
  }

  # Combine all the plots into a single PDF file
  grid_plots <- do.call(grid.arrange, c(plots, ncol = num_cols))

  # Save the combined plots as a PDF file
  ggsave("Combined_Plots_16.pdf", grid_plots, width = num_cols * plot_width, height = num_rows * plot_height)

  # Print the success message
  cat("Combined plots saved as 'Combined_Plots_16.pdf'\n")

### summary for the above

"In the provided code, we are performing the Wilcoxon rank sum test (also known as the Mann-Whitney U test) using the wilcox.test function. The Wilcoxon test is a non-parametric test that compares two independent samples to determine if they come from the same distribution. It does not assume any specific distribution for the data.

Here's a breakdown of the steps involved in the Wilcoxon test within the code:

The code iterates over each histology category.
For each histology category, it subsets the data to include only the relevant columns for that histology.
It creates an empty matrix (p_values) to store the p-values for each cc_ischia_10 category and each column in norm_weights.
It further subsets the data for the current cc_ischia_10 category.
It then iterates over each proportion column and performs a Wilcoxon test between the current cc_ischia_10 group and the groups from other histology categories.
If both groups have enough observations (more than 1), it performs the Wilcoxon test using wilcox.test and stores the resulting p-value in the p_values matrix.
If either group has insufficient observations (less than or equal to 1), it sets the p-value to NA.
Finally, it stores the enrichment results (p-values) for each histology category in the enrichment_results list.
The Wilcoxon test is used to compare the distributions of the group1 and group2 samples for each proportion column. The resulting p-values indicate the statistical significance of the difference between the two groups. The lower the p-value, the stronger the evidence against the null hypothesis, suggesting a significant difference between the distributions.

In summary, the code performs Wilcoxon tests between different groups based on cc_ischia_10 categories and stores the resulting p-values, allowing for the comparison of distributions and identification of statistically significant differences between the groups."

##### Get the P value for the Wilcoxon 

No asterisk: p-value >= 0.05 (not statistically significant)
One asterisk (*): 0.01 <= p-value < 0.05 (statistically significant at the 0.05 level)
Two asterisks (**): 0.001 <= p-value < 0.01 (statistically significant at the 0.01 level)
Three asterisks (***): p-value < 0.001 (statistically significant at the 0.001 level)

For example:

The p-value for "B cells" is ***0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000012, which is less than 0.001, indicating a highly statistically significant result.
The p-value for "C1Q-TAM" is ***0.0078, which is less than 0.01 but greater than 0.001, indicating statistical significance at the 0.01 level.
The p-value for "DCs" is ***0.000000000000000000000000000000000000000000000000000000000000078, which is less than 0.001, indicating a highly statistically significant result.

####

library(openxlsx)

# Merge the norm_weights dataframe with the meta dataframe based on the spot IDs
merged_df <- merge(norm_weights, meta, by = "row.names", all.x = TRUE)

# Subset the merged dataframe to include only the relevant columns
subset_df <- subset(merged_df, select = c("cc_ischia_10", "Histology", colnames(norm_weights)))

# Create an empty list to store the enrichment results
enrichment_results <- list()

# Iterate over each cc_ischia_10 category
for (cc_category in unique(subset_df$cc_ischia_10)) {
  # Subset the data for the current cc_ischia_10 category
  cc_subset <- subset_df[subset_df$cc_ischia_10 == cc_category, ]
  
  # Create a matrix to store the p-values for each Histology category
  p_values <- matrix(NA, nrow = length(unique(cc_subset$Histology)),
                     ncol = ncol(norm_weights),
                     dimnames = list(unique(cc_subset$Histology), colnames(norm_weights)))
  
  # Iterate over each Histology category
  for (histology in unique(cc_subset$Histology)) {
    # Subset the data for the current Histology category
    histology_subset <- cc_subset[cc_subset$Histology == histology, ]
    
    # Iterate over each proportion column and perform Wilcoxon test
    for (col in colnames(norm_weights)) {
      group1 <- histology_subset[, col]
      group2 <- subset_df[subset_df$cc_ischia_10 != cc_category & subset_df$Histology == histology, col]
      
      # Check if both groups have enough observations for the test
      if (length(group1) > 1 && length(group2) > 1) {
        # Perform Wilcoxon test and extract the p-value
        wilcoxon_result <- wilcox.test(group1, group2)
        p_value <- wilcoxon_result$p.value
        
        # Add asterisks to the p-value if it is significant
        if (p_value < 0.05) {
          asterisks <- paste0("***", format(p_value, scientific = FALSE, digits = 2))
        } else {
          asterisks <- format(p_value, scientific = FALSE, digits = 2)
        }
        
        # Store the p-value with asterisks in the matrix
        p_values[histology, col] <- asterisks
      } else {
        # Set the p-value to NA if there are not enough observations
        p_values[histology, col] <- NA
      }
    }
  }
  
  # Store the p-values in the enrichment results list
  enrichment_results[[cc_category]] <- p_values
}

# Save the enrichment results as an Excel sheet
wb <- createWorkbook()
for (cc_category in names(enrichment_results)) {
  sheet <- addWorksheet(wb, cc_category)
  writeData(wb, sheet, enrichment_results[[cc_category]])
}
saveWorkbook(wb, "enrichment_results.xlsx")

#####


######

Relative Abundance 

###


  library(dplyr)
  library(ggplot2)
  library(viridis)
  library(grid)
  library(gridExtra)

assay_matrix <- pdac[["rctd_fullfinal"]]@data
norm_weights <- as.data.frame(t(assay_matrix))
merged_data <- merge(norm_weights, pdac@meta.data, by = 0, sort = FALSE)

# Define the color palette with 12 colors
distinct_palette <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666", "#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f") # 12 colors
distinct_palette <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666", "#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#cab2d6", "#ff7f00", "#b15928", "#33a02c") # 16 colors
 
# Calculate the percentage of spots per cc_ischia_10 and Histology
spots_per_cc_hist <- merged_data %>%
  filter(Histology %in% c("Liver_Mets", "Lymph node", "Normal Pancreas", "PDAC")) %>%
  group_by(cc_ischia_10, Histology) %>%
  summarise(total_spots = n()) %>%
  group_by(cc_ischia_10) %>%
  mutate(percentage = total_spots / sum(total_spots) * 100)

# Create the stacked barplot
plot <- ggplot(spots_per_cc_hist, aes(x = cc_ischia_10, y = percentage, fill = Histology)) +
  geom_bar(stat = "identity") +
  labs(x = "Composition Clusters", y = "Relative Abundance", fill = "Histology") +
  scale_fill_manual(values = distinct_palette) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot as a PDF file
ggsave("Relative_abundance_stacked_barplot1.pdf", plot, width = 8, height = 6)



# Create the summary table
summary_table <- spots_per_cc_hist %>%
  select(cc_ischia_10, Histology, percentage) %>%
  pivot_wider(names_from = Histology, values_from = percentage) %>%
  arrange(cc_ischia_10)


# Save the summary table as a PDF
pdf("TEST2 ", width = 6, height = 8)
table_grob <- tableGrob(summary_table, theme = ttheme_default(base_size = 8))
grid.newpage()
grid.draw(table_grob)
dev.off()



library(dplyr)
library(ggplot2)
library(viridis)
library(tidyr)

# Filter merged_data to include relevant Histology categories
filtered_data <- merged_data %>%
  filter(Histology %in% c("Liver_Mets", "Lymph node", "Normal Pancreas", "PDAC"))

# Calculate the percentage of spots per cc_ischia_10 and Histology
spots_per_cc_hist <- filtered_data %>%
  group_by(cc_ischia_10, Histology) %>%
  summarise(total_spots = n()) %>%
  group_by(cc_ischia_10) %>%
  mutate(percentage = total_spots / sum(total_spots) * 100)

# Create a complete set of cc_ischia_10 and Histology combinations
complete_data <- expand.grid(cc_ischia_10 = unique(spots_per_cc_hist$cc_ischia_10),
                             Histology = c("Liver_Mets", "Lymph node", "Normal Pancreas", "PDAC"))

# Join complete_data with spots_per_cc_hist to fill in missing combinations
spots_per_cc_hist_complete <- complete_data %>%
  left_join(spots_per_cc_hist, by = c("cc_ischia_10", "Histology")) %>%
  replace_na(list(total_spots = 0, percentage = 0))


# Create the heatmap
heatmap_plot <- ggplot(spots_per_cc_hist_complete, aes(x = cc_ischia_10, y = Histology, fill = percentage)) +
  geom_tile() +
  labs(x = "Composition Clusters", y = "Histology", fill = "Relative Abundance") +
  scale_fill_viridis(option = "turbo", guide = guide_colorbar(barwidth = 15, barheight = 1, title.position = "top"), limits = c(0, 100)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 8),
        legend.key.height = unit(1.5, "cm"),
        legend.key.width = unit(0.5, "cm"),
        panel.grid = element_blank(),
        legend.position = "bottom",  # Place the legend below the plot
        legend.box = "horizontal",  # Display the legend in a horizontal layout
        legend.margin = margin(t = 10))  # Adjust the margin around the legend

# Adjust the size and aspect ratio of the heatmap plot
heatmap_plot <- heatmap_plot + coord_fixed(ratio = 0.5)

# Save the heatmap as a PDF file

ggsave("Relative_abundance_heatmap_new_turbo.pdf", heatmap_plot, width = 8, height = 6)

#####

ditto <- dittoBarPlot(pdac_allres, "orig.ident", group.by = "cc_ischia_10")

ggsave("ditto.pdf", ditto, width = 8, height = 6)

######

#### kruskal_stat

# Merge the norm_weights dataframe with the meta dataframe based on the spot IDs
merged_df <- merge(norm_weights, meta, by = "row.names", all.x = TRUE)

# Subset the merged dataframe to include only the relevant columns
subset_df <- subset(merged_df, select = c("cc_ischia_10", "Histology", colnames(norm_weights)))

# Create an empty list to store the enrichment results
enrichment_results <- list()

# Iterate over each cc_ischia_10 category
for (cc_category in unique(subset_df$cc_ischia_10)) {
  # Subset the data for the current cc_ischia_10 category
  cc_subset <- subset_df[subset_df$cc_ischia_10 == cc_category, ]
  
  # Create a matrix to store the Kruskal-Wallis test statistics for each Histology category
  kruskal_stats <- matrix(NA, nrow = length(unique(cc_subset$Histology)),
                          ncol = ncol(norm_weights),
                          dimnames = list(unique(cc_subset$Histology), colnames(norm_weights)))
  
  # Iterate over each Histology category
  for (histology in unique(cc_subset$Histology)) {
    # Subset the data for the current Histology category
    histology_subset <- cc_subset[cc_subset$Histology == histology, ]
    
    # Iterate over each proportion column and perform Kruskal-Wallis test
    for (col in colnames(norm_weights)) {
      group1 <- histology_subset[, col]
      group2 <- subset_df[subset_df$cc_ischia_10 != cc_category & subset_df$Histology == histology, col]
      
      # Check if both groups have enough observations for the test
      if (length(group1) > 1 && length(group2) > 1) {
        # Perform Kruskal-Wallis test and store the test statistic
        kruskal_stat <- kruskal.test(list(group1, group2))$statistic
        kruskal_stats[histology, col] <- kruskal_stat
      } else {
        # Set the test statistic to NA if there are not enough observations
        kruskal_stats[histology, col] <- NA
      }
    }
  }
  
  # Store the Kruskal-Wallis test statistics in the enrichment results list
  enrichment_results[[cc_category]] <- kruskal_stats
}


library(ggplot2)
library(gridExtra)

# Create a list to store the plots
plots <- list()

# Define the histology categories and their corresponding colors
histology_colors <- c("Liver_Mets" = "skyblue", "Lymph node" = "pink","Normal Pancreas" = "orange", "PDAC" = "purple")

# Create box plots for each orig.ident category
for (orig_ident in names(enrichment_results)) {
  # Get the enrichment results for the current orig.ident category
  orig_ident_results <- enrichment_results[[orig_ident]]
  
  # Create a data frame with the p-values for each cc_ischia_10 category
  p_values_df <- data.frame(cc_ischia_10 = rownames(orig_ident_results),
                            p_value = as.numeric(orig_ident_results))
  
  # Remove rows with NA values
  p_values_df <- na.omit(p_values_df)
  
  # Create a new plot for each orig.ident category
  boxplot_plot <- ggplot(p_values_df, aes(x = cc_ischia_10, y = p_value, fill = cc_ischia_10)) +
    geom_boxplot(color = "black", alpha = 0.8, outlier.shape = NA) +
    geom_text(data = subset(p_values_df, p_value < 0.05), aes(label = sprintf("%.2f", p_value)),
              vjust = -0.5, size = 3, color = "black", show.legend = FALSE) +
    labs(title = orig_ident, x = "Origin", y = "Proportions") +
    scale_fill_manual(values = histology_colors, labels = names(histology_colors)) +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          legend.position = "none")  # Hide the legend to avoid duplicate labels
  
  # Add the plot to the list of plots
  plots[[orig_ident]] <- boxplot_plot
}

# Determine the number of plots and calculate the grid layout dimensions
num_plots <- length(plots)
num_cols <- 2
num_rows <- ceiling(num_plots / num_cols)

# Calculate the plot width and height based on desired dimensions
plot_width <- 6
plot_height <- 5

# Adjust the plot size and aspect ratio
for (orig_ident in names(plots)) {
  plots[[orig_ident]] <- plots[[orig_ident]] +
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
    theme(plot.background = element_blank(),
          panel.background = element_blank())
  
  ggsave(paste0(orig_ident, ".pdf"), plots[[orig_ident]], width = plot_width, height = plot_height)
}

# Combine all the plots into a single PDF file
grid_plots <- do.call(grid.arrange, c(plots, ncol = num_cols))

# Save the combined plots as a PDF file
ggsave("Combined_Plots_KK.pdf", grid_plots, width = num_cols * plot_width, height = num_rows * plot_height)

# Print the success message
cat("Combined plots saved as 'Combined_Plots_KK.pdf'\n")

#### Find the p Value for wilcoxon

"find the p-values between the specified histology groups ("Liver_Mets", "Lymph node", "Normal Pancreas", "PDAC") in all cc_ischia_10 categories, we can modify the code as follows in each CELL TYPES"

# Create a vector of the desired histology groups
histology_groups <- c("Liver_Mets", "Lymph node", "Normal Pancreas", "PDAC")

# Create an empty list to store the p-values
p_values_list <- list()

# Iterate over each cc_ischia_10 category
for (cc_category in unique(subset_df$cc_ischia_10)) {
  # Subset the data for the current cc_ischia_10 category
  cc_subset <- subset_df[subset_df$cc_ischia_10 == cc_category, ]
  
  # Create a matrix to store the p-values for the histology groups
  p_values <- matrix(NA, nrow = length(histology_groups),
                     ncol = ncol(norm_weights),
                     dimnames = list(histology_groups, colnames(norm_weights)))
  
  # Iterate over each histology group
  for (histology in histology_groups) {
    # Subset the data for the current histology group
    histology_subset <- cc_subset[cc_subset$Histology == histology, ]
    
    # Iterate over each proportion column and perform Wilcoxon test
    for (col in colnames(norm_weights)) {
      group1 <- histology_subset[, col]
      group2 <- subset_df[subset_df$cc_ischia_10 != cc_category & subset_df$Histology == histology, col]
      
      # Check if both groups have enough observations for the test
      if (length(group1) > 1 && length(group2) > 1) {
        # Perform Wilcoxon test and store the p-value
        wilcoxon_result <- wilcox.test(group1, group2)
        p_value <- wilcoxon_result$p.value
        
        # Store the p-value with asterisks if it is significant
        if (p_value < 0.05) {
          p_value_text <- paste0(format(p_value, scientific = FALSE, digits = 3), "*")
        } else {
          p_value_text <- format(p_value, scientific = FALSE, digits = 3)
        }
        
        p_values[histology, col] <- p_value_text
      } else {
        # Set the p-value to NA if there are not enough observations
        p_values[histology, col] <- NA
      }
    }
  }
  
  # Store the p-values for the current cc_ischia_10 category
  p_values_list[[cc_category]] <- p_values
}

# Load the openxlsx library
library(openxlsx)

# Create a new Excel workbook
wb <- createWorkbook()

# Iterate over each cc_ischia_10 category
for (cc_category in names(p_values_list)) {
  # Create a new worksheet for the current cc_ischia_10 category
  addWorksheet(wb, cc_category)

  # Get the p-values for the current cc_ischia_10 category
  p_values <- p_values_list[[cc_category]]

  # Write the p-values to the worksheet
  writeData(wb, sheet = cc_category, x = p_values)

  # Set the row names in the worksheet
  writeData(wb, sheet = cc_category, x = rownames(p_values), startCol = 1, startRow = 2, colNames = FALSE)

  # Add an asterisk formatting for significant p-values
  for (row in 2:(nrow(p_values) + 1)) {
    for (col in 1:ncol(p_values)) {
      p_value <- p_values[row - 1, col]
      if (!is.na(p_value) && grepl("\\*$", p_value)) {
        addStyle(wb, sheet = cc_category, rows = row, cols = col, 
                 style = createStyle(textDecoration = "italic"))
      }
    }
  }
}

# Save the workbook as an Excel file
saveWorkbook(wb, "p_values_results.xlsx", overwrite = TRUE)


########################
IMPORTANT

To Find the P value between 4 groups

pairwise p-values between different histology groups for the CC10 category. The p-values represent the statistical significance of the differences between the groups based on the Wilcoxon rank sum test.

"Here's a breakdown of the output:

For the CC10 category:

The p-value between "Liver_Mets" and "Lymph node" is approximately 4.27e-21.
The p-value between "Liver_Mets" and "PDAC" is approximately 3.33e-13.
The p-value between "Lymph node" and "PDAC" is approximately 3.72e-03.
The p-values are only calculated for pairs of histology groups where there are sufficient observations (at least 3) in both groups. If there are insufficient observations in either group, the p-value is set to NA to indicate that the test could not be performed.

This output provides you with information on the statistical differences between histology groups within the CC10 category. You can use a similar approach to obtain the pairwise p-values for other cc_ischia_10 categories as well."

###
# Iterate over each cc_ischia_10 category
for (cc_category in unique(subset_df$cc_ischia_10)) {
  # Subset the data for the current cc_ischia_10 category
  cc_subset <- subset_df[subset_df$cc_ischia_10 == cc_category, ]
  
  # Create a matrix to store the p-values
  p_values <- matrix(NA, nrow = length(histology_groups), ncol = length(histology_groups),
                     dimnames = list(histology_groups, histology_groups))
  
  # Iterate over each combination of histology groups
  for (i in 1:length(histology_groups)) {
    for (j in 1:length(histology_groups)) {
      if (i != j) {
        # Select the data for the two histology groups
        group1 <- cc_subset[cc_subset$Histology == histology_groups[i], ]
        group2 <- cc_subset[cc_subset$Histology == histology_groups[j], ]
        
        # Check if either group has fewer than 3 observations
        if (nrow(group1) < 3 || nrow(group2) < 3) {
          p_values[i, j] <- NA  # Set p-value to NA if not enough observations
        } else {
          # Perform the Wilcoxon rank sum test for each numeric column
          for (col in colnames(cc_subset)) {
            if (is.numeric(cc_subset[[col]])) {
              p_values[i, j] <- wilcox.test(group1[[col]], group2[[col]])$p.value
            }
          }
        }
      }
    }
  }
  
  # Store the p-values in the list
  p_values_list[[cc_category]] <- p_values
}

# Create a new workbook
wb <- createWorkbook()

# Iterate over CC categories
for (cc_category in names(p_values_list_modified)) {
  # Create a new worksheet for the CC category
  addWorksheet(wb, sheetName = cc_category)
  
  # Get the p-values matrix for the CC category
  p_values_matrix <- p_values_list_modified[[cc_category]]
  
  # Set the column names and row names of the matrix
  col_names <- colnames(p_values_matrix)
  row_names <- rownames(p_values_matrix)
  
  # Create a data frame to store the p-values
  p_values_df <- data.frame(matrix(NA, nrow = length(row_names), ncol = length(col_names)))
  
  # Set the row names and column names of the data frame
  rownames(p_values_df) <- row_names
  colnames(p_values_df) <- col_names
  
  # Copy the p-values from the matrix to the data frame
  for (i in 1:length(row_names)) {
    for (j in 1:length(col_names)) {
      p_values_df[i, j] <- p_values_matrix[i, j]
    }
  }
  
  # Iterate over rows and columns to add asterisks for significance
  for (i in 1:nrow(p_values_df)) {
    for (j in 1:ncol(p_values_df)) {
      p_value <- p_values_df[i, j]
      if (!is.na(p_value)) {
        if (p_value < 0.001) {
          p_values_df[i, j] <- paste0(format(p_value, scientific = FALSE), "***")
        } else if (p_value < 0.01) {
          p_values_df[i, j] <- paste0(format(p_value, scientific = FALSE), "**")
        } else if (p_value < 0.05) {
          p_values_df[i, j] <- paste0(format(p_value, scientific = FALSE), "*")
        }
      }
    }
  }
  
  # Write the data frame to the worksheet, including row names
  writeData(wb, sheet = cc_category, p_values_df, startRow = 1, startCol = 1, row.names = TRUE)
}

# Save the workbook to a file
saveWorkbook(wb, "p_values.xlsx")

######################




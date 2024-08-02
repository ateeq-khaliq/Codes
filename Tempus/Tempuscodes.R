Created : Ateeq
on:02-13-23
Project: Tempus Spatial and our pdac samples analysis
####
# part1
library(Seurat)

# Reading all the 68 files from tempus

# Set the directory where your data is stored
data_dir <- "/Volumes/G-DRIVE_PRO/Tempus_data/"

# Get a list of folders in the directory that end with "b"
folders <- list.files(data_dir, pattern = "b$", full.names = TRUE)

# Create an empty list to store the Seurat objects
seurat_objects <- list()

# Loop over the folders and read in the data for each one
for (folder in folders) {
  cat("Processing folder:", folder, "\n")
  
  # Read in the filtered_feature_bc_matrix.h5 file
  expr.url <- file.path(folder, "filtered_feature_bc_matrix.h5")
  expr.data <- Read10X_h5(filename =  expr.url)
  
  # Create a Seurat object
  seurat_obj <- CreateSeuratObject(counts = expr.data, project = basename(folder), 
                                    assay = 'Spatial', min.cells = 3, min.features = 200)
  
  # Read in the image data
  img.dir <- file.path(folder, "spatial")
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

# part 2
# reading Our Spatial samples only 13 samples


# Set the directory where your data is stored
data_dir <- "/Users/akhaliq/Desktop/panc_samples/"

# Get a list of folders in the directory that end with "b"
folders <- list.files(data_dir, full.names = TRUE)

# Create an empty list to store the Seurat objects
seurat_objects <- list()

# Loop over the folders and read in the data for each one
for (folder in folders) {
  cat("Processing folder:", folder, "\n")
  
  # Read in the filtered_feature_bc_matrix.h5 file
  expr.url <- file.path(folder, "outs/filtered_feature_bc_matrix.h5")
  expr.data <- Read10X_h5(filename =  expr.url)
  
  # Create a Seurat object
  seurat_obj <- CreateSeuratObject(counts = expr.data, project = basename(folder), 
                                    assay = 'Spatial', min.cells = 3, min.features = 200)
  
  # Read in the image data
  img.dir <- file.path(folder, "outs/spatial")
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
data_dir <- "/Volumes/G-DRIVE_PRO/Tempus_data/"

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


# Set the directory where your data is stored
data_dir <- "/Users/akhaliq/Desktop/panc_samples/"

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

pdac_all <- merge(PDACP_9, c(PDACP_8, PDACP_5,PDACP_4,PDACP_3,PDACP_13,PDACP_12,PDACP_11,PDACP_10,PDACP_1,PDACNP_2,PDACNP_10,PDACNP_1, MV068b, MV067b, MV066b, MV065b, MV064b, MV063b, MV062b, MV061b, MV060b, MV059b, MV058b, MV057b, MV056b, MV055b, MV054b, MV053b, MV052b, MV051b, MV050b, MV049b, MV048b, MV047b, MV046b, MV045b, MV044b, MV043b, MV042b, MV041b, MV040b, MV039b, MV038b, MV037b, MV036b, MV035b, MV034b, MV033b, MV032b, MV031b, MV030b, MV029b, MV027b, MV026b, MV025b, MV024b, MV023b, MV022b, MV021b, MV020b, MV019b, MV018b, MV017b, MV016b, MV015b, MV014b, MV013b, MV012b, MV011b, MV010b, MV009b, MV008b, MV007b, MV006b, MV005b, MV004b, MV003b, MV002b, MV001b))


# part 5
# remove all Except pdac_all

rm(list=setdiff(ls(),"pdac_all"))
gc()


##### This is not needed Part 6 
# part 6

# Add images to the merged object (pdac_all)
# Get list of all directories in image directory
#Add Tempus images
image_dir <- "/Volumes/G-DRIVE_PRO/Tempus_data"
dir_list <- list.files(path = image_dir, full.names = TRUE, recursive = FALSE)

# Filter the list to only include directories that end with "b"
dir_list <- dir_list[grep(pattern = "b$", x = dir_list)]

# Loop through each directory and add images to Seurat object
for (dir_path in dir_list) {
  assay_name <- basename(dir_path)
  image_path <- file.path(dir_path, "spatial")
  image <- Seurat::Read10X_Image(image.dir = image_path)
  Seurat::DefaultAssay(object = image) <- 'Spatial'
  image <- image[colnames(x = pdac_all)]
  pdac_all[[assay_name]] <- image
}

# Add our images
image_dir <- "/Users/akhaliq/Desktop/panc_samples"
dir_list <- list.files(path = image_dir, full.names = TRUE, recursive = FALSE)

# Filter the list to only include directories that end with "b"
#dir_list <- dir_list[grep( x = dir_list)]

# Loop through each directory and add images to Seurat object
for (dir_path in dir_list) {
  assay_name <- basename(dir_path)
  image_path <- file.path(dir_path, "outs/spatial")
  image <- Seurat::Read10X_Image(image.dir = image_path)
  Seurat::DefaultAssay(object = image) <- 'Spatial'
  image <- image[colnames(x = pdac_all)]
  pdac_all[[assay_name]] <- image
}

##############

# part 7

# showed an error "Error in FUN(left, right) : non-numeric argument to binary operator"
# The following code should loop through all image names in the pdac_all Seurat object and apply the as.integer() function to the specified fields for each image.


image_names <- c("MV001b", "MV022b", "MV043b", "MV064b", "MV002b", "MV023b", "MV044b", "MV065b",
                  "MV003b", "MV024b", "MV045b", "MV066b", "MV004b", "MV025b", "MV046b", "MV067b",
                  "MV005b", "MV026b", "MV047b", "MV068b", "MV006b", "MV027b", "MV048b", "PDACNP_1",
                  "MV007b", "MV049b", "PDACNP_10", "MV008b", "MV029b", "MV050b", "PDACNP_2",
                  "MV009b", "MV030b", "MV051b", "PDACP_1", "MV010b", "MV031b", "MV052b", "PDACP_10",
                  "MV011b", "MV032b", "MV053b", "PDACP_11", "MV012b", "MV033b", "MV054b", "PDACP_12",
                  "MV013b", "MV034b", "MV055b", "PDACP_13", "MV014b", "MV035b", "MV056b", "PDACP_3",
                  "MV015b", "MV036b", "MV057b", "PDACP_4", "MV016b", "MV037b", "MV058b", "PDACP_5",
                  "MV017b", "MV038b", "MV059b", "PDACP_8", "MV018b", "MV039b", "MV060b", "PDACP_9",
                  "MV019b", "MV040b", "MV061b", "MV020b", "MV041b", "MV062b", "MV021b", "MV042b", "MV063b")

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


###

library(SpotClean)
library(S4Vectors)
spatial_dir: /Volumes/G-DRIVE_PRO/Tempus_data/MV068b/spatial
pdac_raw: /Volumes/G-DRIVE_PRO/Tempus_data/MV068b/raw_feature_bc_matrix

spatial_clean : /Volumes/G-DRIVE_PRO/Tempus_data/spotclean
pdac_raw <- read10xRaw("/Volumes/G-DRIVE_PRO/Tempus_data/MV068b/raw_feature_bc_matrix")
pdac_sample <- read10xSlide(tissue_csv_file=file.path(spatial_dir,
                                       "tissue_positions_list.csv"),
             tissue_img_file = file.path(spatial_dir,
                                       "tissue_lowres_image.png"),
             scale_factor_file = file.path(spatial_dir,
                                       "scalefactors_json.json"))
pdac_obj <- createSlide(count_mat = pdac_raw, 
                          slide_info = pdac_sample)
visualizeSlide(slide_obj = pdac_obj)


# Decontaminate raw data
decont_obj <- spotclean(pdac_obj)

seurat_obj <- convertToSeurat(decont_obj,image_dir = spatial_clean)



library(SpotClean)
library(S4Vectors)

# Set the base directory
base_dir <- "/Volumes/G-DRIVE_PRO/Tempus_data"

# List all the files and directories in the base directory
all_items <- list.files(base_dir, full.names = TRUE)

# Filter only the directories with the suffix "b"
b_dirs <- all_items[file.info(all_items)$isdir & grepl("b$", all_items)]

# Loop through the "b" directories
for (dir in b_dirs) {
  
  # Define the paths for the input and output directories
  spatial_dir <- file.path(dir, "spatial")
  pdac_raw <- file.path(dir, "raw_feature_bc_matrix")
  spatial_clean <- file.path(base_dir, "spotclean")
  
  # Read the raw data
  pdac_raw <- read10xRaw(pdac_raw)
  
  # Read the spatial data
  pdac_sample <- read10xSlide(tissue_csv_file=file.path(spatial_dir, "tissue_positions_list.csv"),
                              tissue_img_file = file.path(spatial_dir, "tissue_lowres_image.png"),
                              scale_factor_file = file.path(spatial_dir, "scalefactors_json.json"))
  
  # Create the slide object
  pdac_obj <- createSlide(count_mat = pdac_raw, slide_info = pdac_sample)
  
  # Visualize the slide
  #visualizeSlide(slide_obj = pdac_obj)
  
  # Decontaminate raw data
  decont_obj <- spotclean(pdac_obj)
  
  # Convert to Seurat object
  seurat_obj <- convertToSeurat(decont_obj, image_dir = spatial_clean)
  
}


## Seurat Analysis for Individual Patients
#
# This script performs individual Seurat analysis for each patient using the provided data. It subsets the data for each
# patient, performs clustering analysis, saves the results in individual folders named after the patient IDs, and
# generates plots and environment files for each patient. Please make sure to update the file paths in the script
# before running it. Refer to the README for more details.

# Loading required packages
library(ISCHIA)
library(robustbase)
library(data.table)
library(ggplot2)
library(Seurat)
library(dplyr)
library(utils)

# Read the main Seurat object
pdac_allres <- readRDS("/Users/akhaliq/Desktop/spatial_analysis/st_new_seurat/all_new_mod/Pdac_allres_most_updated.rds")

# Get the unique patient IDs
patient_ids <- unique(pdac_allres$patient)

# Iterate over each patient
for (patient_id in patient_ids) {
  print(paste("Processing patient:", patient_id))
  
  # Create a directory for the current patient
  patient_dir <- paste0("/Users/akhaliq/Desktop/spatial_analysis/st_new_seurat/all_new_mod/patientwise/", patient_id)
  dir.create(patient_dir, showWarnings = FALSE)
  
  # Subset the Seurat object for the current patient
  patient_seurat <- subset(pdac_allres, subset = (patient == patient_id))
  
  # Get the assay matrix for the current patient
  assay_matrix <- patient_seurat[["rctd_fullfinal"]]@data
  norm_weights <- as.data.frame(t(assay_matrix))
  
  # Deciding about the k and save the elbow plot in patientwise folder
  pdf(file.path(patient_dir, "elbow.pdf"), height = 6, width = 7)
  Composition.cluster.k(norm_weights, 20)
  dev.off()
  
  # Save the Seurat object for the current patient
  saveRDS(patient_seurat, file.path(patient_dir, paste0("pt_", patient_id, ".rds")))
}

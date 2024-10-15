# Step 1 Formatted GCT File from Seurat Object: 

library(Seurat)
library(Matrix)

# 1. Ensure we're using the correct assay
DefaultAssay(tumor_core) <- "SCT"

# 2. Extract the normalized expression matrix
expression_matrix <- GetAssayData(tumor_core, slot = "data", assay = "SCT")

# 3. Convert to dense matrix if it's sparse
if (inherits(expression_matrix, "dgCMatrix")) {
  expression_matrix <- as.matrix(expression_matrix)
}

# 4. Create the GCT file
write_gct <- function(data, file_name) {
  # Open connection to file
  con <- file(file_name, "w")
  
  # Write header
  writeLines(c("#1.2", 
               paste(nrow(data), ncol(data))), 
             con)
  
  # Write column names
  write.table(t(c("NAME", "Description", colnames(data))), 
              con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Write data rows
  for (i in 1:nrow(data)) {
    write.table(t(c(rownames(data)[i], rownames(data)[i], data[i,])), 
                con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  
  # Close connection
  close(con)
}

# 5. Write the GCT file
write_gct(expression_matrix, "tumor_core_cc7_10.gct")

print("Correctly formatted GCT file 'tumor_core_cc7_10.gct' has been created.")

# 6. Print summary for verification
cat("\nSummary of GCT file:\n")
cat("Number of genes:", nrow(expression_matrix), "\n")
cat("Number of spots:", ncol(expression_matrix), "\n")
cat("File size:", file.size("tumor_core_cc7_10.gct") / 1e6, "MB\n")

# 7. Check the first few lines of the GCT file
cat("\nFirst few lines of the GCT file:\n")
system("head -n 10 tumor_core_cc7_10.gct")

# 8. Check for any NA values
na_count <- sum(is.na(expression_matrix))
cat("\nNumber of NA values in the expression matrix:", na_count, "\n")


####

# Step 2 : Create Matching CLS File for Spatial Data

library(Seurat)

# 1. Get the cell barcodes from the expression matrix
barcodes <- colnames(expression_matrix)

# 2. Get the 'treated' information from the metadata
treatment_status <- tumor_core$treated[barcodes]

# 3. Convert treatment status to numeric (0 for No, 1 for Yes)
phenotypes <- as.integer(treatment_status == "Yes")

# 4. Create the CLS file content
cls_content <- c(
  paste(length(phenotypes), "2", "1"),
  "# Untreated Treated",
  paste(phenotypes, collapse = " ")
)

# 5. Write the CLS file
writeLines(cls_content, "tumor_core_cc7_10.cls")

print("CLS file 'tumor_core_cc7_10.cls' has been created.")

# 6. Print summary for verification
cat("\nSummary:\n")
cat("Total samples:", length(phenotypes), "\n")
cat("Treated samples:", sum(phenotypes == 1), "\n")
cat("Untreated samples:", sum(phenotypes == 0), "\n")

# 7. Print the content of the CLS file for verification
cat("\nContent of the CLS file:\n")
cat(readLines("tumor_core_cc7_10.cls"), sep = "\n")

# 8. Check if all samples have a phenotype
cat("\nAll samples have a phenotype assigned:", 
    all(!is.na(phenotypes)), "\n")


###

#Step 3: Create Correct Ranked Gene List for GSEA

library(data.table)
library(dplyr)

# Read the GCT file
gct_file <- "tumor_core_cc7_10.gct"
gct_data <- fread(gct_file, skip = 2)  # Skip the first two lines of the GCT file

# Extract gene names and expression data
gene_names <- gct_data[[1]]
expression_data <- gct_data[, 3:ncol(gct_data), with = FALSE]  # Skip the first two columns (NAME and Description)

# Calculate the mean expression for each gene
mean_expression <- rowMeans(expression_data, na.rm = TRUE)

# Create a data frame with gene names and mean expression
ranked_genes <- data.frame(gene = gene_names, rank = mean_expression)

# Remove any rows with NA values
ranked_genes <- ranked_genes[complete.cases(ranked_genes), ]

# Sort genes by mean expression in descending order
ranked_genes <- ranked_genes[order(-ranked_genes$rank), ]

# Write the ranked gene list to a file WITHOUT column names
write.table(ranked_genes, file = "tumor_core_ranked_cc7_10.rnk", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

print("Ranked gene list has been created: /N/project/akhaliq/gsea/tumor_core_ranked.rnk")

# Print some statistics
cat("Total number of genes:", nrow(ranked_genes), "\n")
cat("Top 5 genes:\n")
print(head(ranked_genes, 5))
cat("Bottom 5 genes:\n")
print(tail(ranked_genes, 5))

###

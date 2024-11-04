
# Function to perform comprehensive QC on duplicate genes
perform_duplicate_genes_qc <- function(count_matrix) {
    library(dplyr)
    library(ggplot2)
    
    # Create output directory for QC
    dir.create("duplicate_genes_qc", showWarnings = FALSE)
    
    # Add Gene_ID column from rownames
    count_matrix$Gene_ID <- rownames(count_matrix)
    
    # 1. Find and analyze duplicates
    duplicates_analysis <- function(count_matrix) {
        # Find duplicated genes
        dups <- duplicated(count_matrix$Gene_Name) | 
                duplicated(count_matrix$Gene_Name, fromLast = TRUE)
        dup_genes <- unique(count_matrix$Gene_Name[dups])
        
        # Create detailed duplicate analysis
        dup_details <- data.frame(
            Gene_Name = character(),
            Gene_ID = character(),
            Gene_Type = character(),
            Mean_Expression = numeric(),
            Max_Expression = numeric(),
            Non_Zero_Samples = numeric(),
            stringsAsFactors = FALSE
        )
        
        # Get sample columns (all numeric columns)
        sample_cols <- setdiff(colnames(count_matrix), 
                             c("Gene_Name", "Gene_Type", "Gene_ID"))
        
        for(gene in dup_genes) {
            # Get all rows for this gene
            gene_rows <- count_matrix[count_matrix$Gene_Name == gene, ]
            
            # Calculate metrics
            gene_metrics <- data.frame(
                Gene_Name = gene_rows$Gene_Name,
                Gene_ID = gene_rows$Gene_ID,
                Gene_Type = gene_rows$Gene_Type,
                Mean_Expression = rowMeans(gene_rows[, sample_cols]),
                Max_Expression = apply(gene_rows[, sample_cols], 1, max),
                Non_Zero_Samples = rowSums(gene_rows[, sample_cols] > 0)
            )
            
            dup_details <- rbind(dup_details, gene_metrics)
        }
        
        # Print duplicate gene summary
        cat(sprintf("\nFound %d genes with duplicates:\n", length(dup_genes)))
        print(table(dup_details$Gene_Type))
        
        return(list(
            details = dup_details,
            summary = data.frame(
                Total_Duplicates = length(dup_genes),
                Unique_Gene_Types = length(unique(dup_details$Gene_Type)),
                Max_Duplicates_Per_Gene = max(table(dup_details$Gene_Name))
            )
        ))
    }
    
    # 2. Generate resolution strategy
    resolve_duplicates <- function(dup_details) {
        resolution_strategy <- data.frame(
            Gene_Name = character(),
            Selected_ID = character(),
            Strategy = character(),
            Mean_Expression = numeric(),
            stringsAsFactors = FALSE
        )
        
        for(gene in unique(dup_details$Gene_Name)) {
            gene_rows <- dup_details[dup_details$Gene_Name == gene, ]
            
            # Strategy 1: Keep protein coding if available
            if("protein_coding" %in% gene_rows$Gene_Type) {
                protein_coding_rows <- gene_rows[gene_rows$Gene_Type == "protein_coding", ]
                # If multiple protein coding, take highest expression
                selected_row <- protein_coding_rows[which.max(protein_coding_rows$Mean_Expression), ]
                strategy <- "protein_coding"
            }
            # Strategy 2: Highest expression
            else {
                selected_row <- gene_rows[which.max(gene_rows$Mean_Expression), ]
                strategy <- "highest_expression"
            }
            
            resolution_strategy <- rbind(resolution_strategy,
                                      data.frame(
                                          Gene_Name = gene,
                                          Selected_ID = selected_row$Gene_ID,
                                          Strategy = strategy,
                                          Mean_Expression = selected_row$Mean_Expression
                                      ))
        }
        
        return(resolution_strategy)
    }
    
    # 3. Perform analysis
    qc_results <- duplicates_analysis(count_matrix)
    resolution_strategy <- resolve_duplicates(qc_results$details)
    
    # 4. Generate detailed report
    sink("duplicate_genes_qc/duplicate_genes_report.txt")
    cat("Duplicate Genes QC Report\n")
    cat("========================\n\n")
    
    cat("Summary Statistics:\n")
    cat(sprintf("Total number of duplicate genes: %d\n", qc_results$summary$Total_Duplicates))
    cat(sprintf("Number of unique gene types involved: %d\n", qc_results$summary$Unique_Gene_Types))
    cat(sprintf("Maximum duplicates for a single gene: %d\n\n", qc_results$summary$Max_Duplicates_Per_Gene))
    
    cat("Gene Type Distribution in Duplicates:\n")
    print(table(qc_results$details$Gene_Type))
    cat("\n")
    
    cat("Resolution Strategy Summary:\n")
    cat(sprintf("Protein coding priority: %d genes\n", sum(resolution_strategy$Strategy == "protein_coding")))
    cat(sprintf("Highest expression priority: %d genes\n", sum(resolution_strategy$Strategy == "highest_expression")))
    sink()
    
    # 5. Save detailed results
    write.csv(qc_results$details, 
              "duplicate_genes_qc/duplicate_genes_details.csv", 
              row.names = FALSE)
    write.csv(resolution_strategy, 
              "duplicate_genes_qc/resolution_strategy.csv", 
              row.names = FALSE)
    
    # 6. Function to apply resolution strategy
    apply_resolution <- function(original_matrix, resolution_strategy) {
        # Keep only the selected gene IDs
        resolved_matrix <- original_matrix[original_matrix$Gene_ID %in% resolution_strategy$Selected_ID, ]
        return(resolved_matrix)
    }
    
    # 7. Create visualizations
    pdf("duplicate_genes_qc/duplicate_gene_plots.pdf")
    
    # Plot 1: Gene type distribution
    barplot(table(qc_results$details$Gene_Type),
            main="Gene Types in Duplicates",
            las=2)
    
    # Plot 2: Expression distribution of duplicates
    hist(qc_results$details$Mean_Expression,
         main="Expression Distribution of Duplicate Genes",
         xlab="Mean Expression",
         breaks=20)
    
    dev.off()
    
    return(list(
        qc_results = qc_results,
        resolution_strategy = resolution_strategy,
        apply_resolution = apply_resolution
    ))
}
global_counts <- read.csv("/Users/akhaliq/Desktop/asif/bulk_counts/run_20241103_093957/counts/counts_protein_coding.csv",row.names = 1, header = TRUE)

# Usage example:
# qc_results <- perform_duplicate_genes_qc(global_counts)
# resolved_counts <- qc_results$apply_resolution(global_counts, qc_results$resolution_strategy)

# Run the QC analysis
qc_results <- perform_duplicate_genes_qc(global_counts)

# Look at the duplicate genes summary
View(qc_results$qc_results$details)

# Check the resolution strategy
View(qc_results$resolution_strategy)

# Create the resolved matrix (removes duplicates according to the strategy)
resolved_counts <- qc_results$apply_resolution(global_counts, qc_results$resolution_strategy)

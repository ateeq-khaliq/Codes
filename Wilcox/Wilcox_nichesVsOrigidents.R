# Merge the norm_weights dataframe with the meta dataframe based on the spot IDs
merged_df <- merge(norm_weights, meta, by = "row.names", all.x = TRUE)

# Subset the merged dataframe to include only the relevant columns
subset_df <- subset(merged_df, select = c("orig.ident", "cc_ischia_12", colnames(norm_weights)))

# Create an empty list to store the enrichment results
enrichment_results <- list()

# Iterate over each histology category
for (orig.ident in unique(subset_df$orig.ident)) {
  # Subset the data for the current orig.ident category
  orig.ident_subset <- subset_df[subset_df$orig.ident == orig.ident, ]
  
  # Create a matrix to store the p-values for each cc_ischia_12 category
  p_values <- matrix(NA, nrow = length(unique(orig.ident_subset$cc_ischia_12)), 
                     ncol = ncol(norm_weights), 
                     dimnames = list(unique(orig.ident_subset$cc_ischia_12), colnames(norm_weights)))
  
  # Iterate over each cc_ischia_12 category
  for (cc_category in unique(orig.ident_subset$cc_ischia_12)) {
    # Subset the data for the current cc_ischia_12 category
    cc_subset <- orig.ident_subset[orig.ident_subset$cc_ischia_12 == cc_category, ]
    
    # Iterate over each proportion column and perform Wilcoxon test
    for (col in colnames(norm_weights)) {
      group1 <- cc_subset[, col]
      group2 <- subset_df[subset_df$orig.ident != orig.ident & subset_df$cc_ischia_12 == cc_category, col]
      
      # Check if both groups have enough observations for the test
      if (length(group1) > 1 && length(group2) > 1) {
        # Perform Wilcoxon test and store the p-value
        p_value <- wilcox.test(group1, group2)$p.value
        p_values[cc_category, col] <- p_value
      } else {
        # Set the p-value to NA if there are not enough observations
        p_values[cc_category, col] <- NA
      }
    }
  }
  
  # Store the enrichment results for the current orig.ident category
  enrichment_results[[orig.ident]] <- p_values
}

library(ggplot2)
library(patchwork)

# Create a list to store the plots
plots <- list()

# Create box plots for each orig.ident category
for (orig_ident in names(enrichment_results)) {
  # Get the enrichment results for the current orig.ident category
  orig_ident_results <- enrichment_results[[orig_ident]]
  
  # Create a data frame with the p-values for each cc_ischia_12 category
  p_values_df <- data.frame(cc_ischia_12 = rownames(orig_ident_results),
                            p_value = as.numeric(orig_ident_results))
  
  # Remove rows with NA values
  p_values_df <- na.omit(p_values_df)
  
  # Create a new plot for each orig.ident category
  boxplot_plot <- ggplot(p_values_df, aes(x = cc_ischia_12, y = p_value)) +
    geom_boxplot(fill = "skyblue", color = "black", alpha = 0.8, outlier.shape = NA) +  # Modify the fill color here
    labs(title = orig_ident, x = "Composition Clusters", y = "p-value") +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  
  # Add the plot to the list of plots
  plots[[orig_ident]] <- boxplot_plot
}

# Determine the number of pages
num_plots <- length(plots)
plots_per_page <- 10
num_pages <- ceiling(num_plots / plots_per_page)

# Create a list to store the grid plots
grid_plots_list <- list()

# Create box plot grids for each page
for (page in 1:num_pages) {
  start_index <- (page - 1) * plots_per_page + 1
  end_index <- min(start_index + plots_per_page - 1, num_plots)
  
  page_plots <- plots[start_index:end_index]
  
  # Arrange the plots using patchwork
  page_grid <- wrap_plots(page_plots, ncol = 2)
  
  # Add the grid plot to the list
  grid_plots_list[[page]] <- page_grid
}

# Combine all the pages into a single grid layout
grid_plots_all <- wrap_plots(grid_plots_list, nrow = 1)

# Save all the plots as a single PDF file
ggsave("Wilcoxon_CCVsorig.ident_AllPages.pdf", grid_plots_all, width = 14, height = 12)

# Print the success message
cat("All box plots saved as 'Wilcoxon_CCVsorig.ident_AllPages.pdf'\n")


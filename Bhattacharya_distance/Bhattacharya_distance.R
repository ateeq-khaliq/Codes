# Install and load necessary libraries
install.packages("data.table")
install.packages("ggplot2")
install.packages("reshape2")

library(data.table)
library(ggplot2)
library(reshape2)

# Function to calculate Bhattacharyya distance between two distributions
bhattacharyya_distance <- function(x, y) {
  mean_x <- mean(x)
  mean_y <- mean(y)
  cov_x <- var(x)
  cov_y <- var(y)
  bc <- 0.25 * log(0.25 * (cov_x / cov_y + cov_y / cov_x + 2)) + 
        0.25 * ((mean_x - mean_y)^2 / (cov_x + cov_y))
  return(bc)
}

# Read the data
file_path <- '/Users/akhaliq/Desktop/mouse/Caf_metaprograms.csv'
data <- fread(file_path)

# Remove the first column as it is non-numeric
numeric_data <- data[, -1, with = FALSE]

# Calculate Bhattacharyya distance matrix
n <- ncol(numeric_data)
dist_matrix <- matrix(0, n, n)
colnames(dist_matrix) <- colnames(numeric_data)
rownames(dist_matrix) <- colnames(numeric_data)

for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    dist_matrix[i, j] <- bhattacharyya_distance(numeric_data[[i]], numeric_data[[j]])
    dist_matrix[j, i] <- dist_matrix[i, j]
  }
}

# Melt the distance matrix for ggplot
dist_melted <- melt(dist_matrix)

# Plot the Bhattacharyya distance matrix as a heatmap
ggplot(dist_melted, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  labs(title = "Bhattacharyya Distance Matrix",
       x = "Variables",
       y = "Variables",
       fill = "Distance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the distance matrix to a CSV file
write.csv(dist_matrix, '/Users/akhaliq/Desktop/mouse/Bhattacharyya_Distance_Matrix.csv')



###


# Create a PDF file to save the plots
pdf("Comparative_Analysis.pdf", width = 12, height = 10)

# Plot Bhattacharyya distance heatmap
pheatmap(bhattacharyya_matrix, 
         clustering_method = "complete",
         color = colorRampPalette(c("white", "red"))(200),
         display_numbers = TRUE, 
         number_color = "black",
         fontsize = 10,
         fontsize_number = 8,
         border_color = NA,
         main = "Bhattacharyya Distance Matrix")

# Plot Pearson correlation heatmap
pheatmap(pearson_matrix, 
         clustering_method = "complete",
         color = colorRampPalette(c("blue", "white", "red"))(200),
         display_numbers = TRUE, 
         number_color = "black",
         fontsize = 10,
         fontsize_number = 8,
         border_color = NA,
         main = "Pearson Correlation Matrix")

# Scatter plot comparing Bhattacharyya distance and Pearson correlation
ggplot(comparison_df, aes(x = value.x, y = value.y)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = 'lm', color = 'blue') +
  theme_minimal() +
  labs(title = "Comparison between Bhattacharyya Distance and Pearson Correlation",
       x = "Bhattacharyya Distance",
       y = "Pearson Correlation")

# Close the PDF device
dev.off()

# Save the matrices to CSV files
write.csv(bhattacharyya_matrix, '/path/to/your/Bhattacharyya_Distance_Matrix.csv')
write.csv(pearson_matrix, '/path/to/your/Pearson_Correlation_Matrix.csv')




library(ggplot2)
library(gridExtra)
library(factoextra)
library(cluster)
library(fpc)


set.seed(123)
norm_weights <- read.csv("/Users/akhaliq/Desktop/spatial_analysis/st_new_seurat/all_new_mod/norm_weights.csv",header=T,sep=',',row.names=1)

library(ggplot2)
library(gridExtra)
library(cluster)
library(fpc)

# Set the seed for reproducibility
set.seed(123)

# Function to compute total within-cluster sum of squares
wss <- function(k) {
  kmeans(norm_weights, k, nstart = 10)$tot.withinss
}

# Compute and plot wss for k = 1 to k = 20
k.values <- 1:20
wss_values <- sapply(k.values, wss)

# Find the optimal K using the elbow method
diff_wss <- diff(wss_values)
curvature <- diff(diff_wss)
elbow_index <- which.max(curvature) + 1
elbow_x <- k.values[elbow_index]
elbow_y <- wss_values[elbow_index]

# Plot the wss values
plot_wss <- ggplot() +
  geom_line(aes(x = k.values, y = wss_values), color = "blue", size = 1) +
  geom_point(aes(x = k.values, y = wss_values), color = "blue", size = 3) +
  geom_vline(xintercept = elbow_x, color = "red", linetype = "dashed") +
  geom_hline(yintercept = elbow_y, color = "red", linetype = "dashed") +
  labs(x = "Number of Clusters (K)", y = "Total Within-Cluster Sum of Squares",
       title = "Elbow Method: Determining Optimal K") +
  theme_minimal()

# Compute Calinski-Harabasz Index for k = 2 to k = 20
calinski_harabasz_index <- function(k) {
  km.res <- kmeans(norm_weights, k, nstart = 10)
  cluster_stats <- cluster.stats(dist(norm_weights), km.res$cluster)
  cluster_stats$ch
}

calinski_harabasz_values <- sapply(k.values, calinski_harabasz_index)

# Find the optimal K using the Calinski-Harabasz Index
optimal_k_calinski <- which.max(calinski_harabasz_values)

# Plot the Calinski-Harabasz Index
plot_calinski_harabasz <- ggplot() +
  geom_line(aes(x = k.values, y = calinski_harabasz_values), color = "red", size = 1) +
  geom_point(aes(x = k.values, y = calinski_harabasz_values), color = "red", size = 3) +
  geom_vline(xintercept = optimal_k_calinski, color = "blue", linetype = "dashed") +
  labs(x = "Number of Clusters (K)", y = "Calinski-Harabasz Index",
       title = "Calinski-Harabasz Index: Determining Optimal K") +
  theme_minimal()

# Compute and plot Gap statistic for k = 1 to k = 20
gap_stat <- function(k) {
  km.res <- kmeans(norm_weights, k, nstart = 10)
  cluster_stats <- cluster.stats(dist(norm_weights), km.res$cluster)
  cluster_stats$gap
}

gap_stat_values <- sapply(k.values, gap_stat)

# Find the optimal K using the Gap statistic
optimal_k_gap <- which.max(gap_stat_values)

# Plot the Gap statistic
plot_gap_stat <- ggplot() +
  geom_line(aes(x = k.values, y = gap_stat_values), color = "green", size = 1) +
  geom_point(aes(x = k.values, y = gap_stat_values), color = "green", size = 3) +
  geom_vline(xintercept = optimal_k_gap, color = "blue", linetype = "dashed") +
  labs(x = "Number of Clusters (K)", y = "Gap Statistic",
       title = "Gap Statistic: Determining Optimal K") +
  theme_minimal()

# Combine all plots into one grid
grid <- grid.arrange(plot_wss, plot_calinski_harabasz, plot_gap_stat, nrow = 1)

# Save all plots in one PDF
pdf("/Users/akhaliq/Desktop/spatial_analysis/st_new_seurat/all_new_mod/ischia/cluster_validation_plots.pdf")
print(grid)
dev.off()

# Print the optimal K values for each method
optimal_k_values <- c(elbow_x, optimal_k_calinski, optimal_k_gap)
cat("Optimal K values:\n")
cat("Elbow Method:", elbow_x, "\n")
cat("Calinski-Harabasz Index:", optimal_k_calinski, "\n")
cat("Gap Statistic:", optimal_k_gap, "\n")


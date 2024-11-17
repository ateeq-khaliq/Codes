# Part 3: Usage Example
# Save this as "3_run_analysis.R"

# Source the functions
source("1_statistical_analysis.R")
source("2_visualization.R")

# Create proportion data from your Seurat object
proportion_data <- mycaf@meta.data %>%
  group_by(sample_id, treatment_group) %>%
  count(assigned_MP) %>%
  group_by(sample_id, treatment_group) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

# Run statistical analysis
results <- run_statistical_analysis(
  data = proportion_data,
  cell_type_col = "assigned_MP",
  treatment_col = "treatment_group",
  proportion_col = "proportion"
)

# Save statistical results
write.csv(results, "statistical_results.csv", row.names = FALSE)

# Create and save all plots
plots <- create_statistical_plots(results, proportion_data)
save_plots(plots, prefix = "mp_analysis")

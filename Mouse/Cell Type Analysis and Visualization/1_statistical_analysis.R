
# Part 1: Statistical Analysis
# Save this as "1_statistical_analysis.R"

# Load required packages
library(dplyr)
library(tidyr)

# Function to perform statistical test with exact p-values
perform_stat_test <- function(data, mp_type, group1, group2) {
  test_data <- data %>%
    filter(
      assigned_MP == mp_type,
      treatment_group %in% c(group1, group2)
    )
  
  # Get group values
  g1_vals <- test_data$proportion[test_data$treatment_group == group1]
  g2_vals <- test_data$proportion[test_data$treatment_group == group2]
  
  # Perform both Wilcoxon and t-test
  wilcox_test <- wilcox.test(
    g1_vals, g2_vals,
    exact = TRUE,
    alternative = "two.sided"
  )
  
  t_test <- t.test(
    g1_vals, g2_vals,
    var.equal = FALSE
  )
  
  # Calculate means and SDs
  g1_mean <- mean(g1_vals)
  g2_mean <- mean(g2_vals)
  g1_sd <- sd(g1_vals)
  g2_sd <- sd(g2_vals)
  
  # Calculate effect size (Cohen's d)
  pooled_sd <- sqrt(((length(g1_vals)-1)*g1_sd^2 + (length(g2_vals)-1)*g2_sd^2) /
                     (length(g1_vals) + length(g2_vals) - 2))
  cohens_d <- abs(g1_mean - g2_mean) / pooled_sd
  
  return(data.frame(
    assigned_MP = mp_type,
    group1 = group1,
    group2 = group2,
    group1_mean = g1_mean,
    group1_sd = g1_sd,
    group2_mean = g2_mean,
    group2_sd = g2_sd,
    mean_difference = g2_mean - g1_mean,
    fold_change = g2_mean / g1_mean,
    wilcox_p = wilcox_test$p.value,
    t_test_p = t_test$p.value,
    cohens_d = cohens_d
  ))
}

# Main analysis function
run_statistical_analysis <- function(data, cell_type_col = "assigned_MP", 
                                   treatment_col = "treatment_group", 
                                   proportion_col = "proportion") {
  # Rename columns to match function expectations
  data <- data %>%
    rename(
      assigned_MP = !!cell_type_col,
      treatment_group = !!treatment_col,
      proportion = !!proportion_col
    )
  
  # Get unique values
  mp_types <- unique(data$assigned_MP)
  treatment_groups <- unique(data$treatment_group)
  
  # Run all pairwise comparisons
  all_stats <- list()
  for(mp in mp_types) {
    for(i in 1:(length(treatment_groups)-1)) {
      for(j in (i+1):length(treatment_groups)) {
        all_stats[[length(all_stats) + 1]] <- perform_stat_test(
          data, mp, 
          treatment_groups[i], 
          treatment_groups[j]
        )
      }
    }
  }
  
  # Combine and process results
  comprehensive_stats <- bind_rows(all_stats) %>%
    group_by(assigned_MP) %>%
    mutate(
      wilcox_p_adj = p.adjust(wilcox_p, method = "BH"),
      t_test_p_adj = p.adjust(t_test_p, method = "BH"),
      min_p = pmin(wilcox_p, t_test_p),
      min_p_adj = p.adjust(min_p, method = "BH"),
      significance = case_when(
        min_p_adj < 0.001 ~ "***",
        min_p_adj < 0.01 ~ "**",
        min_p_adj < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    ) %>%
    ungroup() %>%
    arrange(assigned_MP, min_p) %>%
    mutate(across(where(is.numeric), ~round(., 4)))
  
  return(comprehensive_stats)
}

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

all <- readRDS("/Users/akhaliq/Desktop/mouse/Cleaned_data/All_celltypes_merged.rds")

# Set cell type order
cell_type_order <- c("B Cells", "CAFs", "DC", "Endothelial Cells", 
                     "Epithelial cells", "Myeloids", "T cells", "Tumor Cells")

# Calculate raw counts
cell_counts <- data.frame(
  CellType = all$new_celltype,
  Treatment = all$treatment_group
) %>%
  group_by(Treatment, CellType) %>%
  summarise(
    count = n(),
    .groups = 'drop'
  )

# Calculate total cells per treatment
treatment_totals <- cell_counts %>%
  group_by(Treatment) %>%
  summarise(
    total = sum(count),
    .groups = 'drop'
  )

# Calculate proportions
proportions <- cell_counts %>%
  complete(Treatment, CellType) %>%  # Ensure all combinations exist
  replace_na(list(count = 0)) %>%    # Replace NA counts with 0
  left_join(treatment_totals, by = "Treatment") %>%
  mutate(proportion = (count / total) * 100)

# Save raw counts and proportions
write.csv(cell_counts %>%
            pivot_wider(names_from = Treatment, 
                       values_from = count, 
                       values_fill = 0), 
          "cell_counts.csv", row.names = FALSE)

write.csv(proportions %>%
            select(CellType, Treatment, proportion) %>%
            pivot_wider(names_from = Treatment, 
                       values_from = proportion), 
          "cell_proportions.csv", row.names = FALSE)

# Calculate DMP
dmp_data <- proportions %>%
  select(CellType, Treatment, proportion) %>%
  pivot_wider(
    id_cols = CellType,
    names_from = Treatment,
    values_from = proportion,
    values_fill = 0
  ) %>%
  mutate(
    Folfiri = Folfiri - Control,
    Proagio = Proagio - Control,
    `Proagio+Folfiri` = `Proagio+Folfiri` - Control
  ) %>%
  select(CellType, Folfiri, Proagio, `Proagio+Folfiri`)

# Statistical testing
stat_results <- lapply(unique(all$new_celltype), function(cell_type) {
  lapply(c("Folfiri", "Proagio", "Proagio+Folfiri"), function(treat) {
    control_data <- c(
      sum(all$new_celltype[all$treatment_group == "Control"] == cell_type),
      sum(all$treatment_group == "Control")
    )
    treat_data <- c(
      sum(all$new_celltype[all$treatment_group == treat] == cell_type),
      sum(all$treatment_group == treat)
    )
    
    test <- fisher.test(matrix(c(
      control_data[1], control_data[2] - control_data[1],
      treat_data[1], treat_data[2] - treat_data[1]
    ), nrow = 2))
    
    data.frame(
      CellType = cell_type,
      Treatment = treat,
      p.value = test$p.value
    )
  }) %>% bind_rows()
}) %>% bind_rows()

# Save statistical results
write.csv(stat_results, "statistical_results.csv", row.names = FALSE)

# Create heatmap data
heatmap_data <- dmp_data %>%
  pivot_longer(
    -CellType,
    names_to = "Treatment",
    values_to = "DMP"
  ) %>%
  mutate(
    Treatment = factor(Treatment, levels = c("Proagio+Folfiri", "Proagio", "Folfiri")),
    CellType = factor(CellType, levels = cell_type_order)
  )

# Save DMP values
write.csv(dmp_data, "dmp_values.csv", row.names = FALSE)

# Create proportion bar plot
prop_bars <- proportions %>%
  group_by(CellType) %>%
  summarise(proportion = mean(proportion))

p1 <- ggplot(prop_bars, aes(x = CellType, y = proportion)) +
  geom_bar(stat = "identity", fill = "grey70") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    plot.margin = margin(b = -10)
  ) +
  labs(y = "Proportion(%)") +
  scale_y_continuous(limits = c(0, max(prop_bars$proportion) * 1.2))

# Create heatmap
p2 <- ggplot(heatmap_data, aes(x = CellType, y = Treatment)) +
  geom_tile(aes(fill = DMP), color = "white", linewidth = 0.5) +
  geom_text(data = stat_results %>% filter(p.value < 0.001),
            aes(x = CellType, y = Treatment),
            label = "*",
            size = 5,
            color = "black") +
  scale_fill_gradient2(
    low = "#4575B4",
    mid = "white",
    high = "#D73027",
    midpoint = 0,
    limits = c(-0.2, 0.2),
    name = "DMP",
    breaks = seq(-0.2, 0.2, by = 0.1),
    oob = scales::squish
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "right"
  ) +
  labs(x = "Cell Type", y = "")

# Combine plots
combined_plot <- plot_grid(
  p1, p2,
  ncol = 1,
  align = "v",
  rel_heights = c(0.3, 0.7)
)

# Save plot in different formats
pdf("treatment_changes_plot.pdf", width = 12, height = 8)
print(combined_plot)
dev.off()

png("treatment_changes_plot.png", width = 12, height = 8, units = "in", res = 300)
print(combined_plot)
dev.off()

# Create summary results table
summary_results <- heatmap_data %>%
  left_join(stat_results, by = c("CellType", "Treatment")) %>%
  mutate(
    Significant = ifelse(p.value < 0.001, "Yes", "No"),
    Effect = case_when(
      DMP > 0.1 ~ "Strong Increase",
      DMP > 0 ~ "Mild Increase",
      DMP < -0.1 ~ "Strong Decrease",
      DMP < 0 ~ "Mild Decrease",
      TRUE ~ "No Change"
    )
  )

# Save summary results
write.csv(summary_results, "summary_results.csv", row.names = FALSE)

# Print summary of key findings
print("Summary of key findings:")
summary_results %>%
  filter(abs(DMP) > 0.1) %>%  # Show only substantial changes
  arrange(desc(abs(DMP))) %>%
  select(CellType, Treatment, DMP, Significant, Effect) %>%
  print(n = 100)


  ###


  # Statistical Methods Used
1. **Fisher's Exact Test**
   - Used to test for significant differences in cell type proportions between treatment and control groups
   - Chosen because it's appropriate for comparing proportions and handles small counts well
   - Significance threshold: p < 0.001 (marked with asterisks *)
   - More conservative than chi-square test, better for uneven group sizes

2. **DMP (Difference in Mean Proportion) Calculation**
   - Formula: Treatment_proportion - Control_proportion
   - Measures absolute change in percentage points
   - Range in visualization: -0.2 to 0.2 (or -20% to +20%)

# Key Findings by Cell Type

## Tumor Cells (Major Findings)
- Strongest reductions observed:
  - Folfiri: -16.1% (p < 0.001)
  - Proagio+Folfiri: -13.7% (p < 0.001)
- Paradoxical increase with Proagio alone: +2.36% (p < 0.001)
- Suggests Folfiri-containing treatments are more effective at reducing tumor burden

## Cancer-Associated Fibroblasts (CAFs)
- Divergent responses:
  - Increased with Folfiri: +7.08% (p < 0.001)
  - Decreased with Proagio: -5.06% (p < 0.001)
  - Decreased with Proagio+Folfiri: -3.27% (p < 0.001)
- Suggests Proagio may target CAF population

## Immune Cell Populations
1. Myeloids
   - Increased across all treatments:
     - Proagio+Folfiri: +9.05% (p < 0.001)
     - Folfiri: +6.02% (p < 0.001)
     - Proagio: +2.53% (p < 0.001)
   - Suggests enhanced immune cell recruitment

2. T cells
   - Variable responses:
     - Strong increase with Proagio+Folfiri: +5.10% (p < 0.001)
     - Moderate increase with Folfiri: +0.955% (p < 0.001)
     - Slight decrease with Proagio: -0.102% (not significant)

3. B Cells
   - Consistent increases:
     - Folfiri: +1.79% (p < 0.001)
     - Proagio+Folfiri: +1.31% (p < 0.001)
     - Proagio: +0.145% (p < 0.001)

## Minor Cell Populations
1. Dendritic Cells (DC)
   - Small but significant increases:
     - Proagio+Folfiri: +0.355% (p < 0.001)
     - Folfiri: +0.179% (p < 0.001)
     - Proagio: +0.0688% (p < 0.001)

2. Endothelial Cells
   - Modest changes:
     - Proagio+Folfiri: +1.08% (p < 0.001)
     - Folfiri: +0.112% (p < 0.001)
     - Proagio: no change

# Treatment-Specific Effects

## Proagio+Folfiri
1. Most consistent significant changes across cell types
2. Strong reduction in tumor cells (-13.7%)
3. Highest increase in myeloids (+9.05%)
4. Best T cell response (+5.10%)

## Folfiri Alone
1. Strongest tumor cell reduction (-16.1%)
2. Highest CAF increase (+7.08%)
3. Strong myeloid recruitment (+6.02%)

## Proagio Alone
1. Unique CAF reduction (-5.06%)
2. Modest myeloid increase (+2.53%)
3. Less effective for tumor cell reduction

# Clinical Implications

1. Combination Benefit
   - Proagio+Folfiri shows balanced effects:
     - Strong tumor reduction
     - Enhanced immune cell recruitment
     - Moderate CAF reduction

2. Treatment Selection
   - Folfiri: Best for direct tumor reduction
   - Proagio: Better for targeting CAFs
   - Combination: Most balanced approach

3. Immune Response
   - All treatments enhance immune cell presence
   - Strongest immune response with combination therapy
   - Suggests potential for immunotherapy combinations

# Limitations and Considerations

1. **Statistical Analysis**
   - Very small p-values due to large sample sizes
   - Biological significance may differ from statistical significance
   - Effect sizes should be considered alongside p-values

2. **Cell Type Interactions**
   - Analysis doesn't capture dynamic interactions
   - Temporal changes not captured
   - Spatial relationships not considered


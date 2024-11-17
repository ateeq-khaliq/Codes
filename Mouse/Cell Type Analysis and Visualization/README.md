
# **Cell Type Analysis and Visualization**

This repository contains R scripts to perform statistical analysis and visualization of cell type proportions across treatment groups. The workflow is designed for data derived from single-cell or spatial transcriptomics experiments, particularly leveraging metadata from Seurat objects.

## **Overview**

### **Scripts:**
1. **`1_statistical_analysis.R`**  
   Performs statistical comparisons of cell type proportions between treatment groups. Outputs include p-values (from Wilcoxon and t-tests), adjusted p-values, fold changes, mean differences, and Cohen's d (effect size).

2. **`2_visualization.R`**  
   Creates visualizations to represent statistical results, including boxplots with significance annotations and other customizable plots.

3. **`3_run_analysis.R`**  
   Demonstrates how to use the statistical and visualization functions on a real dataset. Prepares input data, runs statistical analysis, and generates plots for reporting.

---

## **Requirements**

### **R Version:**
- The scripts require R (≥ 4.0.0).

### **Required R Packages:**
- **Data Manipulation**:  
  - `dplyr`
  - `tidyr`
- **Statistical Analysis**:  
  - Base R functions (`wilcox.test`, `t.test`, etc.)
- **Visualization**:  
  - `ggplot2`
  - `viridis`
  - `patchwork`

### **Data Requirements:**
- Input data should include:
  - **Sample IDs** (`sample_id`)
  - **Treatment Groups** (`treatment_group`)
  - **Cell Types** (`assigned_MP`)
  - **Proportions** of each cell type for each sample.

---

## **Usage**

### **1. Statistical Analysis**
The script `1_statistical_analysis.R` performs the following:
- Pairwise comparisons of cell type proportions between treatment groups using:
  - Wilcoxon rank-sum test
  - Student’s t-test
- Calculation of:
  - **Cohen's d** (effect size)
  - **Fold changes**
  - Adjusted p-values (Benjamini-Hochberg correction)

To use this script:
```r
source("1_statistical_analysis.R")

# Example:
results <- run_statistical_analysis(data = your_data)
```

### **2. Visualization**
The script `2_visualization.R` generates plots such as:
- Boxplots with significance annotations
- Heatmaps and effect size visualizations (extendable)

To use this script:
```r
source("2_visualization.R")

# Example:
plots <- create_statistical_plots(results, your_data)
save_plots(plots, prefix = "analysis_results")
```

### **3. Complete Workflow**
The script `3_run_analysis.R` integrates data preparation, analysis, and visualization:
- Extracts proportion data from a Seurat object.
- Runs statistical analysis using `run_statistical_analysis`.
- Generates and saves plots using `create_statistical_plots`.

To run the workflow:
```r
source("3_run_analysis.R")
```

---

## **Outputs**

1. **Statistical Results**:
   - A CSV file containing:
     - Means, standard deviations, fold changes, and Cohen's d for each cell type and treatment group comparison.
     - p-values and adjusted p-values for significance.

2. **Plots**:
   - Individual PDFs for:
     - Boxplots with significance annotations.
     - Heatmaps (optional).
     - Effect size plots.
   - A combined summary plot.

---

## **File Structure**

```plaintext
.
├── 1_statistical_analysis.R    # Functions for statistical analysis
├── 2_visualization.R           # Functions for visualization
├── 3_run_analysis.R            # Example usage and full workflow
├── statistical_results.csv     # Example output (if run)
└── README.md                   # This file
```

---

## **Key Features**

- **Comprehensive Statistical Analysis**: Includes p-values, fold changes, and Cohen's d to quantify and interpret differences between groups.
- **Customizable Visualizations**: Easily modify visual styles and annotations to fit your needs.
- **Workflow Flexibility**: Designed for integration with Seurat and scalable to large datasets.

---

## **Example Data**
- This repository assumes you have metadata from a Seurat object with the following structure:
  ```r
  your_data <- data.frame(
    sample_id = c("Sample1", "Sample2"),
    treatment_group = c("Treated", "Untreated"),
    assigned_MP = c("CellType1", "CellType2"),
    proportion = c(0.5, 0.6)
  )
  ```

---

## **Contributing**
Contributions are welcome! Please feel free to submit a pull request or open an issue to report bugs or suggest enhancements.

---

## **License**
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---


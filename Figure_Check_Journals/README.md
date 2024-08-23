# Figure Check Script for Scientific Journal Submissions

## Overview

This Python script automates the process of checking scientific figures for compliance with journal submission guidelines, specifically tailored for Nature Genetics. It helps researchers ensure their figures meet the required standards before submission, potentially saving time and reducing the likelihood of revision requests due to technical issues with figures.

## Why This Script?

Scientific journals often have strict requirements for figure submissions, including specifications for:

- Color mode (typically RGB)
- Resolution (usually minimum 300 DPI)
- Dimensions (e.g., maximum width for single or double column figures)
- Font type and size

Manually checking each figure for all these criteria can be time-consuming and error-prone. This script automates the process, providing a quick and reliable way to verify figure compliance.

## How This Script Helps

1. **Automated Checking**: Examines all figures in a specified directory for compliance with journal guidelines.
2. **Comprehensive Analysis**: Checks multiple criteria including RGB mode, resolution, dimensions, font type, and font size.
3. **Support for Multiple Formats**: Works with both image files (PNG, JPG, TIFF, BMP) and PDFs.
4. **Detailed Reporting**: Generates a Word document report summarizing the results for easy review and sharing.
5. **Time-Saving**: Reduces the time needed to manually inspect each figure.
6. **Error Reduction**: Minimizes the risk of overlooking issues that could lead to revision requests.

## Features

- Checks image files for RGB mode, resolution, and dimensions.
- Analyzes PDFs for embedded image resolution, dimensions, font type, and font size.
- Generates a detailed report in DOCX format, highlighting any issues found.
- Provides clear pass/fail indicators for each criterion checked.

## Requirements

- Python 3.6+
- Required Python packages: 
  - Pillow
  - PyMuPDF
  - python-docx

## Installation

1. Clone this repository:
   ```
   git clone https://github.com/yourusername/figure-check-script.git
   ```
2. Install required packages:
   ```
   pip install Pillow PyMuPDF python-docx
   ```

## Usage

1. Place all your figure files (PDFs and images) in a single directory.
2. Run the script from the command line:
   ```
   python figure_check_report.py /path/to/your/figures
   ```
3. The script will generate a report named `Figure_Check_Report.docx` in the same directory as your figures.

## Output

The generated report includes:
- A summary of the checks performed
- A table with results for each figure
- Explanatory notes about the criteria and how to interpret the results

## Limitations

- The script may not catch all possible issues, and manual review is still recommended.
- RGB mode checking for PDFs may require manual verification in some cases.
- The script is tailored for Nature Genetics guidelines; adjustments may be needed for other journals.

## Contributing

Contributions to improve the script or extend its functionality are welcome. Please feel free to submit pull requests or open issues for bugs and feature requests.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Developed to assist researchers in streamlining their journal submission process.
- Inspired by the submission guidelines of Nature Genetics and similar scientific journals.

# Wide XLSX to Multi-Page PDF Converter

This Python script converts wide Excel (XLSX) files into multi-page PDF documents, optimized for readability and efficient use of page space.

## Features

- Converts XLSX files to PDF format
- Handles wide spreadsheets by splitting columns across multiple pages
- Uses A3 landscape orientation for maximum column visibility
- Repeats header rows on each page for easy reference
- Adjusts font sizes to fit more data
- Provides detailed console output for tracking conversion progress

## Requirements

- Python 3.6 or higher
- openpyxl library
- reportlab library

## Installation

1. Ensure you have Python 3.6 or higher installed on your system.
2. Install the required libraries using pip:

```
pip install openpyxl reportlab
```

## Usage

1. Save the script as `wide_xlsx_to_pdf_converter.py` in your desired directory.
2. Open the script in a text editor and modify the following lines at the bottom of the file:

```python
xlsx_file = "path/to/your/input/file.xlsx"
pdf_file = "path/to/your/output/file.pdf"
```

Replace these paths with the actual path to your input XLSX file and the desired path for your output PDF file.

3. Open a terminal or command prompt.
4. Navigate to the directory containing the script.
5. Run the script using:

```
python wide_xlsx_to_pdf_converter.py
```

## Customization

You can customize the script by modifying the following variables:

- `max_columns_per_page`: Change this value to adjust the number of columns per page (default is 30).
- Table styles: Modify the `TableStyle` settings to change colors, fonts, or spacing.
- Page size: Change `landscape(A3)` to a different page size if needed.

## Troubleshooting

If you encounter any issues:

1. Ensure all dependencies are correctly installed.
2. Check that the input file path is correct and the file exists.
3. Make sure you have write permissions for the output PDF location.
4. If the script fails, check the console output for error messages and traceback information.

## Contributing

Feel free to fork this repository and submit pull requests with any enhancements.

## License

This script is released under the MIT License. See the LICENSE file for more details.

import openpyxl
from reportlab.lib import colors
from reportlab.lib.pagesizes import A3, landscape
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, PageBreak
from reportlab.lib.units import inch
import os
import sys
import traceback

def xlsx_to_pdf(xlsx_file, pdf_file):
    try:
        print(f"Starting conversion of {xlsx_file}")

        # Check if input file exists
        if not os.path.exists(xlsx_file):
            raise FileNotFoundError(f"Input file not found: {xlsx_file}")

        # Load the workbook and select the active worksheet
        print("Loading workbook...")
        wb = openpyxl.load_workbook(xlsx_file)
        ws = wb.active
        print(f"Workbook loaded. Active sheet: {ws.title}")

        # Get the data from the worksheet
        print("Extracting data from worksheet...")
        data = list(ws.values)
        print(f"Data extracted. Total rows: {len(data)}, Total columns: {len(data[0])}")

        # Create the PDF
        print(f"Creating PDF: {pdf_file}")
        doc = SimpleDocTemplate(pdf_file, pagesize=landscape(A3))
        elements = []

        # Function to create table
        def create_table(data_chunk, font_size):
            table = Table(data_chunk, repeatRows=1)
            style = TableStyle([
                ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
                ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
                ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
                ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                ('FONTSIZE', (0, 0), (-1, 0), font_size),
                ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
                ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
                ('TEXTCOLOR', (0, 1), (-1, -1), colors.black),
                ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
                ('FONTNAME', (0, 1), (-1, -1), 'Helvetica'),
                ('FONTSIZE', (0, 1), (-1, -1), font_size - 2),
                ('TOPPADDING', (0, 1), (-1, -1), 3),
                ('BOTTOMPADDING', (0, 1), (-1, -1), 3),
                ('GRID', (0, 0), (-1, -1), 0.5, colors.black)
            ])
            table.setStyle(style)
            return table

        # First page with at least 15 columns
        first_page_columns = min(20, len(data[0]))  # Try to fit up to 20 columns on first page
        first_page_data = [row[:first_page_columns] for row in data]
        table = create_table(first_page_data, 8)  # Larger font for first page
        elements.append(table)

        # Remaining columns on subsequent pages
        if len(data[0]) > first_page_columns:
            elements.append(PageBreak())
            remaining_columns = len(data[0]) - first_page_columns
            columns_per_page = 30  # Adjust as needed

            for i in range(first_page_columns, len(data[0]), columns_per_page):
                chunk = [row[i:i+columns_per_page] for row in data]
                table = create_table(chunk, 6)  # Smaller font for subsequent pages
                elements.append(table)
                if i + columns_per_page < len(data[0]):
                    elements.append(PageBreak())

        # Build the PDF
        print("Building PDF...")
        doc.build(elements)
        print(f"PDF created successfully: {pdf_file}")

    except Exception as e:
        print(f"An error occurred: {str(e)}")
        print("Traceback:")
        traceback.print_exc()

# Usage
xlsx_file = "/Users/akhaliq/Desktop/new_manuscript/Editors_final/Editors_check_V2/Suppli_tables/suppli_tables_v1/Supplementary Data 5.xlsx"
pdf_file = "/Users/akhaliq/Desktop/new_manuscript/Editors_final/Editors_check_V2/Suppli_tables/suppli_tables_v1/Supplementary Data 5.pdf"

print("Starting conversion process...")
xlsx_to_pdf(xlsx_file, pdf_file)
print("Process completed.")

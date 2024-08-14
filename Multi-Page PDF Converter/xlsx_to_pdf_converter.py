import openpyxl
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter, landscape
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle
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
        data = []
        for row in ws.iter_rows(values_only=True):
            data.append(row)
        print(f"Data extracted. Total rows: {len(data)}")

        # Create the PDF
        print(f"Creating PDF: {pdf_file}")
        doc = SimpleDocTemplate(pdf_file, pagesize=landscape(letter))
        elements = []

        # Create the table and add style
        print("Creating table and applying style...")
        table = Table(data)
        style = TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 14),
            ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
            ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
            ('TEXTCOLOR', (0, 1), (-1, -1), colors.black),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('FONTNAME', (0, 1), (-1, -1), 'Helvetica'),
            ('FONTSIZE', (0, 1), (-1, -1), 12),
            ('TOPPADDING', (0, 1), (-1, -1), 6),
            ('BOTTOMPADDING', (0, 1), (-1, -1), 6),
            ('GRID', (0, 0), (-1, -1), 1, colors.black)
        ])
        table.setStyle(style)

        # Add the table to the elements list
        elements.append(table)

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
pdf_file = "/Users/akhaliq/Desktop/Supplementary_Data_5.pdf"

print("Starting conversion process...")
xlsx_to_pdf(xlsx_file, pdf_file)
print("Process completed.")

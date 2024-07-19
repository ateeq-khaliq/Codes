import os
import base64
from jinja2 import Template
from PIL import Image
from pdf2image import convert_from_path
from plot_config import plot_files, output_dir, output_html, title, lab_name, proprietary_text
from plot_descriptions import get_plot_details

# Function to convert PDF to PNG
def convert_pdf_to_png(pdf_path, output_dir):
    base_name = os.path.splitext(os.path.basename(pdf_path))[0]
    png_path = os.path.join(output_dir, f"{base_name}.png")
    
    # Convert PDF to PNG
    pages = convert_from_path(pdf_path, 300)  # 300 DPI
    pages[0].save(png_path, 'PNG')
    
    return png_path

# Function to encode image to base64
def encode_image(image_path):
    with open(image_path, "rb") as image_file:
        encoded_string = base64.b64encode(image_file.read()).decode('utf-8')
    return f"data:image/png;base64,{encoded_string}"

# Create output directory
os.makedirs(output_dir, exist_ok=True)

# Convert PDFs to PNGs and encode all images
encoded_plots = []
for plot_file in plot_files:
    if plot_file.lower().endswith('.pdf'):
        png_path = convert_pdf_to_png(plot_file, output_dir)
    else:
        png_path = plot_file
    
    encoded_image = encode_image(png_path)
    filename = os.path.basename(plot_file)
    plot_title = os.path.splitext(filename)[0]
    
    plot_details = get_plot_details(plot_title)
    
    encoded_plots.append({
        'path': encoded_image,
        'title': plot_title,
        'heading': plot_details['heading'],
        'description': plot_details['description']
    })

# HTML template
html_template = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{{ title }}</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet">
    <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/js/bootstrap.bundle.min.js"></script>
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0-beta3/css/all.min.css">
    <style>
        @import url('https://fonts.googleapis.com/css2?family=Roboto:wght@300;400;700&display=swap');
        
        body, html {
            height: 100%;
            margin: 0;
            font-family: 'Roboto', sans-serif;
            background: linear-gradient(135deg, #1a2a6c, #b21f1f, #fdbb2d);
            background-size: 400% 400%;
            animation: gradientBG 15s ease infinite;
        }
        
        @keyframes gradientBG {
            0% { background-position: 0% 50%; }
            50% { background-position: 100% 50%; }
            100% { background-position: 0% 50%; }
        }
        
        .container-fluid {
            height: 100%;
            padding: 20px;
        }
        
        .header {
            text-align: center;
            margin-bottom: 30px;
            color: #ffffff;
        }
        
        h1 {
            padding: 20px 0 10px;
            font-weight: bold;
            text-shadow: 2px 2px 4px rgba(0,0,0,0.5);
            font-size: 2.5em;
            margin-bottom: 10px;
        }
        
        .lab-name {
            font-size: 1.8em;
            font-weight: 500;
            text-shadow: 1px 1px 2px rgba(0,0,0,0.5);
            margin-bottom: 10px;
        }
        
        .proprietary {
            font-style: italic;
            font-size: 0.9em;
            opacity: 0.8;
        }
        
        .content-wrapper {
            display: flex;
            height: calc(100vh - 200px);
            background-color: rgba(255,255,255,0.1);
            border-radius: 15px;
            overflow: hidden;
        }
        
        .nav-tabs {
            flex-basis: 200px;
            flex-shrink: 0;
            flex-direction: column;
            border: none;
            background-color: rgba(255,255,255,0.2);
            padding: 20px 0;
            overflow-y: auto;
        }
        
        .nav-tabs .nav-link {
            color: #ffffff;
            border: none;
            padding: 10px 15px;
            margin: 5px 10px;
            border-radius: 10px;
            transition: all 0.3s ease;
            text-align: left;
            white-space: nowrap;
            overflow: hidden;
            text-overflow: ellipsis;
            font-size: 0.9em;
        }
        
        .nav-tabs .nav-link:hover {
            background-color: rgba(255,255,255,0.3);
            transform: translateX(5px);
        }
        
        .nav-tabs .nav-link.active {
            background-color: #ffffff;
            color: #1a2a6c;
            font-weight: bold;
            box-shadow: 0 0 15px rgba(255,255,255,0.5);
        }
        
        .tab-content {
            flex-grow: 1;
            background-color: rgba(255,255,255,0.9);
            border-radius: 0 15px 15px 0;
            padding: 20px;
            overflow-y: auto;
        }
        
        .tab-pane {
            height: 100%;
            text-align: center;
        }
        
        .plot-container {
            position: relative;
            height: calc(100% - 150px);
            display: flex;
            justify-content: center;
            align-items: center;
        }
        
        .plot-image {
            max-width: 90%;
            max-height: 90%;
            object-fit: contain;
            box-shadow: 0 4px 8px rgba(0,0,0,0.1);
            border-radius: 10px;
            transition: all 0.3s ease;
        }
        
        .plot-image:hover {
            transform: scale(1.02);
            box-shadow: 0 8px 16px rgba(0,0,0,0.2);
        }
        
        .plot-heading {
            font-size: 1.8em;
            color: #1a2a6c;
            margin-bottom: 10px;
            font-weight: bold;
        }
        
        .plot-description {
            font-size: 1em;
            color: #333;
            margin-bottom: 20px;
            max-width: 800px;
            margin-left: auto;
            margin-right: auto;
        }
        
        .fullscreen-overlay {
            display: none;
            position: fixed;
            top: 0;
            left: 0;
            width: 100%;
            height: 100%;
            background-color: rgba(0,0,0,0.9);
            z-index: 1000;
            opacity: 0;
            transition: opacity 0.3s ease;
        }
        
        .fullscreen-content {
            width: 100%;
            height: 100%;
            display: flex;
            justify-content: center;
            align-items: center;
        }
        
        .fullscreen-content img {
            max-width: 95%;
            max-height: 95%;
            object-fit: contain;
            box-shadow: 0 0 30px rgba(255,255,255,0.1);
            border-radius: 10px;
        }
        
        .close-button {
            position: absolute;
            top: 20px;
            right: 30px;
            color: #f1f1f1;
            font-size: 40px;
            font-weight: bold;
            transition: 0.3s;
        }
        
        .close-button:hover,
        .close-button:focus {
            color: #bbb;
            text-decoration: none;
            cursor: pointer;
        }
        
        .controls {
            position: absolute;
            bottom: 20px;
            left: 50%;
            transform: translateX(-50%);
            display: flex;
            gap: 20px;
        }
        
        .control-button {
            background-color: rgba(255,255,255,0.2);
            border: none;
            color: white;
            padding: 10px 20px;
            border-radius: 25px;
            cursor: pointer;
            transition: all 0.3s ease;
        }
        
        .control-button:hover {
            background-color: rgba(255,255,255,0.3);
            transform: translateY(-2px);
        }
        
        @media (max-width: 768px) {
            .content-wrapper {
                flex-direction: column;
            }
            
            .nav-tabs {
                flex-basis: auto;
                flex-direction: row;
                overflow-x: auto;
                overflow-y: hidden;
            }
            
            .nav-tabs .nav-link {
                white-space: nowrap;
            }
            
            .tab-content {
                border-radius: 0 0 15px 15px;
            }
        }
    </style>
</head>
<body>
    <div class="container-fluid d-flex flex-column">
        <div class="header">
            <h1><i class="fas fa-dna"></i> {{ title }}</h1>
            <div class="lab-name">{{ lab_name }}</div>
            <div class="proprietary">{{ proprietary_text }}</div>
        </div>
        
        <div class="content-wrapper">
            <ul class="nav nav-tabs" id="plotTabs" role="tablist">
            {% for plot in plots %}
                <li class="nav-item" role="presentation">
                    <button class="nav-link {% if loop.first %}active{% endif %}" 
                            id="plot-{{ loop.index }}-tab" 
                            data-bs-toggle="tab" 
                            data-bs-target="#plot-{{ loop.index }}" 
                            type="button" 
                            role="tab" 
                            aria-controls="plot-{{ loop.index }}" 
                            aria-selected="{% if loop.first %}true{% else %}false{% endif %}"
                            title="{{ plot.title }}">
                        <i class="fas fa-chart-bar"></i> Plot {{ loop.index }}
                    </button>
                </li>
            {% endfor %}
            </ul>
            
            <div class="tab-content" id="plotTabsContent">
            {% for plot in plots %}
                <div class="tab-pane fade {% if loop.first %}show active{% endif %} h-100" 
                     id="plot-{{ loop.index }}" 
                     role="tabpanel" 
                     aria-labelledby="plot-{{ loop.index }}-tab">
                    <h2 class="plot-heading">{{ plot.heading }}</h2>
                    <p class="plot-description">{{ plot.description }}</p>
                    <div class="plot-container">
                        <img src="{{ plot.path }}" alt="{{ plot.title }}" class="plot-image" onclick="openFullscreen(this.src, {{ loop.index0 }})">
                    </div>
                </div>
            {% endfor %}
            </div>
        </div>
    </div>

    <div id="fullscreenOverlay" class="fullscreen-overlay">
        <span class="close-button" onclick="closeFullscreen()">&times;</span>
        <div class="fullscreen-content">
            <img id="fullscreenImage" src="" alt="Fullscreen Image">
        </div>
        <div class="controls">
            <button class="control-button" onclick="prevImage()"><i class="fas fa-chevron-left"></i> Previous</button>
            <button class="control-button" onclick="nextImage()">Next <i class="fas fa-chevron-right"></i></button>
        </div>
    </div>

    <script>
        let currentImageIndex = 0;
        const images = [
            {% for plot in plots %}
                "{{ plot.path }}"{% if not loop.last %},{% endif %}
            {% endfor %}
        ];

        function openFullscreen(imgSrc, index) {
            document.getElementById('fullscreenImage').src = imgSrc;
            document.getElementById('fullscreenOverlay').style.display = 'block';
            setTimeout(() => {
                document.getElementById('fullscreenOverlay').style.opacity = '1';
            }, 50);
            currentImageIndex = index;
        }

        function closeFullscreen() {
            document.getElementById('fullscreenOverlay').style.opacity = '0';
            setTimeout(() => {
                document.getElementById('fullscreenOverlay').style.display = 'none';
            }, 300);
        }

        function prevImage() {
            currentImageIndex = (currentImageIndex - 1 + images.length) % images.length;
            document.getElementById('fullscreenImage').src = images[currentImageIndex];
        }

        function nextImage() {
            currentImageIndex = (currentImageIndex + 1) % images.length;
            document.getElementById('fullscreenImage').src = images[currentImageIndex];
        }

        // Enable keyboard navigation
        document.addEventListener('keydown', function(e) {
            if (document.getElementById('fullscreenOverlay').style.display === 'block') {
                if (e.key === 'ArrowLeft') prevImage();
                if (e.key === 'ArrowRight') nextImage();
                if (e.key === 'Escape') closeFullscreen();
            }
        });
    </script>
</body>
</html>
"""

# Render the template
template = Template(html_template)
html_content = template.render(
    plots=encoded_plots,
    title=title,
    lab_name=lab_name,
    proprietary_text=proprietary_text
)

# Write the HTML file
with open(output_html, 'w') as f:
    f.write(html_content)

print(f"Self-contained HTML file has been generated at: {output_html}")
import os

# Snakemake input and output
input_files = snakemake.input
output_file = snakemake.output.html

# 分类对应结构（按路径关键词归类）
section_map = {
    "before_fpkm_filter.png": ("Expression Filtering", "Before FPKM Filter"),
    "after_fpkm_filter.png": ("Expression Filtering", "After FPKM Filter"),
    "gene_boxplot_FPKM.png": ("Expression Filtering", "Gene FPKM Boxplot"),

    "exon_number.png": ("Transcript Features", "Exon Number Distribution"),
    "gene-length-distrubution.png": ("Transcript Features", "Gene Length Distribution"),

    "gene_conservation_ecdf.png": ("Conservation Analysis", "Conservation ECDF"),

    "Fickett_score.png": ("Sequence Features", "Fickett Score"),
    "Max_ORF_length.png": ("Sequence Features", "Max ORF Length"),
    "Max_ORF_coverage.png": ("Sequence Features", "Max ORF Coverage"),
}

# 分组整理
image_sections = {}
for file in input_files:
    filename = os.path.basename(file)
    section, title = section_map.get(filename, ("Other", filename))
    image_sections.setdefault(section, []).append((title, file))

# 构建 HTML
html = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>NovelScan Analysis Report</title>
    <style>
        body {
            font-family: 'Helvetica', sans-serif;
            margin: 40px;
            background-color: #f9f9f9;
            color: #333;
        }
        h1 {
            color: #2c3e50;
            border-bottom: 2px solid #ccc;
            padding-bottom: 10px;
        }
        h2 {
            color: #34495e;
            margin-top: 40px;
        }
        .section {
            margin-bottom: 50px;
        }
        .image-container {
            margin-top: 20px;
            border: 1px solid #ccc;
            background: #fff;
            padding: 10px;
            box-shadow: 0 2px 5px rgba(0,0,0,0.1);
        }
        img {
            width: 100%;
            height: auto;
        }
        .toc {
            background-color: #ecf0f1;
            padding: 15px;
            border: 1px solid #bdc3c7;
            margin-bottom: 30px;
        }
        .toc a {
            text-decoration: none;
            color: #2980b9;
            display: block;
            margin: 5px 0;
        }
        .toc a:hover {
            text-decoration: underline;
        }
    </style>
</head>
<body>
    <h1>NovelScan Analysis Report</h1>
    
    <div class="toc">
        <h3>Table of Contents</h3>
"""

# 添加目录
for section in image_sections:
    html += f'<a href="#{section.replace(" ", "_")}">{section}</a>\n'

html += "</div>\n"

# 添加每节内容
for section, images in image_sections.items():
    html += f'<div class="section" id="{section.replace(" ", "_")}">\n'
    html += f'  <h2>{section}</h2>\n'
    for title, path in images:
        html += f'  <div class="image-container">\n'
        html += f'    <h3>{title}</h3>\n'
        html += f'    <img src="{path}" alt="{title}">\n'
        html += f'  </div>\n'
    html += '</div>\n'

html += """
</body>
</html>
"""

# 写入 HTML 输出文件
with open(output_file, "w") as f:
    f.write(html)

Genome = config["Genome"]
metadata = config["metadata"]
minGene = config["minGene"]
maxGene = config["maxGene"]
pctMT = config["pctMT"]
PCAnum = config["PCAnum"]
fpkm_filter = config["fpkm_filter"]
conservation_filter = config["conservation_filter"]
G0file = config["G0file"]
specie = config["specie"]
SingleR_ref = config["SingleR"]
group_chose = config["group_chose"]
resolution = config["resolution"]
celltypist_model = config["Celltypist"]
OutputDir=config["outdir"]

rule all:
	input:
		expand("{OutputDir}/snakemake/result/11_data_analysis/scRNA.rds",OutputDir=OutputDir),
		expand("{OutputDir}/snakemake/result/11_data_analysis/cluster_diff/unannotated_marker.csv",OutputDir=OutputDir),
		expand("{OutputDir}/snakemake/result/11_data_analysis/conservation/ConservedScore.txt",OutputDir=OutputDir),
		expand("{OutputDir}/snakemake/result/11_data_analysis/cluster/tSNE_{group_chose}.pdf",group_chose=group_chose,OutputDir=OutputDir),
		expand("{OutputDir}/snakemake/result/11_data_analysis/celltype/scRNA_celltype.rds",OutputDir=OutputDir),
		expand("{OutputDir}/snakemake/result/11_data_analysis/Tissue_specificity/celltype_specificity.pdf",OutputDir=OutputDir),
		expand("{OutputDir}/snakemake/result/11_data_analysis/Tissue_specificity/tissue_heatmap.pdf",OutputDir=OutputDir),
		expand("{OutputDir}/snakemake/result/11_data_analysis/Tissue_specificity/celltype_heatmap.pdf",OutputDir=OutputDir),
		expand("{OutputDir}/snakemake/result/11_data_analysis/cluster/lncRNA/tSNE_group.pdf",OutputDir=OutputDir),
		expand("{OutputDir}/snakemake/result/11_data_analysis/cluster/protein/tSNE_group.pdf",OutputDir=OutputDir),
		expand("{OutputDir}/snakemake/result/11_data_analysis/cluster/tsne_lncRNA.pdf",OutputDir=OutputDir),
		expand("snakemake/result/CellInfo_report.html")
rule cell_cluster:
	input:
		"{OutputDir}/snakemake/result/11_data_analysis/finally_tpm.txt",
		"{OutputDir}/snakemake/result/10_featurecount/all.txt"
	output:
		"{OutputDir}/snakemake/result/11_data_analysis/QC/vlnplot_before_qc.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/QC/vlnplot_after_qc.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/VariableFeatures.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/pca.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/cell_cluster.csv",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/tSNE.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/UMAP.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/scRNA.rds",
		"{OutputDir}/snakemake/result/11_data_analysis/QC/pearplot.pdf"
	message:
		"QC and cluster"
	conda:
		"env/R_Library.yaml"
	shell:
		"Rscript snakemake/script/cell_cluster.R -n {minGene} -p {maxGene} -q {pctMT} -t {PCAnum} -r {resolution} -a {input[0]} -b {input[1]} -j {metadata} -c {output[0]} -k {output[8]} -d {output[1]} -e {output[2]} -f {output[3]} -g {output[4]} -l {output[5]} -i {output[6]} -m {output[7]}"
rule cell_marker_gene:
	input:
		"{OutputDir}/snakemake/result/11_data_analysis/scRNA.rds",
		"{OutputDir}/snakemake/result/10_featurecount/unannotated_lncRNA.csv",
		"{OutputDir}/snakemake/result/10_featurecount/unannotated_protein.csv" 
	output:
		"{OutputDir}/snakemake/result/11_data_analysis/cluster_diff/diff_genes_wilcox.csv",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster_diff/top10_diff_genes_wilcox.csv",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster_diff/unannotated_marker_vloplot.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster_diff/unannotated_marker_featureplot.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster_diff/unannotated_marker.csv",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster_diff/top10_marker.pdf"
	message:
		"cluster marker gene identify"
	conda:
		"env/R_Library.yaml"
	shell:
		"Rscript snakemake/script/cell_marker_gene.R -a {input[0]} -f {input[1]} -b {output[0]} -c {output[1]} -d {output[2]} -e {output[3]} -g {output[4]} -i {output[5]} -k {input[2]}"	
rule conservedscore:
	input:
		"{OutputDir}/snakemake/result/11_data_analysis/conservation/ConservedScore.txt",
		"{OutputDir}/snakemake/result/10_featurecount/all.txt",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster_diff/unannotated_marker.csv",
		"{OutputDir}/snakemake/result/11_data_analysis/scRNA.rds"
	output:
		"{OutputDir}/snakemake/result/11_data_analysis/conservation/gene_conservation_density.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/conservation/gene_conservation_ecdf.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/conservation/unannotated_score.txt",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster_diff/unannotated_marker_vloplot_conservation.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster_diff/unannotated_marker_featureplot_conservation.pdf"
	message:
		"calculate gene conservation"
	conda:
		"env/R_Library.yaml"
	shell:
		"""
		Rscript snakemake/script/conservation.R -a {input[0]} -b {input[1]} -c {output[0]} -d {output[1]} -e {wildcards.OutputDir}/snakemake/result/11_data_analysis/conservation/all_gene.gtf -f {wildcards.OutputDir}/snakemake/result/11_data_analysis/conservation/unannotated_score.txt && \
		awk -F, '{{print $1}}' {input[2]} | uniq > {wildcards.OutputDir}/snakemake/result/11_data_analysis/conservation/unannotated_id.txt && \
		cat {wildcards.OutputDir}/snakemake/result/11_data_analysis/conservation/unannotated_id.txt |sed 's/\"//g'|uniq > {wildcards.OutputDir}/snakemake/result/11_data_analysis/conservation/unannotated_marker_id.txt && \
		sed -i '1d' {wildcards.OutputDir}/snakemake/result/11_data_analysis/conservation/unannotated_marker_id.txt && \
		sed -i 's/-/_/g' {wildcards.OutputDir}/snakemake/result/11_data_analysis/conservation/unannotated_marker_id.txt && \
		grep -F -f {wildcards.OutputDir}/snakemake/result/11_data_analysis/conservation/unannotated_marker_id.txt {output[2]} > {wildcards.OutputDir}/snakemake/result/11_data_analysis/conservation/unannotated_marker_score.txt && \
		Rscript snakemake/script/unannotated_score.R -a {input[3]} -d {output[3]} -e {output[4]} -i {input[2]} -I {output[2]} -k {wildcards.OutputDir}/snakemake/result/11_data_analysis/cluster_diff/unannotated_marker_conservation.csv -f {fpkm_filter} -c {conservation_filter} && \
		awk -F, '{{print $NF}}' {wildcards.OutputDir}/snakemake/result/11_data_analysis/cluster_diff/unannotated_marker_conservation.csv | uniq > {wildcards.OutputDir}/snakemake/result/11_data_analysis/unannotated_marker_conservation.txt && \
		sed -i '1d' {wildcards.OutputDir}/snakemake/result/11_data_analysis/unannotated_marker_conservation.txt
		"""
rule repeatmasker:
	input:
		"{OutputDir}/snakemake/result/11_data_analysis/unannotated_lncRNA_finally.fa",
		"{OutputDir}/snakemake/result/11_data_analysis/unannotated_protein_finally.fa",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster_diff/unannotated_marker.csv",
		"{OutputDir}/snakemake/result/11_data_analysis/unannotated_lncRNA_finally.gtf",
		"{OutputDir}/snakemake/result/11_data_analysis/unannotated_protein_finally.gtf"
	output:
		"{OutputDir}/snakemake/result/11_data_analysis/repeakmasker/unannotated_marker_finally.fa"
	message:
		"unannotated marker gene repeatmasker"
	shell:
		"""
		cat {input[0]} {input[1]} > {wildcards.OutputDir}/snakemake/result/11_data_analysis/repeakmasker/unannotated_gene_finallly.fa && \
		cat {input[3]} {input[4]} > {wildcards.OutputDir}/snakemake/result/11_data_analysis/repeakmasker/unannotated_gene_finallly.gtf && \
		awk -F, '{{print $1}}' {input[2]} > {wildcards.OutputDir}/snakemake/result/11_data_analysis/repeakmasker/unannotated_gene.txt && \
		cat {wildcards.OutputDir}/snakemake/result/11_data_analysis/repeakmasker/unannotated_gene.txt |sed 's/\"//g'|uniq > {wildcards.OutputDir}/snakemake/result/11_data_analysis/repeakmasker/unannotated_gene_1.txt && \
		sed -i '1d' {wildcards.OutputDir}/snakemake/result/11_data_analysis/repeakmasker/unannotated_gene_1.txt && \
		sed -i 's/-/_/g' {wildcards.OutputDir}/snakemake/result/11_data_analysis/repeakmasker/unannotated_gene_1.txt && \
		if [[ ! -f "{wildcards.OutputDir}/snakemake/result/11_data_analysis/repeakmasker/unannotated_marker_finallly.gtf" ]];then echo "file no exist";else rm -f {wildcards.OutputDir}/snakemake/result/11_data_analysis/repeakmasker/unannotated_marker_finallly.gtf;fi && \
		for trs in `cat {wildcards.OutputDir}/snakemake/result/11_data_analysis/repeakmasker/unannotated_gene_1.txt`;do grep "$trs" {wildcards.OutputDir}/snakemake/result/11_data_analysis/repeakmasker/unannotated_gene_finallly.gtf >> {wildcards.OutputDir}/snakemake/result/11_data_analysis/repeakmasker/unannotated_marker_finallly.gtf;done && \
		awk '{{print $10}}' {wildcards.OutputDir}/snakemake/result/11_data_analysis/repeakmasker/unannotated_marker_finallly.gtf > {wildcards.OutputDir}/snakemake/result/11_data_analysis/repeakmasker/unannotated_transcript_id.txt && \
		cat {wildcards.OutputDir}/snakemake/result/11_data_analysis/repeakmasker/unannotated_transcript_id.txt |sed 's/\"//g'|uniq > {wildcards.OutputDir}/snakemake/result/11_data_analysis/repeakmasker/unannotated_transcript_id_1.txt && \
		sed 's/;$//' {wildcards.OutputDir}/snakemake/result/11_data_analysis/repeakmasker/unannotated_transcript_id_1.txt|uniq > {wildcards.OutputDir}/snakemake/result/11_data_analysis/repeakmasker/unannotated_transcript_id_2.txt && \
		python snakemake/script/extract_fasta.py {wildcards.OutputDir}/snakemake/result/11_data_analysis/repeakmasker/unannotated_gene_finallly.fa {wildcards.OutputDir}/snakemake/result/11_data_analysis/repeakmasker/unannotated_transcript_id_2.txt {output} && \
		if [[ ! -d "{wildcards.OutputDir}/snakemake/result/11_data_analysis/repeakmasker/output/" ]]; then mkdir {wildcards.OutputDir}/snakemake/result/11_data_analysis/repeakmasker/output/;else echo "file exist";fi && \
		RepeatMasker -species "{specie}" -poly -html -gff -dir {wildcards.OutputDir}/snakemake/result/11_data_analysis/repeakmasker/output {output} 
		"""
rule cell_metadata:
	input:
		"{OutputDir}/snakemake/result/11_data_analysis/scRNA.rds"
	output:
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/pca_{group_chose}.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/tSNE_{group_chose}.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/UMAP_{group_chose}.pdf"
	conda:
		"env/R_Library.yaml"
	shell:
		"Rscript snakemake/script/cell_metadata.R -f {output[0]} -l {output[1]} -i {output[2]} -m {input} -j {group_chose}"
if config["celltypeMethod"] == 'true':
	rule celltype:
		input:
			"{OutputDir}/snakemake/result/11_data_analysis/scRNA.rds",
			"{OutputDir}/snakemake/result/10_featurecount/all.txt"
		output:
			"{OutputDir}/snakemake/result/11_data_analysis/celltype/celltype.csv",
			"{OutputDir}/snakemake/result/11_data_analysis/celltype/tSNE_celltype.pdf",
			"{OutputDir}/snakemake/result/11_data_analysis/celltype/UMAP_celltype.pdf",
			"{OutputDir}/snakemake/result/11_data_analysis/celltype/scRNA_celltype.rds"
		message:
			"identity cell type"
		conda:
			"env/celltypist.yaml"
		shell:
			"Rscript snakemake/script/celltypist.R -a {input[0]} -b {celltypist_model} -d {output[0]} -e {output[1]} -f {output[2]} -g {input[1]} -m {output[3]} -t {OutputDir}"
else:
	rule celltype:
		input:
			"{OutputDir}/snakemake/result/11_data_analysis/scRNA.rds",
			"{OutputDir}/snakemake/result/10_featurecount/all.txt"
		output:
			"{OutputDir}/snakemake/result/11_data_analysis/celltype/celltype.csv",
			"{OutputDir}/snakemake/result/11_data_analysis/celltype/tSNE_celltype.pdf",
			"{OutputDir}/snakemake/result/11_data_analysis/celltype/UMAP_celltype.pdf",
			"{OutputDir}/snakemake/result/11_data_analysis/celltype/scRNA_celltype.rds"
		conda:
			"env/R_Library.yaml"
		message:
			"identity cell type"
		shell:
			"Rscript snakemake/script/celltype.R -a {input[0]} -b {SingleR_ref}  -d {output[0]} -e {output[1]} -f {output[2]} -g {input[1]} -m {output[3]} -t {OutputDir}"
rule tissue_celltype_specificity:
	input:
		"{OutputDir}/snakemake/result/11_data_analysis/finally_tpm.txt",
		"{OutputDir}/snakemake/result/11_data_analysis/celltype/scRNA_celltype.rds"    
	output:
		"{OutputDir}/snakemake/result/11_data_analysis/Tissue_specificity/celltype_specificity.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/Tissue_specificity/tissue_specificity.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/Tissue_specificity/celltype_specificity_ecdf.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/Tissue_specificity/tissue_specificity_ecdf.pdf"
	message:
		"calculate tissue/cell sepcificity score"
	conda:
		"env/TissueSpecificity.yaml"
	shell:
		"Rscript snakemake/script/tissue_specificity.R -a {input[0]}  -c {input[1]} -i {output[0]} -I {output[1]} -o {output[2]} -O {output[3]} -d {group_chose}" 
rule tissue_specificity_heatmap:
	input:
		"{OutputDir}/snakemake/result/11_data_analysis/finally_tpm.txt",
		"{OutputDir}/snakemake/result/11_data_analysis/celltype/scRNA_celltype.rds",
		"{OutputDir}/snakemake/result/11_data_analysis/Tissue_specificity/celltype_specificity.pdf"
	output:
		"{OutputDir}/snakemake/result/11_data_analysis/Tissue_specificity/tissue_heatmap.pdf"
	message:
		"plot tissue specificity heatmap"
	conda:
		"env/TissueSpecificity.yaml"
	shell:
		"Rscript snakemake/script/tissue_specificity_heatmap.R -a {input[0]} -c {input[1]} -o {output} -b {group_chose}"
rule celltype_specificity_heatmap:
	input:
		"{OutputDir}/snakemake/result/11_data_analysis/finally_tpm.txt",
		"{OutputDir}/snakemake/result/11_data_analysis/celltype/scRNA_celltype.rds",
		"{OutputDir}/snakemake/result/11_data_analysis/Tissue_specificity/celltype_specificity.pdf"
	output:
		"{OutputDir}/snakemake/result/11_data_analysis/Tissue_specificity/celltype_heatmap.pdf"
	message:
		"plot celltype specificity heatmap"
	conda:
		"env/TissueSpecificity.yaml"
	shell:
		"Rscript snakemake/script/celltype_specificity_heatmap.R -a {input[0]} -c {input[1]} -o {output}"
rule lncRNA_cluster:
	input:
		"{OutputDir}/snakemake/result/11_data_analysis/finally_tpm.txt",
		"{OutputDir}/snakemake/result/11_data_analysis/celltype/scRNA_celltype.rds",
		"{OutputDir}/snakemake/result/11_data_analysis/Tissue_specificity/celltype_heatmap.pdf"
	output:
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/lncRNA/VariableFeatures.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/lncRNA/pca.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/lncRNA/cell_cluster.csv",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/lncRNA/tSNE.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/lncRNA/UMAP.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/lncRNA/top10_marker.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/lncRNA/tSNE_group.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/lncRNA/UMAP_group.pdf"
	message:
		"cluster data by lncRNA"
	conda:
		"env/R_Library.yaml"
	shell:
		"Rscript snakemake/script/lncRNA_cluster.R -a {input[0]} -b {group_chose} -j {resolution} -c {input[1]} -k {PCAnum} -e {output[0]} -f {output[1]} -g {output[2]} -l {output[3]} -i {output[4]} -m {output[5]} -n {output[6]} -o {output[7]} -t {OutputDir}"
rule protein_coding_cluster:
	input:
		"{OutputDir}/snakemake/result/11_data_analysis/finally_tpm.txt",
		"{OutputDir}/snakemake/result/11_data_analysis/celltype/scRNA_celltype.rds",
		"{OutputDir}/snakemake/result/11_data_analysis/Tissue_specificity/celltype_heatmap.pdf"
	output:
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/protein/VariableFeatures.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/protein/pca.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/protein/cell_cluster.csv",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/protein/tSNE.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/protein/UMAP.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/protein/top10_marker.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/protein/tSNE_group.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/protein/UMAP_group.pdf"
	message:
		"cluster data by protein coding gene"
	conda:
		"env/R_Library.yaml"
	shell:
		"Rscript snakemake/script/protein_cluster.R -a {input[0]} -b {group_chose} -j {resolution} -c {input[1]} -k {PCAnum} -e {output[0]} -f {output[1]} -g {output[2]} -l {output[3]} -i {output[4]} -m {output[5]} -n {output[6]} -o {output[7]} -t {OutputDir}"
rule lnc_protein:
	input:
		"{OutputDir}/snakemake/result/11_data_analysis/celltype/scRNA_celltype.rds",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/lncRNA/cell_cluster.csv",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/protein/cell_cluster.csv"
	output:
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/tsne_lncRNA.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/umap_lncRNA.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/tsne_protein.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/umap_protein.pdf"
	conda:
		"env/R_Library.yaml"
	shell:
		"Rscript snakemake/script/lncRNA_protein_cluster.R -a {input[0]} -b {input[1]} -c {input[2]} -e {output[0]} -f {output[1]} -n {output[2]} -o {output[3]}"
rule merge_QC_pdf:
	input:
		"{OutputDir}/snakemake/result/11_data_analysis/QC/vlnplot_before_qc.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/QC/vlnplot_after_qc.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/VariableFeatures.pdf"
	output:
		"{OutputDir}/snakemake/result/11_data_analysis/QC.pdf"
	shell:
		"""
		mutool merge -o {output} {input[0]} {input[1]} {input[2]}
		"""
rule merge_cluster_pdf:
	input:
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/pca.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/tSNE.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/UMAP.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster_diff/top10_marker.pdf",
		expand("{OutputDir}/snakemake/result/11_data_analysis/cluster/pca_{group_chose}.pdf",group_chose=group_chose,OutputDir=OutputDir),
		expand("{OutputDir}/snakemake/result/11_data_analysis/cluster/tsne_{group_chose}.pdf",group_chose=group_chose,OutputDir=OutputDir),
		expand("{OutputDir}/snakemake/result/11_data_analysis/cluster/umap_{group_chose}.pdf",group_chose=group_chose,OutputDir=OutputDir)
	output:
		"{OutputDir}/snakemake/result/11_data_analysis/cluster.pdf"
	shell:
		"""
		mutool merge -o {output} {input[0]} {input[1]} {input[2]} {input[3]} {input[4]} {input[5]} {input[6]}
		"""
rule celltype_pdf:
	input:
		"{OutputDir}/snakemake/result/11_data_analysis/celltype/tSNE_celltype.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/celltype/UMAP_celltype.pdf"
	output:
		"{OutputDir}/snakemake/result/11_data_analysis/celltype.pdf"
	shell:
		"""
		mutool merge -o {output} {input[0]} {input[1]}
		"""	
rule tissue_specificity_pdf:
	input:
		"{OutputDir}/snakemake/result/11_data_analysis/Tissue_specificity/celltype_specificity.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/Tissue_specificity/tissue_specificity.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/Tissue_specificity/celltype_specificity_ecdf.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/Tissue_specificity/tissue_specificity_ecdf.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/Tissue_specificity/tissue_protein_coding_heatmap.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/Tissue_specificity/tissue_annotated_lncRNA_heatmap.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/Tissue_specificity/tissue_unannotated_lncRNA_heatmap.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/Tissue_specificity/tissue_unannotated_protein_heatmap.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/Tissue_specificity/celltype_protein_coding_heatmap.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/Tissue_specificity/celltype_annotated_lncRNA_heatmap.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/Tissue_specificity/celltype_unannotated_lncRNA_heatmap.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/Tissue_specificity/celltype_unannotated_protein_heatmap.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/conservation/gene_conservation_density.pdf",
        "{OutputDir}/snakemake/result/11_data_analysis/conservation/gene_conservation_ecdf.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster_diff/unannotated_marker_vloplot_conservation.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster_diff/unannotated_marker_featureplot_conservation.pdf"
	output:
		"{OutputDir}/snakemake/result/11_data_analysis/uannotated_gene_feature.pdf"	
	shell:
		"""
		mutool merge -o {output} {input[0]} {input[1]} {input[2]} {input[3]} {input[4]} {input[5]} {input[6]} {input[7]} {input[8]} {input[9]} {input[10]} {input[11]} {input[12]} {input[13]} {input[14]} {input[15]}
		"""	
rule lncRNA_pdf:        
	input:
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/lncRNA/VariableFeatures.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/lncRNA/pca.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/lncRNA/tSNE.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/lncRNA/UMAP.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/lncRNA/top10_marker.pdf",
		expand("{OutputDir}/snakemake/result/11_data_analysis/cluster/lncRNA/tSNE_{group_chose}.pdf",group_chose=group_chose,OutputDir=OutputDir),
		expand("{OutputDir}/snakemake/result/11_data_analysis/cluster/lncRNA/UMAP_{group_chose}.pdf",group_chose=group_chose,OutputDir=OutputDir)
	output:
		"{OutputDir}/snakemake/result/11_data_analysis/lncRNA_cluster.pdf"
	shell:
		"""
		mutool merge -o {output} {input[0]} {input[1]} {input[2]} {input[3]} {input[4]} {input[5]} {input[6]} 
		"""
rule protein_pdf:
	input:
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/protein/VariableFeatures.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/protein/pca.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/protein/tSNE.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/protein/UMAP.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/protein/top10_marker.pdf",
		expand("{OutputDir}/snakemake/result/11_data_analysis/cluster/protein/tSNE_{group_chose}.pdf",group_chose=group_chose,OutputDir=OutputDir),
		expand("{OutputDir}/snakemake/result/11_data_analysis/cluster/protein/UMAP_{group_chose}.pdf",group_chose=group_chose,OutputDir=OutputDir)
	output:
		"{OutputDir}/snakemake/result/11_data_analysis/protein_cluster.pdf"
	shell:
		"""
		mutool merge -o {output} {input[0]} {input[1]} {input[2]} {input[3]} {input[4]} {input[5]} {input[6]}
		"""	
rule make_report:
	input:
		"{OutputDir}/snakemake/result/11_data_analysis/celltype/scRNA_celltype.rds",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/tsne_lncRNA.pdf",
		"{OutputDir}/snakemake/result/10_featurecount/unannotated_lncRNA.csv",
		"{OutputDir}/snakemake/result/11_data_analysis/conservation/ConservedScore.txt",
		"{OutputDir}/snakemake/result/10_featurecount/unannotated_protein.csv",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/protein/cell_cluster.csv",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster/lncRNA/cell_cluster.csv",
		"{OutputDir}/snakemake/result/11_data_analysis/Tissue_specificity/celltype_specificity_ecdf.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/Tissue_specificity/tissue_heatmap.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/Tissue_specificity/celltype_heatmap.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/conservation/gene_conservation_ecdf.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/cluster_diff/top10_marker.pdf"
	output:
		"{OutputDir}/snakemake/result/CellInfo_report.html"
	conda:
		"env/R_Library.yaml"
	params:
		finally_tpm="{OutputDir}/snakemake/result/11_data_analysis/finally_tpm.txt",
		all_id="{OutputDir}/snakemake/result/10_featurecount/all.txt",
		scrna="{OutputDir}/snakemake/result/11_data_analysis/celltype/scRNA_celltype.rds",
		unannotatedLncRNA="{OutputDir}/snakemake/result/10_featurecount/unannotated_lncRNA.csv",
		unannotatedProtein="{OutputDir}/snakemake/result/10_featurecount/unannotated_protein.csv",
		ConservedScore="{OutputDir}/snakemake/result/11_data_analysis/conservation/ConservedScore.txt",
		gtf_all="{OutputDir}/snakemake/result/10_featurecount/all.gtf",
		proteinCluster="{OutputDir}/snakemake/result/11_data_analysis/cluster/protein/cell_cluster.csv",
		lncRNACluster="{OutputDir}/snakemake/result/11_data_analysis/cluster/lncRNA/cell_cluster.csv",
		group_chose={group_chose},
		lncRNA_rds="{OutputDir}/snakemake/result/11_data_analysis/cluster/lncRNA/lncRNA.rds",
		protein_rds="{OutputDir}/snakemake/result/11_data_analysis/cluster/protein/protein.rds",
		AllCluster="{OutputDir}/snakemake/result/11_data_analysis/cluster/cell_cluster.csv",
		Outputdir={OutputDir}
	script:
		"script/CellInfo_group.Rmd"
rule makr_report2:
	input:
		expand("{OutputDir}/snakemake/result/CellInfo_report.html",OutputDir=OutputDir)
	output:
		"snakemake/result/CellInfo_report.html"
	shell:
		"cp -r {OutputDir}/snakemake/result/CellInfo_report.html {output}"
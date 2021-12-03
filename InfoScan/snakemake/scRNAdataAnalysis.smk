configfile: "snakemake/config/DataAnalysis_config.yml"
(SAMPLES,READS,) = glob_wildcards("snakemake/result/result/2.fastp_output/{sample}_clean_{read}.fq.gz")
READS=["1","2"]
Genome = config["Genome"]
metadata = config["metadata"]
minGene = config["minGene"]
maxGene = config["maxGene"]
pctMT = config["pctMT"]
PCAnum = config["PCAnum"]
pctERCC = config["pctERCC"]
fpkm_filter = config["fpkm_filter"]
conservation_filter = config["conservation_filter"]
phastcons = config["phastcons"]
G0file = config["G0file"]
specie = config["specie"]
SingleR_ref1 = config["SingleR_ref1"]
SingleR_ref2 = config["SingleR_ref2"]
group_chose = config["group_chose"]
resolution = config["resolution"]
rule all:
	input:
		expand("snakemake/result/11_data_analysis/scRNA.rds"),
		expand("snakemake/result/11_data_analysis/cluster_diff/unannotated_marker.csv"),
		expand("snakemake/result/11_data_analysis/conservation/ConservedScore.txt"),
		expand("snakemake/result/11_data_analysis/expression.txt"),
		expand("snakemake/result/11_data_analysis/repeakmasker/unannotated_marker_finally.fa"),
		expand("snakemake/result/11_data_analysis/cluster/tsne_{group_chose}.pdf",group_chose=group_chose),
		expand("snakemake/result/11_data_analysis/celltype/scRNA_celltype.rds"),
		expand("snakemake/result/11_data_analysis/Tissue_specificity/celltype_specificity.pdf"),
		expand("snakemake/result/11_data_analysis/Tissue_specificity/tissue_unannotated_lncRNA_heatmap.pdf"),
		expand("snakemake/result/11_data_analysis/Tissue_specificity/celltype_unannotated_lncRNA_heatmap.pdf"),
		expand("snakemake/result/11_data_analysis/cluster/lncRNA/tSNE_{group_chose}.pdf",group_chose=group_chose),
		expand("snakemake/result/11_data_analysis/cluster/protein/tSNE_{group_chose}.pdf",group_chose=group_chose),
		expand("snakemake/result/11_data_analysis/QC.pdf"),
		expand("snakemake/result/11_data_analysis/cluster.pdf"),
		expand("snakemake/result/11_data_analysis/celltype.pdf"),
		expand("snakemake/result/11_data_analysis/uannotated_gene_feature.pdf"),
		expand("snakemake/result/11_data_analysis/lncRNA_cluster.pdf"),
		expand("snakemake/result/11_data_analysis/protein_cluster.pdf")

rule cell_cluster:
	input:
		"snakemake/result/11_data_analysis/finally_fpkm.txt",
		"snakemake/result/10_featurecount/all.txt"
	output:
		"snakemake/result/11_data_analysis/QC/vlnplot_before_qc.pdf",
		"snakemake/result/11_data_analysis/QC/vlnplot_after_qc.pdf",
		"snakemake/result/11_data_analysis/cluster/VariableFeatures.pdf",
		"snakemake/result/11_data_analysis/cluster/pca.pdf",
		"snakemake/result/11_data_analysis/cluster/cell_cluster.csv",
		"snakemake/result/11_data_analysis/cluster/tSNE.pdf",
		"snakemake/result/11_data_analysis/cluster/UMAP.pdf",
		"snakemake/result/11_data_analysis/scRNA.rds",
		"snakemake/11_data_analysis/QC/pearplot.pdf"
	message:
		"QC and cluster"
	shell:
		"Rscript snakemake/script/cell_cluster.R -n {minGene} -p {maxGene} -q {pctMT} -t {PCAnum} -s {pctERCC} -r {resolution} -a {input[0]} -b {input[1]} -j {metadata} -c {output[0]} -k {output[8]} -d {output[1]} -e {output[2]} -f {output[3]} -g {output[4]} -l {output[5]} -i {output[6]} -m {output[7]}"
rule cell_marker_gene:
	input:
		"snakemake/result/11_data_analysis/scRNA.rds",
		"snakemake/result/10_featurecount/unannotated_lncRNA.csv",
		"snakemake/result/10_featurecount/unannotated_protein.csv" 
	output:
		"snakemake/result/11_data_analysis/cluster_diff/diff_genes_wilcox.csv",
		"snakemake/result/11_data_analysis/cluster_diff/top10_diff_genes_wilcox.csv",
		"snakemake/result/11_data_analysis/cluster_diff/unannotated_marker_vloplot.pdf",
		"snakemake/result/11_data_analysis/cluster_diff/unannotated_marker_featureplot.pdf",
		"snakemake/result/11_data_analysis/cluster_diff/unannotated_marker.csv",
		"snakemake/result/11_data_analysis/cluster_diff/top10_marker.pdf"
	message:
		"cluster marker gene identify"
	shell:
		"Rscript snakemake/script/cell_marker_gene.R -a {input[0]} -f {input[1]} -b {output[0]} -c {output[1]} -d {output[2]} -e {output[3]} -g {output[4]} -i {output[5]} -k {input[2]}"	
rule conservedscore:
	input:
		"snakemake/result/11_data_analysis/unannotated_lncRNA_finally.gtf",
		"snakemake/result/11_data_analysis/unannotated_protein_finally.gtf",
		"snakemake/result/11_data_analysis/annotated_lncRNA_finally.gtf",
		"snakemake/result/11_data_analysis/protein_coding_finally.gtf",
		"snakemake/result/10_featurecount/all.txt",
		"snakemake/result/11_data_analysis/cluster_diff/unannotated_marker.csv"
	output:
		"snakemake/result/11_data_analysis/conservation/all_transcript.bed",
		"snakemake/result/11_data_analysis/conservation/ConservedScore.txt",
		"snakemake/result/11_data_analysis/conservation/gene_conservation_density.pdf",
		"snakemake/result/11_data_analysis/conservation/gene_conservation_ecdf.pdf",
		"snakemake/result/11_data_analysis/conservation/unannotated_score.txt"
	message:
		"calculate gene conservation"
	shell:
		"""
		cat {input[0]} {input[1]} {input[2]} {input[3]} > snakemake/result/11_data_analysis/conservation/all_gene.gtf && \
        	snakemake/script/gtfToGenePred snakemake/result/11_data_analysis/conservation/all_gene.gtf snakemake/result/11_data_analysis/conservation/all_gene.genePred && \
		snakemake/script/genePredToBed snakemake/result/11_data_analysis/conservation/all_gene.genePred {output[0]} && \
        	snakemake/script/bigWigAverageOverBed {phastcons} {output[0]}  {output[1]} && \
		Rscript snakemake/script/conservation.R -a {output[1]} -b {input[4]} -c {output[2]} -d {output[3]} -e snakemake/result/11_data_analysis/conservation/all_gene.gtf -f snakemake/result/11_data_analysis/conservation/unannotated_score.txt && \
		awk -F, '{{print $1}}' {input[5]} | uniq > snakemake/result/11_data_analysis/conservation/unannotated_id.txt && \
		cat snakemake/result/11_data_analysis/conservation/unannotated_id.txt |sed 's/\"//g'|uniq > snakemake/result/11_data_analysis/conservation/unannotated_marker_id.txt && \
		sed -i '1d' snakemake/result/11_data_analysis/conservation/unannotated_marker_id.txt && \
		sed -i 's/-/_/g' snakemake/result/11_data_analysis/conservation/unannotated_marker_id.txt && \
		grep -F -f snakemake/result/11_data_analysis/conservation/unannotated_marker_id.txt {output[4]} > snakemake/result/11_data_analysis/conservation/unannotated_marker_score.txt && \
		Rscript snakemake/script/unannotated_score.R -i {input[5]} -I {output[4]} -k snakemake/result/11_data_analysis/cluster_diff/unannotated_marker_conservation.csv -f {fpkm_filter} -c {conservation_filter}
		"""
rule coexpression:
	input:
		"snakemake/result/11_data_analysis/finally_fpkm.txt",
		"snakemake/result/10_featurecount/all.txt",
		"snakemake/result/11_data_analysis/cluster_diff/unannotated_marker.csv",
		"snakemake/result/11_data_analysis/conservation/ConservedScore.txt"
	output:
		"snakemake/result/11_data_analysis/expression.txt",
		"snakemake/result/11_data_analysis/RowNumber.txt"
	message:
		"the coexpression gene of unannotated gene"
	shell:
		"""
		Rscript snakemake/script/coexpression.R -a {input[0]} -c {input[2]} -i {output[0]} -I {output[1]} && \
		if [[ ! -d "snakemake/result/11_data_analysis/coexpression/" ]]; then mkdir snakemake/result/11_data_analysis/coexpression/;else rm -rf snakemake/result/11_data_analysis/coexpression/*;fi && \
		for i in `cat {output[1]}`; do snakemake/script/coexpressionFDRwithLineNum/coExpressionFDR -p 1 -q 0.05 -n $i -o snakemake/result/11_data_analysis/coexpression/$i.txt {output[0]};done && \
		Rscript snakemake/script/extract_enrich_id.R -a snakemake/result/11_data_analysis/coexpression/ && \
		mkdir snakemake/result/11_data_analysis/coexpression/GO_BP && \
		mkdir snakemake/result/11_data_analysis/coexpression/GO_CC && \
		mkdir snakemake/result/11_data_analysis/coexpression/GO_MF && \
		mkdir snakemake/result/11_data_analysis/coexpression/KEGG && \
		cd snakemake/result/11_data_analysis/coexpression/ && \
		find . -type f -name "TCONS*.txt" > name.txt && \
		cd ../../../.. && \
		for i in `cat snakemake/result/11_data_analysis/coexpression/name.txt`;do snakemake/script/genesetEnrichment/geneSetEnrichment -b {G0file} -g snakemake/script/pathway/{Genome}_c5.go.bp.v7.4.symbols.gmt -s snakemake/result/11_data_analysis/coexpression/$i -o snakemake/result/11_data_analysis/coexpression/GO_BP/$i;done && \

		for i in `cat snakemake/result/11_data_analysis/coexpression/name.txt`;do snakemake/script/genesetEnrichment/geneSetEnrichment -b {G0file} -g snakemake/script/pathway/{Genome}_c5.go.cc.v7.4.symbols.gmt -s snakemake/result/11_data_analysis/coexpression/$i -o snakemake/result/11_data_analysis/coexpression/GO_CC/$i;done && \
		
		for i in `cat snakemake/result/11_data_analysis/coexpression/name.txt`;do snakemake/script/genesetEnrichment/geneSetEnrichment -b {G0file} -g snakemake/script/pathway/{Genome}_c5.go.mf.v7.4.symbols.gmt -s snakemake/result/11_data_analysis/coexpression/$i -o snakemake/result/11_data_analysis/coexpression/GO_MF/$i;done && \
		
		for i in `cat snakemake/result/11_data_analysis/coexpression/name.txt`;do snakemake/script/genesetEnrichment/geneSetEnrichment -b {G0file} -g snakemake/script/pathway/{Genome}_c2.cp.kegg.v7.4.symbols.gmt -s snakemake/result/11_data_analysis/coexpression/$i -o snakemake/result/11_data_analysis/coexpression/KEGG/$i;done
		"""
rule repeatmasker:
	input:
		"snakemake/result/11_data_analysis/unannotated_lncRNA_finally.fa",
		"snakemake/result/11_data_analysis/unannotated_protein_finally.fa",
		"snakemake/result/11_data_analysis/cluster_diff/unannotated_marker.csv",
		"snakemake/result/11_data_analysis/unannotated_lncRNA_finally.gtf",
		"snakemake/result/11_data_analysis/unannotated_protein_finally.gtf"
	output:
		"snakemake/result/11_data_analysis/repeakmasker/unannotated_marker_finally.fa"
	message:
		"unannotated marker gene repeatmasker"
	shell:
		"""
		cat {input[0]} {input[1]} > snakemake/result/11_data_analysis/repeakmasker/unannotated_gene_finallly.fa && \
		cat {input[3]} {input[4]} > snakemake/result/11_data_analysis/repeakmasker/unannotated_gene_finallly.gtf && \
		awk -F, '{{print $1}}' {input[2]} > snakemake/result/11_data_analysis/repeakmasker/unannotated_gene.txt && \
		cat snakemake/result/11_data_analysis/repeakmasker/unannotated_gene.txt |sed 's/\"//g'|uniq > snakemake/result/11_data_analysis/repeakmasker/unannotated_gene_1.txt && \
		sed -i '1d' snakemake/result/11_data_analysis/repeakmasker/unannotated_gene_1.txt && \
		sed -i 's/-/_/g' snakemake/result/11_data_analysis/repeakmasker/unannotated_gene_1.txt && \
		if [[ ! -f "snakemake/result/11_data_analysis/repeakmasker/unannotated_marker_finallly.gtf" ]];then echo "file no exist";else rm -f snakemake/result/11_data_analysis/repeakmasker/unannotated_marker_finallly.gtf;fi && \
		for trs in `cat snakemake/result/11_data_analysis/repeakmasker/unannotated_gene_1.txt`;do grep "$trs" snakemake/result/11_data_analysis/repeakmasker/unannotated_gene_finallly.gtf >> snakemake/result/11_data_analysis/repeakmasker/unannotated_marker_finallly.gtf;done && \
		awk '{{print $10}}' snakemake/result/11_data_analysis/repeakmasker/unannotated_marker_finallly.gtf > snakemake/result/11_data_analysis/repeakmasker/unannotated_transcript_id.txt && \
		cat snakemake/result/11_data_analysis/repeakmasker/unannotated_transcript_id.txt |sed 's/\"//g'|uniq > snakemake/result/11_data_analysis/repeakmasker/unannotated_transcript_id_1.txt && \
		sed 's/;$//' snakemake/result/11_data_analysis/repeakmasker/unannotated_transcript_id_1.txt|uniq > snakemake/result/11_data_analysis/repeakmasker/unannotated_transcript_id_2.txt && \
		python snakemake/script/extract_fasta.py snakemake/result/11_data_analysis/repeakmasker/unannotated_gene_finallly.fa snakemake/result/11_data_analysis/repeakmasker/unannotated_transcript_id_2.txt {output} && \
		if [[ ! -d "snakemake/result/11_data_analysis/repeakmasker/output/" ]]; then mkdir snakemake/result/11_data_analysis/repeakmasker/output/;else echo "file exist";fi && \
		RepeatMasker -species "{specie}" -poly -html -gff -dir snakemake/result/11_data_analysis/repeakmasker/output {output} 
		"""
rule cell_metadata:
	input:
		"snakemake/result/11_data_analysis/scRNA.rds"
	output:
		"snakemake/result/11_data_analysis/cluster/pca_{group_chose}.pdf",
		"snakemake/result/11_data_analysis/cluster/tsne_{group_chose}.pdf",
		"snakemake/result/11_data_analysis/cluster/umap_{group_chose}.pdf"
	shell:
		"Rscript snakemake/script/cell_metadata.R -f {output[0]} -l {output[1]} -i {output[2]} -m {input} -j {group_chose}"
rule celltype:
	input:
		"snakemake/result/11_data_analysis/scRNA.rds",
		"snakemake/result/11_data_analysis/conservation/unannotated_score.txt"

	output:
		"snakemake/result/11_data_analysis/celltype/celltype_singleR.csv",
		"snakemake/result/11_data_analysis/celltype/tSNE_celltype.pdf",
		"snakemake/result/11_data_analysis/celltype/UMAP_celltype.pdf",
		"snakemake/result/11_data_analysis/celltype/scRNA_celltype.rds"
	message:
		"identity cell type"
	shell:
		"Rscript snakemake/script/celltype.R -a {input[0]} -b {SingleR_ref1} -c {SingleR_ref2} -d {output[0]} -e {output[1]} -f {output[2]} -g {input[1]} -m {output[3]}"
rule tissue_celltype_specificity:
	input:
		"snakemake/result/11_data_analysis/finally_fpkm.txt",
		"snakemake/result/11_data_analysis/celltype/scRNA_celltype.rds"    
	output:
		"snakemake/result/11_data_analysis/Tissue_specificity/celltype_specificity.pdf",
		"snakemake/result/11_data_analysis/Tissue_specificity/tissue_specificity.pdf",
		"snakemake/result/11_data_analysis/Tissue_specificity/celltype_specificity_ecdf.pdf",
		"snakemake/result/11_data_analysis/Tissue_specificity/tissue_specificity_ecdf.pdf"
	message:
		"calculate tissue/cell sepcificity score"
	shell:
		"Rscript snakemake/script/tissue_specificity.R -a {input[0]}  -c {input[1]} -i {output[0]} -I {output[1]} -o {output[2]} -O {output[3]} -d {group_chose}" 
rule tissue_specificity_heatmap:
	input:
		"snakemake/result/11_data_analysis/finally_fpkm.txt",
		"snakemake/result/11_data_analysis/celltype/scRNA_celltype.rds"
	output:
		"snakemake/result/11_data_analysis/Tissue_specificity/tissue_protein_coding_heatmap.pdf",
		"snakemake/result/11_data_analysis/Tissue_specificity/tissue_annotated_lncRNA_heatmap.pdf",
		"snakemake/result/11_data_analysis/Tissue_specificity/tissue_unannotated_lncRNA_heatmap.pdf",
		"snakemake/result/11_data_analysis/Tissue_specificity/tissue_unannotated_protein_heatmap.pdf"
	message:
		"plot tissue specificity heatmap"
	shell:
		"Rscript snakemake/script/tissue_specificity_heatmap.R -a {input[0]} -c {input[1]} -i {output[0]} -I {output[1]} -o {output[2]} -O {output[3]} -b {group_chose}"
rule celltype_specificity_heatmap:
	input:
		"snakemake/result/11_data_analysis/finally_fpkm.txt",
		"snakemake/result/11_data_analysis/celltype/scRNA_celltype.rds"
	output:
		"snakemake/result/11_data_analysis/Tissue_specificity/celltype_protein_coding_heatmap.pdf",
		"snakemake/result/11_data_analysis/Tissue_specificity/celltype_annotated_lncRNA_heatmap.pdf",
		"snakemake/result/11_data_analysis/Tissue_specificity/celltype_unannotated_lncRNA_heatmap.pdf",
		"snakemake/result/11_data_analysis/Tissue_specificity/celltype_unannotated_protein_heatmap.pdf"
	message:
		"plot tissue specificity heatmap"
	shell:
		"Rscript snakemake/script/celltype_specificity_heatmap.R -a {input[0]} -c {input[1]} -i {output[0]} -I {output[1]} -o {output[2]} -O {output[3]}"
rule lncRNA_cluster:
	input:
		"snakemake/result/11_data_analysis/finally_fpkm.txt",
		"snakemake/result/11_data_analysis/celltype/scRNA_celltype.rds"
	output:
		"snakemake/result/11_data_analysis/cluster/lncRNA/VariableFeatures.pdf",
		"snakemake/result/11_data_analysis/cluster/lncRNA/pca.pdf",
		"snakemake/result/11_data_analysis/cluster/lncRNA/cell_cluster.csv",
		"snakemake/result/11_data_analysis/cluster/lncRNA/tSNE.pdf",
		"snakemake/result/11_data_analysis/cluster/lncRNA/UMAP.pdf",
		"snakemake/result/11_data_analysis/cluster/lncRNA/top10_marker.pdf",
		"snakemake/result/11_data_analysis/cluster/lncRNA/tSNE_{group_chose}.pdf",
		"snakemake/result/11_data_analysis/cluster/lncRNA/UMAP_{group_chose}.pdf"
	message:
		"cluster data by lncRNA"
	shell:
		"Rscript snakemake/script/lncRNA_cluster.R -a {input[0]} -b {group_chose} -j {resolution} -c {input[1]} -k {PCAnum} -e {output[0]} -f {output[1]} -g {output[2]} -l {output[3]} -i {output[4]} -m {output[5]} -n {output[6]} -o {output[7]}"
rule protein_coding_cluster:
	input:
		"snakemake/result/11_data_analysis/finally_fpkm.txt",
		"snakemake/result/11_data_analysis/celltype/scRNA_celltype.rds"
	output:
		"snakemake/result/11_data_analysis/cluster/protein/VariableFeatures.pdf",
		"snakemake/result/11_data_analysis/cluster/protein/pca.pdf",
		"snakemake/result/11_data_analysis/cluster/protein/cell_cluster.csv",
		"snakemake/result/11_data_analysis/cluster/protein/tSNE.pdf",
		"snakemake/result/11_data_analysis/cluster/protein/UMAP.pdf",
		"snakemake/result/11_data_analysis/cluster/protein/top10_marker.pdf",
		"snakemake/result/11_data_analysis/cluster/protein/tSNE_{group_chose}.pdf",
		"snakemake/result/11_data_analysis/cluster/protein/UMAP_{group_chose}.pdf"
	message:
		"cluster data by protein coding gene"
	shell:
		"Rscript snakemake/script/protein_cluster.R -a {input[0]} -b {group_chose} -j {resolution} -c {input[1]} -k {PCAnum} -e {output[0]} -f {output[1]} -g {output[2]} -l {output[3]} -i {output[4]} -m {output[5]} -n {output[6]} -o {output[7]}"
rule merge_QC_pdf:
	input:
		"snakemake/result/11_data_analysis/QC/vlnplot_before_qc.pdf",
		"snakemake/result/11_data_analysis/QC/vlnplot_after_qc.pdf",
		"snakemake/result/11_data_analysis/cluster/VariableFeatures.pdf"
	output:
		"snakemake/result/11_data_analysis/QC.pdf"
	shell:
		"""
		mutool merge -o {output} {input[0]} {input[1]} {input[2]} && \
		mupdf-x11 {output} &> /dev/null
		"""
rule merge_cluster_pdf:
	input:
		"snakemake/result/11_data_analysis/cluster/pca.pdf",
		"snakemake/result/11_data_analysis/cluster/tSNE.pdf",
		"snakemake/result/11_data_analysis/cluster/UMAP.pdf",
		"snakemake/result/11_data_analysis/cluster_diff/top10_marker.pdf",
		expand("snakemake/result/11_data_analysis/cluster/pca_{group_chose}.pdf",group_chose=group_chose),
		expand("snakemake/result/11_data_analysis/cluster/tsne_{group_chose}.pdf",group_chose=group_chose),
		expand("snakemake/result/11_data_analysis/cluster/umap_{group_chose}.pdf",group_chose=group_chose)
	output:
		"snakemake/result/11_data_analysis/cluster.pdf"
	shell:
		"""
		mutool merge -o {output} {input[0]} {input[1]} {input[2]} {input[3]} {input[4]} {input[5]} && \
		mupdf-x11 {output} &> /dev/null
		"""
rule celltype_pdf:
	input:
		"snakemake/result/11_data_analysis/celltype/tSNE_celltype.pdf",
		"snakemake/result/11_data_analysis/celltype/UMAP_celltype.pdf"
	output:
		"snakemake/result/11_data_analysis/celltype.pdf"
	shell:
		"""
		mutool merge -o {output} {input[0]} {input[1]} && \
		mupdf-x11 {output} &> /dev/null
		"""	
rule tissue_specificity_pdf:
	input:
		"snakemake/result/11_data_analysis/Tissue_specificity/celltype_specificity.pdf",
		"snakemake/result/11_data_analysis/Tissue_specificity/tissue_specificity.pdf",
		"snakemake/result/11_data_analysis/Tissue_specificity/celltype_specificity_ecdf.pdf",
		"snakemake/result/11_data_analysis/Tissue_specificity/tissue_specificity_ecdf.pdf",
		"snakemake/result/11_data_analysis/Tissue_specificity/tissue_protein_coding_heatmap.pdf",
		"snakemake/result/11_data_analysis/Tissue_specificity/tissue_annotated_lncRNA_heatmap.pdf",
		"snakemake/result/11_data_analysis/Tissue_specificity/tissue_unannotated_lncRNA_heatmap.pdf",
		"snakemake/result/11_data_analysis/Tissue_specificity/tissue_unannotated_protein_heatmap.pdf",
		"snakemake/result/11_data_analysis/Tissue_specificity/celltype_protein_coding_heatmap.pdf",
		"snakemake/result/11_data_analysis/Tissue_specificity/celltype_annotated_lncRNA_heatmap.pdf",
		"snakemake/result/11_data_analysis/Tissue_specificity/celltype_unannotated_lncRNA_heatmap.pdf",
		"snakemake/result/11_data_analysis/Tissue_specificity/celltype_unannotated_protein_heatmap.pdf",
		"snakemake/result/11_data_analysis/conservation/gene_conservation_density.pdf",
        "snakemake/result/11_data_analysis/conservation/gene_conservation_ecdf.pdf",
		"snakemake/result/11_data_analysis/cluster_diff/unannotated_marker_vloplot.pdf",
		"snakemake/result/11_data_analysis/cluster_diff/unannotated_marker_featureplot.pdf"
	output:
		"snakemake/result/11_data_analysis/uannotated_gene_feature.pdf"	
	shell:
		"""
		mutool merge -o {output} {input[0]} {input[1]} {input[2]} {input[3]} {input[4]} {input[5]} {input[6]} {input[7]} {input[8]} {input[9]} {input[10]} {input[11]} {input[12]} {input[13]} {input[14]} {input[15]} && \
		mupdf-x11 {output} &> /dev/null
		"""	
rule lncRNA_pdf:        
	input:
		"snakemake/result/11_data_analysis/cluster/lncRNA/VariableFeatures.pdf",
		"snakemake/result/11_data_analysis/cluster/lncRNA/pca.pdf",
		"snakemake/result/11_data_analysis/cluster/lncRNA/tSNE.pdf",
		"snakemake/result/11_data_analysis/cluster/lncRNA/UMAP.pdf",
		"snakemake/result/11_data_analysis/cluster/lncRNA/top10_marker.pdf",
		expand("snakemake/result/11_data_analysis/cluster/lncRNA/tSNE_{group_chose}.pdf",group_chose=group_chose),
		expand("snakemake/result/11_data_analysis/cluster/lncRNA/UMAP_{group_chose}.pdf",group_chose=group_chose)
	output:
		"snakemake/result/11_data_analysis/lncRNA_cluster.pdf"
	shell:
		"""
		mutool merge -o {output} {input[0]} {input[1]} {input[2]} {input[3]} {input[4]} {input[5]} {input[6]} && \
		mupdf-x11 {output} &> /dev/null
		"""
rule protein_pdf:
	input:
		"snakemake/result/11_data_analysis/cluster/protein/VariableFeatures.pdf",
		"snakemake/result/11_data_analysis/cluster/protein/pca.pdf",
		"snakemake/result/11_data_analysis/cluster/protein/tSNE.pdf",
		"snakemake/result/11_data_analysis/cluster/protein/UMAP.pdf",
		"snakemake/result/11_data_analysis/cluster/protein/top10_marker.pdf",
		expand("snakemake/result/11_data_analysis/cluster/protein/tSNE_{group_chose}.pdf",group_chose=group_chose),
		expand("snakemake/result/11_data_analysis/cluster/protein/UMAP_{group_chose}.pdf",group_chose=group_chose)
	output:
		"snakemake/result/11_data_analysis/protein_cluster.pdf"
	shell:
		"""
		mutool merge -o {output} {input[0]} {input[1]} {input[2]} {input[3]} {input[4]} {input[5]} {input[6]} && \
		mupdf-x11 {output} &> /dev/null
		"""	

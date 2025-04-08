Genome = config["Genome"]
G0file = config["G0file"]
Gtfile = config["gtf"]
lncRNA_id = config["gene_id"]
gene_set = config["gene_set"]
OutputDir = config["outputDir"]
rule all:
	input:
		expand("{OutputDir}/snakemake/result/11_data_analysis/unannotated_transcript/RowNumber.txt",OutputDir=OutputDir),
		expand("{OutputDir}/snakemake/result/11_data_analysis/unannotated_transcript/{lncRNA_id}/{lncRNA_id}_structure.pdf",OutputDir=OutputDir,lncRNA_id=lncRNA_id),
		expand("{OutputDir}/snakemake/result/11_data_analysis/unannotated_transcript/{lncRNA_id}/{lncRNA_id}.pdf",OutputDir=OutputDir,lncRNA_id=lncRNA_id)
rule coexpression:
	input:
		"{OutputDir}/snakemake/result/11_data_analysis/finally_tpm.txt"
	output:
		"{OutputDir}/snakemake/result/11_data_analysis/unannotated_transcript/expression.txt",
		"{OutputDir}/snakemake/result/11_data_analysis/unannotated_transcript/RowNumber.txt"
	conda:
		"env/R_Library.yaml"
	message:
		"the coexpression gene of unannotated gene"
	shell:
		"""
		Rscript snakemake/script/coexpression_2.R -a {input[0]} -c {lncRNA_id} -i {output[0]} -I {output[1]} && \
		if [[ ! -d "{OutputDir}/snakemake/result/11_data_analysis/unannotated_transcript/{lncRNA_id}/" ]]; then mkdir {OutputDir}/snakemake/result/11_data_analysis/unannotated_transcript/{lncRNA_id}/;fi && \
		for i in `cat {output[1]}`; do snakemake/script/coexpressionFDRwithLineNum/coExpressionFDR -p 1 -q 0.05 -n $i -o {OutputDir}/snakemake/result/11_data_analysis/unannotated_transcript/{lncRNA_id}/$i.txt {output[0]};done && \
		Rscript snakemake/script/extract_enrich_id.R -a {OutputDir}/snakemake/result/11_data_analysis/unannotated_transcript/{lncRNA_id}/ && \
		snakemake/script/genesetEnrichment/geneSetEnrichment -b {G0file} -g {gene_set} -s {OutputDir}/snakemake/result/11_data_analysis/unannotated_transcript/{lncRNA_id}/{lncRNA_id}.CoexGene -o {OutputDir}/snakemake/result/11_data_analysis/unannotated_transcript/{lncRNA_id}/{lncRNA_id}.file && \
		Rscript snakemake/script/Function_plot.R -i {OutputDir}/snakemake/result/11_data_analysis/unannotated_transcript/{lncRNA_id}/{lncRNA_id}.file -o {OutputDir}/snakemake/result/11_data_analysis/unannotated_transcript/{lncRNA_id}/{lncRNA_id}_function.pdf
		"""
rule transcript_structure:
	input:
		"{OutputDir}/snakemake/result/9_unannotated_lncRNA/unLncRNA.unannotated_lncRNA_4.gtf.tmap",
		"{OutputDir}/snakemake/result/8_unannotated_protein_coding/unMRNA.unannotated_protein.gtf.tmap",
		"{OutputDir}/snakemake/result/11_data_analysis/unannotated_lncRNA_finally.gtf",
		"{OutputDir}/snakemake/result/11_data_analysis/unannotated_protein_finally.gtf",
		"{OutputDir}/snakemake/result/11_data_analysis/celltype/scRNA_celltype.rds"
	output:
		"{OutputDir}/snakemake/result/11_data_analysis/unannotated_transcript/{lncRNA_id}/{lncRNA_id}_structure.pdf"
	conda:
		"env/R_Library.yaml"
	shell:
		"""
		Rscript snakemake/script/transcript_structure.R -a {input[0]} -b {input[1]} -c {input[2]} -d {input[3]} -e {Gtfile} -f {lncRNA_id} -g {input[4]} -o {output}
		"""
rule merge_pdf:
	input:
		"{OutputDir}/snakemake/result/11_data_analysis/unannotated_transcript/{lncRNA_id}/{lncRNA_id}_structure.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/unannotated_transcript/RowNumber.txt"
	output:
		"{OutputDir}/snakemake/result/11_data_analysis/unannotated_transcript/{lncRNA_id}/{lncRNA_id}.pdf"
	shell:
		"""
		mutool merge -o {output} {OutputDir}/snakemake/result/11_data_analysis/unannotated_transcript/{lncRNA_id}/{lncRNA_id}_function.pdf {input[0]} && \
		mupdf-x11 {output} &> /dev/null
		"""

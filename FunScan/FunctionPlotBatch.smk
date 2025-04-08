Gtfile = config["GTFfile"]
gene_list=config["gene_list"]
Genome=config["Genome"]
G0file=config["G0file"]
OutputDir=config["outputdir"]
rule all:
	input:
		expand("{OutputDir}/snakemake/result/11_data_analysis/expression.txt",OutputDir=OutputDir),
		expand("{OutputDir}/snakemake/result/11_data_analysis/transcript_scture",OutputDir=OutputDir)
rule coexpression:
	input:
		"{OutputDir}/snakemake/result/11_data_analysis/finally_tpm.txt"
	output:
		"{OutputDir}/snakemake/result/11_data_analysis/expression.txt",
		"{OutputDir}/snakemake/result/11_data_analysis/coexpression/Row.Number"
	conda:
		"env/R_Library.yaml"
	message:
		"the coexpression gene of unannotated gene"
	shell:
		"""
		Rscript snakemake/script/coexpression.R -a {input[0]} -c {gene_list} -i {output[0]} -I {output[1]} && \
		if [[  -d "{wildcards.OutputDir}/snakemake/result/11_data_analysis/coexpression/" ]]; then rm -rf {wildcards.OutputDir}/snakemake/result/11_data_analysis/coexpression/*.txt; \
				rm -rf {wildcards.OutputDir}/snakemake/result/11_data_analysis/coexpression/*.CoexGene;fi && \
		for i in `cat {output[1]}`; do snakemake/script/coexpressionFDRwithLineNum/coExpressionFDR -p 1 -q 0.05 -n $i -o {wildcards.OutputDir}/snakemake/result/11_data_analysis/coexpression/$i.txt {output[0]};done && \
		Rscript snakemake/script/extract_enrich_id.R -a {wildcards.OutputDir}/snakemake/result/11_data_analysis/coexpression/ && \
		ls {wildcards.OutputDir}/snakemake/result/11_data_analysis/coexpression/*CoexGene > {wildcards.OutputDir}/snakemake/result/11_data_analysis/coexpression/name.list && \
		for i in `cat {wildcards.OutputDir}/snakemake/result/11_data_analysis/coexpression/name.list`;do file_name=$(basename $i);file_name2=${{file_name%%.*}}; \
			if [[ ! -d "{wildcards.OutputDir}/snakemake/result/11_data_analysis/coexpression/$file_name2" ]]; then mkdir {wildcards.OutputDir}/snakemake/result/11_data_analysis/coexpression/$file_name2/;fi; \
			snakemake/script/genesetEnrichment/geneSetEnrichment -b {G0file} -g snakemake/script/pathway/{Genome}_c5.go.bp.v7.4.symbols.gmt \
			-s {wildcards.OutputDir}/snakemake/result/11_data_analysis/coexpression/$file_name -o {wildcards.OutputDir}/snakemake/result/11_data_analysis/coexpression/$file_name2/GO_BP.txt;done && \

		for i in `cat {wildcards.OutputDir}/snakemake/result/11_data_analysis/coexpression/name.list`;do file_name=$(basename $i);file_name2=${{file_name%%.*}}; \
			if [[ ! -d "{wildcards.OutputDir}/snakemake/result/11_data_analysis/coexpression/$file_name2" ]]; then mkdir {wildcards.OutputDir}/snakemake/result/11_data_analysis/coexpression/$file_name2/;fi; \
			snakemake/script/genesetEnrichment/geneSetEnrichment -b {G0file} -g snakemake/script/pathway/{Genome}_c5.go.cc.v7.4.symbols.gmt \
			-s {wildcards.OutputDir}/snakemake/result/11_data_analysis/coexpression/$file_name -o {wildcards.OutputDir}/snakemake/result/11_data_analysis/coexpression/$file_name2/GO_CC.txt;done && \
		
		for i in `cat {wildcards.OutputDir}/snakemake/result/11_data_analysis/coexpression/name.list`;do file_name=$(basename $i);file_name2=${{file_name%%.*}}; \
			if [[ ! -d "{wildcards.OutputDir}/snakemake/result/11_data_analysis/coexpression/$file_name2" ]]; then mkdir {wildcards.OutputDir}/snakemake/result/11_data_analysis/coexpression/$file_name2/;fi; \
			snakemake/script/genesetEnrichment/geneSetEnrichment -b {G0file} -g snakemake/script/pathway/{Genome}_c5.go.mf.v7.4.symbols.gmt \
			-s {wildcards.OutputDir}/snakemake/result/11_data_analysis/coexpression/$file_name -o {wildcards.OutputDir}/snakemake/result/11_data_analysis/coexpression/$file_name2/GO_MF.txt;done && \
		
		for i in `cat {wildcards.OutputDir}/snakemake/result/11_data_analysis/coexpression/name.list`;do file_name=$(basename $i);file_name2=${{file_name%%.*}}; \
			if [[ ! -d "{wildcards.OutputDir}/snakemake/result/11_data_analysis/coexpression/$file_name2" ]]; then mkdir {wildcards.OutputDir}/snakemake/result/11_data_analysis/coexpression/$file_name2/;fi; \
			snakemake/script/genesetEnrichment/geneSetEnrichment -b {G0file} -g snakemake/script/pathway/{Genome}_c2.cp.kegg.v7.4.symbols.gmt \
			-s {wildcards.OutputDir}/snakemake/result/11_data_analysis/coexpression/$file_name -o {wildcards.OutputDir}/snakemake/result/11_data_analysis/coexpression/$file_name2/KEGG.txt;done	
		"""
rule transcripts_features:
	input:
		"{OutputDir}/snakemake/result/9_unannotated_lncRNA/unLncRNA.unannotated_lncRNA_4.gtf.tmap",
		"{OutputDir}/snakemake/result/8_unannotated_protein_coding/unMRNA.unannotated_protein.gtf.tmap",
		"{OutputDir}/snakemake/result/11_data_analysis/unannotated_lncRNA_finally.gtf",
		"{OutputDir}/snakemake/result/11_data_analysis/unannotated_protein_finally.gtf",
		"{OutputDir}/snakemake/result/11_data_analysis/celltype/scRNA_celltype.rds"
	output:
		directory("{OutputDir}/snakemake/result/11_data_analysis/transcript_scture")
	conda:
		"env/R_Library.yaml"
	shell:
		"""
		if [[ ! -d "{OutputDir}/snakemake/result/11_data_analysis/transcript_scture/" ]]; \
			then mkdir {OutputDir}/snakemake/result/11_data_analysis/transcript_scture;fi && \
		for i in `cat {gene_list}`;do Rscript snakemake/script/transcript_structure.R -a {input[0]} \
			-b {input[1]} -c {input[2]} \
			-d {input[3]} -e {Gtfile} -f $i -g {input[4]} \
			-o {output}/$i-structure.pdf;done
		"""
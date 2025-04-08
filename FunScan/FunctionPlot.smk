lncRNAID = config["lncRNA_id"]
outputFile = config["outputfile"]
OutputDir=config["outputdir"]
rule all:
	input:
		expand("{OutputDir}/{outputFile}",outputFile=outputFile,OutputDir=OutputDir)
rule function_plot:
	input:
		expand("{OutputDir}/snakemake/result/11_data_analysis/finally_tpm.txt",OutputDir=OutputDir)
	output:
		expand("{OutputDir}/{outputFile}",outputFile=outputFile,OutputDir=OutputDir)
	conda:
		"env/R_Library.yaml"
	shell:
		"""
		Rscript snakemake/script/Function_plot2.R -i {lncRNAID} -o {OutputDir}/{outputFile} -d {OutputDir} && \
		mupdf-x11 {OutputDir}/{outputFile}
		"""
rule makr_report:
	input:
		"snakemake/result/11_data_analysis/celltype/scRNA_celltype.rds"
	output:
		"snakemake/result/FunScan_report.html"
	conda:
		"env/R_Library.yaml"
	params:
		finally_tpm="snakemake/result/11_data_analysis/finally_tpm.txt",
		all_id="snakemake/result/10_featurecount/all.txt"
	script:
		"script/FunScan.Rmd"
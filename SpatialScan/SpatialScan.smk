rawData = config["rawData"]
Celltype = config["Celltype"]
Gene = config["Gene"]
OutputDir=config["outputDir"]
SpatialData=config["rdsData"]
rule all:
	input:
		expand("{OutputDir}/snakemake/result/11_data_analysis/Spatial/SpatialScan_report.html",OutputDir=OutputDir),
		expand("snakemake/result/SpatialScan_report.html")
rule celltrek:
	input:
		"{OutputDir}/snakemake/result/11_data_analysis/celltype/scRNA_celltype.rds"
	output:
		"{OutputDir}/snakemake/result/11_data_analysis/Spatial/Spatial.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/Spatial/sc_st_coembedding.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/Spatial/sc_st_SpatialDimPlot.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/Spatial/heatmap.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/Spatial/celltype_featureplot_module.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/Spatial/celltype_spatialplot_module.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/Spatial/gene_spatialplot.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/Spatial/gene_spatialplot_module.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/Spatial/nCount_Spatial.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/Spatial/umap.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/Spatial/top10_genes.pdf"
	conda:
		"env/R_Library.yaml"
	shell:
		"Rscript snakemake/script/celltrek.R -a {input} -b {ST} -c {output[0]} -d {output[1]} -e {output[2]} -f {Celltype} \
		-g {output[3]} -m {output[4]} -i {output[5]} -j {Gene} -k {output[6]} -l {output[7]}"
if config["DataType"] == 'true':
	rule make_report:
		input:
			"{OutputDir}/snakemake/result/11_data_analysis/celltype/scRNA_celltype.rds"
		output:
			"{OutputDir}/snakemake/result/11_data_analysis/Spatial/SpatialScan_report.html"
		conda:
			"env/R_Library.yaml"
		params:
			Outputdir={OutputDir},
			rawdata={rawData},
			scRNAData="{OutputDir}/snakemake/result/11_data_analysis/celltype/scRNA_celltype.rds",
			STpdf="{OutputDir}/snakemake/result/11_data_analysis/Spatial/Spatial.pdf",
			Coembedding="{OutputDir}/snakemake/result/11_data_analysis/Spatial/sc_st_coembedding.pdf",
			SpatialDimPlot="{OutputDir}/snakemake/result/11_data_analysis/Spatial/sc_st_SpatialDimPlot.pdf",
			Celltype={Celltype},
			Heatmap="{OutputDir}/snakemake/result/11_data_analysis/Spatial/heatmap.pdf",
			CelltypeFeatureplot="{OutputDir}/snakemake/result/11_data_analysis/Spatial/celltype_featureplot_module.pdf",
			CelltypeSpatialplot="{OutputDir}/snakemake/result/11_data_analysis/Spatial/celltype_spatialplot_module.pdf",
			Gene={Gene},
			GeneSpatial="{OutputDir}/snakemake/result/11_data_analysis/Spatial/gene_spatialplot.pdf",
			GeneSpatialModule="{OutputDir}/snakemake/result/11_data_analysis/Spatial/gene_spatialplot_module.pdf",
			nCountSpatialBefore="{OutputDir}/snakemake/result/11_data_analysis/Spatial/nCountSpatialBefore.pdf",
			nFeatureSpatialBefore="{OutputDir}/snakemake/result/11_data_analysis/Spatial/nFeatureSpatialBefore.pdf",
			percentMTbefore="{OutputDir}/snakemake/result/11_data_analysis/Spatial/percentMTbefore.pdf",
			nCountSpatialAfter="{OutputDir}/snakemake/result/11_data_analysis/Spatial/nCountSpatialAfter.pdf",
			nFeatureSpatialAfter="{OutputDir}/snakemake/result/11_data_analysis/Spatial/nFeatureSpatialAfter.pdf",
			percentMTAfter="{OutputDir}/snakemake/result/11_data_analysis/Spatial/percentMTAfter.pdf",
			minCount=config["minCount"],
			minFeature=config["minFeature"],
			maxMT=config["maxMT"],
			umap="{OutputDir}/snakemake/result/11_data_analysis/Spatial/umap.pdf",
			top10genes="{OutputDir}/snakemake/result/11_data_analysis/Spatial/top10_genes.pdf"
		script:
			"script/SpatialScan.Rmd"
else:
	rule make_report:
		input:
			"{OutputDir}/snakemake/result/11_data_analysis/celltype/scRNA_celltype.rds"
		output:
			"{OutputDir}/snakemake/result/11_data_analysis/Spatial/SpatialScan_report.html"
		conda:
			"env/R_Library.yaml"
		params:
			Outputdir={OutputDir},
			Spatialdata={SpatialData},
			scRNAData="{OutputDir}/snakemake/result/11_data_analysis/celltype/scRNA_celltype.rds",
			STpdf="{OutputDir}/snakemake/result/11_data_analysis/Spatial/Spatial.pdf",
			Coembedding="{OutputDir}/snakemake/result/11_data_analysis/Spatial/sc_st_coembedding.pdf",
			SpatialDimPlot="{OutputDir}/snakemake/result/11_data_analysis/Spatial/sc_st_SpatialDimPlot.pdf",
			Celltype={Celltype},
			Heatmap="{OutputDir}/snakemake/result/11_data_analysis/Spatial/heatmap.pdf",
			CelltypeFeatureplot="{OutputDir}/snakemake/result/11_data_analysis/Spatial/celltype_featureplot_module.pdf",
			CelltypeSpatialplot="{OutputDir}/snakemake/result/11_data_analysis/Spatial/celltype_spatialplot_module.pdf",
			Gene={Gene},
			GeneSpatial="{OutputDir}/snakemake/result/11_data_analysis/Spatial/gene_spatialplot.pdf",
			GeneSpatialModule="{OutputDir}/snakemake/result/11_data_analysis/Spatial/gene_spatialplot_module.pdf"
		script:
			"script/SpatialScan2.Rmd"
rule makr_report2:
	input:
		expand("{OutputDir}/snakemake/result/11_data_analysis/Spatial/SpatialScan_report.html",OutputDir=OutputDir)
	output:
		"snakemake/result/SpatialScan_report.html"
	shell:
		"cp -r {OutputDir}/snakemake/result/11_data_analysis/Spatial/SpatialScan_report.html {output}"
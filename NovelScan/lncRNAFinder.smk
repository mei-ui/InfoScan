#configfile: "snakemake/config/lncRNAFinder_config.yml"
GTF = config["gtf"]
FA = config["fa"]
GENCODE_protein_gtf = config["protein"]
GENCODE_lncRNA_gene_gtf = config["lncRNA"]
ERCC = config["ercc"]
CPATlogitModel = config["CPATlogitModel"]
CPATHexamer = config["CPATHexamer"]
CPATcutoff = config["CPATcutoff"]
specie = config["genome"]
model = config["model"]
len_filter = config["length_filter"] 
fpkm_threshold = config["fpkm_filter"]
exon_filter=config["exon_filter"]
phastcons = config["phastcons"]
GENCODE_pseudogene = config["GENCODE_pseudogene"]
GENCODE_other_rna = config["GENCODE_other_rna"]
SAMPLES = config["sample"]
OutputDir = config["outdir"]
workDir = config["appDir"]
#(SAMPLES,) = glob_wildcards("{OutputDir}/snakemake/result/3_HISAT2_aligned/{sample}.bam")
rule all:
	input:
		expand("{OutputDir}/snakemake/result/6_protein_coding/protein_coding.gtf",OutputDir = OutputDir),
		expand("{OutputDir}/snakemake/result/7_annotated_lncRNA/annotated_lncRNA.gtf",OutputDir = OutputDir),
		expand("{OutputDir}/snakemake/result/8_unannotated_protein_coding/overlap/overlap.txt",OutputDir = OutputDir),
		expand("{OutputDir}/snakemake/result/8_unannotated_protein_coding/unannotated_protein.gtf",OutputDir = OutputDir),
		expand("{OutputDir}/snakemake/result/9_unannotated_lncRNA/overlap/overlap.txt",OutputDir = OutputDir),
		expand("{OutputDir}/snakemake/result/9_unannotated_lncRNA/unannotated_lncRNA.gtf",OutputDir = OutputDir),
		expand("{OutputDir}/snakemake/result/10_featurecount/unannotated_transcript.gtf",OutputDir = OutputDir),
		expand("{OutputDir}/snakemake/result/10_featurecount/fpkm_matrix.txt",OutputDir = OutputDir),
		expand("{OutputDir}/snakemake/result/11_data_analysis/finally_fpkm.txt",OutputDir = OutputDir),
		expand("{OutputDir}/snakemake/result/11_data_analysis/gene_FPKM_distrubution.pdf",OutputDir = OutputDir),
		expand("{OutputDir}/snakemake/result/11_data_analysis/gene-length-distrubution.pdf",OutputDir = OutputDir),
		expand("{OutputDir}/snakemake/result/11_data_analysis/exon_number.pdf",OutputDir = OutputDir),
		expand("{OutputDir}/snakemake/result/11_data_analysis/unannotated_lncRNA_finally.fa",OutputDir = OutputDir),
		expand("{OutputDir}/snakemake/result/11_data_analysis/Sequence_features/features.txt",OutputDir = OutputDir),
		expand("{OutputDir}/snakemake/result/11_data_analysis/conservation/gene_conservation_density.pdf",OutputDir = OutputDir),
		expand("{OutputDir}/snakemake/result/12_pseudogene/pseudogene.gtf",OutputDir = OutputDir),
		expand("{OutputDir}/snakemake/result/13_other_rnas/other_rna.gtf",OutputDir = OutputDir),
		expand("{OutputDir}/snakemake/result/NovelScan_report_pdf.html",OutputDir = OutputDir),
		expand("{OutputDir}/snakemake/result/NovelScan_report.html",OutputDir = OutputDir)

rule protein_coding:
	input:
		"{OutputDir}/snakemake/result/5_stringtie_merge/stringtie_merged.gtf",
		"{OutputDir}/snakemake/result/5_stringtie_merge/stringtie_merged_transcript.gtf"
	output:
		temp("{OutputDir}/snakemake/result/6_protein_coding/protein_coding_1.gtf"),
		"{OutputDir}/snakemake/result/6_protein_coding/protein_name_1.txt",
		"{OutputDir}/snakemake/result/6_protein_coding/protein_coding.gtf"
	message:
		"identity annotated protein coding gene"
	conda:
		"env/InfoScan.yaml"
	params:
		workDir = config["appDir"]
	shell:
		"""
		bedtools intersect -a {GENCODE_protein_gtf} -b {input[1]} -f 0.9 > {output[0]} && \
		{params.workDir}/script/extract_transcript_name.sh {output[0]} {output[1]} && \
		python {params.workDir}/script/extract_gtf.py -i {output[1]} -f {GENCODE_protein_gtf} -o {output[2]}
		"""
rule annotated_lncRNA:
	input:
		"{OutputDir}/snakemake/result/5_stringtie_merge/stringtie_merged.gtf",
		"{OutputDir}/snakemake/result/5_stringtie_merge/stringtie_merged_transcript.gtf"
	output:
		temp("{OutputDir}/snakemake/result/7_annotated_lncRNA/annotated_lncRNA_1.gtf"),
		"{OutputDir}/snakemake/result/7_annotated_lncRNA/annotated_lncRNA_name_1.txt",
		"{OutputDir}/snakemake/result/7_annotated_lncRNA/annotated_lncRNA.gtf"
	message:
		"identity annotated lncRNA"
	conda:
		"env/InfoScan.yaml"
	params:
		workDir=config["appDir"]
	shell:
		"""
		bedtools intersect -a {GENCODE_lncRNA_gene_gtf} -b {input[1]} -f 0.9 > {output[0]} && \
		{params.workDir}/script/extract_transcript_name.sh {output[0]} {output[1]} && \
		python {params.workDir}/script/extract_gtf.py -i {output[1]} -f {GENCODE_lncRNA_gene_gtf} -o {output[2]}
		"""
rule pseudogene:
	input:
		"{OutputDir}/snakemake/result/5_stringtie_merge/stringtie_merged.gtf",
		"{OutputDir}/snakemake/result/5_stringtie_merge/stringtie_merged_transcript.gtf"
	output:
		temp("{OutputDir}/snakemake/result/12_pseudogene/pseudo_1.gtf"),
		"{OutputDir}/snakemake/result/12_pseudogene/pseudo_name.txt",
		"{OutputDir}/snakemake/result/12_pseudogene/pseudogene.gtf"
	message:
		"identity pseudo gene"
	conda:
		"env/InfoScan.yaml"
	params:
		workDir=config["appDir"]
	shell:
		"""
		bedtools intersect -a {GENCODE_pseudogene} -b {input[1]} -f 0.9 > {output[0]} && \
		{params.workDir}/script/extract_transcript_name.sh {output[0]} {output[1]} && \
		python {params.workDir}/script/extract_gtf.py -i {output[1]} -f {GENCODE_pseudogene} -o {output[2]}
		"""	
rule other_rna:
	input:
		"{OutputDir}/snakemake/result/5_stringtie_merge/stringtie_merged.gtf",
		"{OutputDir}/snakemake/result/5_stringtie_merge/stringtie_merged_transcript.gtf"
	output:
		temp("{OutputDir}/snakemake/result/13_other_rnas/other_rna_1.gtf"),
		"{OutputDir}/snakemake/result/13_other_rnas/other_rna_name.txt",
		"{OutputDir}/snakemake/result/13_other_rnas/other_rna.gtf"
	message:
		"identity other rnas"
	params:
		workDir=config["appDir"]
	conda:
		"env/InfoScan.yaml"
	shell:
		"""
		bedtools intersect -a {GENCODE_other_rna} -b {input[1]} -f 0.9 > {output[0]} && \
		{params.workDir}/script/extract_transcript_name.sh {output[0]} {output[1]} && \
		python {params.workDir}/script/extract_gtf.py -i {output[1]} -f {GENCODE_other_rna} -o {output[2]}
		"""	
rule unannotated_lncRNA:
	input:
		"{OutputDir}/snakemake/result/5_stringtie_merge/stringtie_merged.gtf",
		"{OutputDir}/snakemake/result/5_stringtie_merge/stringtie_merged_transcript.gtf"
	output:
		"{OutputDir}/snakemake/result/9_unannotated_lncRNA/cufcomp.stringtie_merged.gtf.tmap",
		lncScan_svm="{OutputDir}/snakemake/result/9_unannotated_lncRNA/lncScan/lncScan_svm.txt",
		lncScan_XGB="{OutputDir}/snakemake/result/9_unannotated_lncRNA/lncScan/lncScan_XGB.txt",
		CPAT="{OutputDir}/snakemake/result/9_unannotated_lncRNA/CPAT/CPAT.txt",
		lncRNA="{OutputDir}/snakemake/result/9_unannotated_lncRNA/lncRNA_2.gtf",
		fa="{OutputDir}/snakemake/result/9_unannotated_lncRNA/lncRNA_1.fa"
	conda:
		"env/InfoScan.yaml"
	message:
		"identity unannotated lncRNA"
	params:
		workDir=config["appDir"]
	shell:
		"""
		gffcompare -r {GTF} -o {OutputDir}/snakemake/result/9_unannotated_lncRNA/cufcomp {input[0]} && \
		mv {OutputDir}/snakemake/result/5_stringtie_merge/cufcomp.stringtie_merged.gtf.tmap  {output[0]} && \
		python {params.workDir}/script/transcript_filter.py {OutputDir}/snakemake/result/9_unannotated_lncRNA/cufcomp.annotated.gtf {output[lncRNA]} --length {len_filter} && \
		awk '{{print $10}}' {output[lncRNA]} > {OutputDir}/snakemake/result/9_unannotated_lncRNA/lncRNA_1.txt && \
		sed -i "" 's/"//g; s/;//g' {OutputDir}/snakemake/result/9_unannotated_lncRNA/lncRNA_1.txt && \
		gffread -w {output[fa]} -g {FA} {output[lncRNA]} && \
		python {params.workDir}/script/GtfToBed.py {output[lncRNA]} {OutputDir}/snakemake/result/9_unannotated_lncRNA/lncRNA_2.bed && \
		{params.workDir}/script/bigWigAverageOverBed {phastcons} {OutputDir}/snakemake/result/9_unannotated_lncRNA/lncRNA_2.bed {OutputDir}/snakemake/result/9_unannotated_lncRNA/lncRNA_2_conservation.txt && \
		python {params.workDir}/script/Feature_extract.py -a {output[fa]} -d {OutputDir}/snakemake/result/9_unannotated_lncRNA/lncRNA_2_conservation.txt -o {OutputDir}/snakemake/result/9_unannotated_lncRNA/feature.txt && \

		{params.workDir}/script/lncScan.R {specie} {OutputDir}/snakemake/result/9_unannotated_lncRNA/feature.txt {OutputDir}/snakemake/result/9_unannotated_lncRNA/lncScan/result_svm.txt {params.workDir}/script svm && \
		cat {OutputDir}/snakemake/result/9_unannotated_lncRNA/lncScan/result_svm.txt | grep 'non_coding' | cut -f 1 > {output[lncScan_svm]} && \

		{params.workDir}/script/lncScan.R {specie} {OutputDir}/snakemake/result/9_unannotated_lncRNA/feature.txt {OutputDir}/snakemake/result/9_unannotated_lncRNA/lncScan/result_XGB.txt {params.workDir}/script XGB && \
		cat {OutputDir}/snakemake/result/9_unannotated_lncRNA/lncScan/result_XGB.txt | grep 'non_coding' | cut -f 1 > {output[lncScan_XGB]} && \
		
		{params.workDir}/CPAT/bin/cpat.py -g  {output[fa]} --antisense -d {CPATlogitModel} -x {CPATHexamer} -o {OutputDir}/snakemake/result/9_unannotated_lncRNA/CPAT/CPAToutput && \
		Rscript {params.workDir}/script/cpat_noncoding.R -i {OutputDir}/snakemake/result/9_unannotated_lncRNA/CPAT/CPAToutput.ORF_prob.tsv -p {CPATcutoff} -I {OutputDir}/snakemake/result/9_unannotated_lncRNA/CPAT/CPAToutput.no_ORF.txt -o {output[CPAT]}
		"""
rule lnc_overlap:
	input:
		lncScan_svm="{OutputDir}/snakemake/result/9_unannotated_lncRNA/lncScan/lncScan_svm.txt",
		CPAT="{OutputDir}/snakemake/result/9_unannotated_lncRNA/CPAT/CPAT.txt",
		lncScan_XGB="{OutputDir}/snakemake/result/9_unannotated_lncRNA/lncScan/lncScan_XGB.txt"
	output:
		veen="{OutputDir}/snakemake/result/9_unannotated_lncRNA/overlap/veen-pdf.pdf",
		overlap="{OutputDir}/snakemake/result/9_unannotated_lncRNA/overlap/overlap.txt"
	params:
		workDir=config["appDir"]
	shell:
		"python {params.workDir}/script/overlap.py -a {input[lncScan_svm]} -b {input[CPAT]} -c {input[lncScan_XGB]} -d {output[veen]} -e {output[overlap]}"
rule unannotated_lncRNA_2:
	input:
		"{OutputDir}/snakemake/result/9_unannotated_lncRNA/lncRNA_2.gtf",
		"{OutputDir}/snakemake/result/9_unannotated_lncRNA/overlap/overlap.txt"
	output:	
		"{OutputDir}/snakemake/result/9_unannotated_lncRNA/unannotated_lncRNA.gtf"
	params:
		workDir=config["appDir"]
	shell:
		"""
		python {params.workDir}/script/extract_gtf.py -f {input[0]} -i {input[1]} -o {OutputDir}/snakemake/result/9_unannotated_lncRNA/lncRNA_3.gtf && \
		gffcompare -r {GENCODE_protein_gtf} -o {OutputDir}/snakemake/result/9_unannotated_lncRNA/unLncRNA {OutputDir}/snakemake/result/9_unannotated_lncRNA/lncRNA_3.gtf && \
		awk '{{if ($3=="u" || $3=="i"|| $3=="j"|| $3=="x" ||$3=="o") print $5}}' {OutputDir}/snakemake/result/9_unannotated_lncRNA/unLncRNA.lncRNA_3.gtf.tmap|uniq > {OutputDir}/snakemake/result/9_unannotated_lncRNA/unlncrna.txt && \
		python {params.workDir}/script/extract_gtf.py -f {OutputDir}/snakemake/result/9_unannotated_lncRNA/lncRNA_3.gtf -i {OutputDir}/snakemake/result/9_unannotated_lncRNA/unlncrna.txt -o {output[0]}
		"""
rule unannotated_protein:
	input:
		"{OutputDir}/snakemake/result/5_stringtie_merge/stringtie_merged.gtf",
		"{OutputDir}/snakemake/result/5_stringtie_merge/stringtie_merged_transcript.gtf",
		"{OutputDir}/snakemake/result/9_unannotated_lncRNA/cufcomp.stringtie_merged.gtf.tmap"
	output:
		"{OutputDir}/snakemake/result/8_unannotated_protein_coding/cufcomp.stringtie_merged.gtf.tmap",
		gtf="{OutputDir}/snakemake/result/8_unannotated_protein_coding/filter_2.gtf",
		fa="{OutputDir}/snakemake/result/8_unannotated_protein_coding/filter_2.fa",
		lncScan_svm="{OutputDir}/snakemake/result/8_unannotated_protein_coding/lncScan/lncScan_svm.txt",
		lncScan_XGB="{OutputDir}/snakemake/result/8_unannotated_protein_coding/lncScan/lncScan_XGB.txt",
		CPAT="{OutputDir}/snakemake/result/8_unannotated_protein_coding/CPAT/CPAT.txt"
#		unannotated_protein_1_fa="{OutputDir}/snakemake/result/8_unannotated_protein_coding/unannotated_protein_1.fa"
	message:
		"identity unannotated protein"
	params:
		workDir=config["appDir"]
	shell:
		"""
		cp {input[2]} {output[0]} && \
		python {params.workDir}/script/transcript_filter.py {OutputDir}/snakemake/result/9_unannotated_lncRNA/cufcomp.annotated.gtf {output[gtf]} --length {len_filter} && \
		awk '{{print $10}}' {output[gtf]} > {OutputDir}/snakemake/result/8_unannotated_protein_coding/filter_2.txt && \
		sed -i "" 's/"//g; s/;//g' {OutputDir}/snakemake/result/8_unannotated_protein_coding/filter_2.txt && \
		gffread -w {output[fa]} -g {FA} {output[gtf]} && \
		python {params.workDir}/script/GtfToBed.py {output[gtf]} {OutputDir}/snakemake/result/8_unannotated_protein_coding/filter_2.bed && \
		{params.workDir}/script/bigWigAverageOverBed {phastcons} {OutputDir}/snakemake/result/8_unannotated_protein_coding/filter_2.bed {OutputDir}/snakemake/result/8_unannotated_protein_coding/filter_2_conservation.txt && \
		python {params.workDir}/script/Feature_extract.py -a {output[fa]} -d {OutputDir}/snakemake/result/8_unannotated_protein_coding/filter_2_conservation.txt -o {OutputDir}/snakemake/result/8_unannotated_protein_coding/feature.txt && \

		{params.workDir}/script/lncScan.R {specie} {OutputDir}/snakemake/result/8_unannotated_protein_coding/feature.txt {OutputDir}/snakemake/result/8_unannotated_protein_coding/lncScan/result_svm.txt {params.workDir}/script svm && \
		awk -F'\t' 'NR>1 && $2=="coding" {{print $1}}' {OutputDir}/snakemake/result/8_unannotated_protein_coding/lncScan/result_svm.txt > {output[lncScan_svm]} && \

		{params.workDir}/script/lncScan.R {specie} {OutputDir}/snakemake/result/8_unannotated_protein_coding/feature.txt {OutputDir}/snakemake/result/8_unannotated_protein_coding/lncScan/result_XGB.txt {params.workDir}/script XGB && \
		awk -F'\t' 'NR>1 && $2=="coding" {{print $1}}' {OutputDir}/snakemake/result/8_unannotated_protein_coding/lncScan/result_XGB.txt > {output[lncScan_XGB]} && \
		
		{params.workDir}/CPAT/bin/cpat.py -g {output[fa]} --antisense -d {CPATlogitModel} -x {CPATHexamer} -o {OutputDir}/snakemake/result/8_unannotated_protein_coding/CPAT/CPAToutput && \
		Rscript {params.workDir}/script/cpat_coding.R -i {OutputDir}/snakemake/result/8_unannotated_protein_coding/CPAT/CPAToutput.ORF_prob.tsv -p {CPATcutoff} -o {output[CPAT]}							
		"""
rule mRNA_overlap:
	input:
		lncScan_svm="{OutputDir}/snakemake/result/8_unannotated_protein_coding/lncScan/lncScan_svm.txt",
		CPAT="{OutputDir}/snakemake/result/8_unannotated_protein_coding/CPAT/CPAT.txt",
		lncScan_XGB="{OutputDir}/snakemake/result/8_unannotated_protein_coding/lncScan/lncScan_XGB.txt",
		gtf="{OutputDir}/snakemake/result/8_unannotated_protein_coding/filter_2.gtf"
	output:
		veen="{OutputDir}/snakemake/result/8_unannotated_protein_coding/overlap/veen-pdf.pdf",
		overlap="{OutputDir}/snakemake/result/8_unannotated_protein_coding/overlap/overlap.txt",
		unannotated_protein="{OutputDir}/snakemake/result/8_unannotated_protein_coding/unannotated_protein.gtf"
	params:
		workDir=config["appDir"]
	shell:
		"""
		python {params.workDir}/script/overlap.py -a {input[lncScan_svm]} -b {input[CPAT]} -c {input[lncScan_XGB]} -d {output[veen]} -e {output[overlap]} && \
		python {params.workDir}/script/extract_gtf.py -i {output[overlap]} -f {input[gtf]} -o {output[unannotated_protein]} && \
		gffcompare -r {GENCODE_protein_gtf} -o {OutputDir}/snakemake/result/8_unannotated_protein_coding/unMRNA {output[unannotated_protein]}
		"""
rule merge_gtf:
	input:
		"{OutputDir}/snakemake/result/9_unannotated_lncRNA/unannotated_lncRNA.gtf",
		"{OutputDir}/snakemake/result/8_unannotated_protein_coding/unannotated_protein.gtf",
		"{OutputDir}/snakemake/result/7_annotated_lncRNA/annotated_lncRNA.gtf",
		"{OutputDir}/snakemake/result/6_protein_coding/protein_coding.gtf",
		"{OutputDir}/snakemake/result/12_pseudogene/pseudogene.gtf",
		"{OutputDir}/snakemake/result/13_other_rnas/other_rna.gtf"
	output:
		"{OutputDir}/snakemake/result/10_featurecount/all.gtf",
		"{OutputDir}/snakemake/result/10_featurecount/all.txt",
		"{OutputDir}/snakemake/result/10_featurecount/unannotated_transcript.gtf",
		"{OutputDir}/snakemake/result/10_featurecount/annotated_transcript.gtf"
	params:
		workDir = config["appDir"]
	shell:
		"""
		python {params.workDir}/script/merge_gtf.py -a {input[3]} -b {input[2]} -c {input[0]} -d {input[1]} -e {ERCC} -f {input[4]} -g {input[5]} -o {output[1]} && \
		cat {input[0]} {input[1]} {input[2]} {input[3]} {input[4]} {input[5]} {ERCC} > {output[0]} && \
		cat {input[0]} {input[1]} > {output[2]} && \
		cat {input[2]} {input[3]} {ERCC} {input[4]} {input[5]} > {output[3]}
		"""
rule featurecounts_unannotated:
	input:
		gtf="{OutputDir}/snakemake/result/10_featurecount/unannotated_transcript.gtf"
	output:
		"{OutputDir}/snakemake/result/10_featurecount/unannotated_count.txt",
		"{OutputDir}/snakemake/result/10_featurecount/unannotated_count.txt.summary"
	conda:
		"env/InfoScan.yaml"
	shell:
		"featureCounts -Q 10 -t transcript -g transcript_id -f -p -T 16 -a {input[gtf]} -o {output[0]} {OutputDir}/snakemake/result/3_HISAT2_aligned/*.bam"
rule featurecounts_annotated:
	input:
		gtf="{OutputDir}/snakemake/result/10_featurecount/annotated_transcript.gtf"
		
	output:
		"{OutputDir}/snakemake/result/10_featurecount/annotated_count.txt",
		"{OutputDir}/snakemake/result/10_featurecount/annotated_count.txt.summary"
	conda:
		"env/InfoScan.yaml"
	message:
		"featurecounts to quantity"
	shell:
		"featureCounts -Q 10 -g gene_name -p -T 16 -a {input[gtf]} -o {output[0]} {OutputDir}/snakemake/result/3_HISAT2_aligned/*.bam" 
rule filter_fpkm:
	input:
		"{OutputDir}/snakemake/result/10_featurecount/annotated_count.txt",
		"{OutputDir}/snakemake/result/10_featurecount/unannotated_count.txt",
		"{OutputDir}/snakemake/result/10_featurecount/all.txt"
	output:
		"{OutputDir}/snakemake/result/10_featurecount/fpkm_matrix.txt"
	conda:
		"env/R_Library.yaml"
	message:
		"filter expression matrix by fpkm "
	params:
		workDir=config["appDir"]
	shell:
		"Rscript {params.workDir}/script/filter_by_fpkm.R -i {input[0]} -I {input[1]} -q {fpkm_threshold} -b {input[2]} -o {output[0]}"
rule filter_exon_number:
	input:
		"{OutputDir}/snakemake/result/10_featurecount/fpkm_matrix.txt",
		"{OutputDir}/snakemake/result/10_featurecount/all.gtf"
	output:
		"{OutputDir}/snakemake/result/11_data_analysis/exon_number.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/protein_coding_finally.gtf",
		"{OutputDir}/snakemake/result/11_data_analysis/annotated_lncRNA_finally.gtf",
		"{OutputDir}/snakemake/result/11_data_analysis/unannotated_lncRNA_finally.gtf",
		"{OutputDir}/snakemake/result/11_data_analysis/unannotated_protein_finally.gtf",
		"{OutputDir}/snakemake/result/11_data_analysis/finally_fpkm.txt",
		exon_number_png = "{OutputDir}/snakemake/result/11_data_analysis/exon_number.png"
	conda:
		"env/R_Library.yaml"
	params:
 		workDir=config["appDir"]
	shell:
		"Rscript {params.workDir}/script/exon_number.R -a {input[0]} -c {input[1]} -o {output[0]} -d {output[1]} -e {output[2]} -f {output[3]} -g {output[4]} -q {exon_filter}  -O {output[5]} -p {output[exon_number_png]}"
rule length_distribution:
	input:
		"{OutputDir}/snakemake/result/11_data_analysis/finally_fpkm.txt",
		"{OutputDir}/snakemake/result/10_featurecount/all.txt"
	output:
		"{OutputDir}/snakemake/result/11_data_analysis/gene-length-distrubution.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/before_fpkm_filter.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/after_fpkm_filter.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/gene-length-distrubution.png",
		"{OutputDir}/snakemake/result/11_data_analysis/before_fpkm_filter.png",
		"{OutputDir}/snakemake/result/11_data_analysis/after_fpkm_filter.png"
	conda:
		"env/R_Library.yaml"
	params:
		workDir=config["appDir"]
	shell:
		"Rscript {params.workDir}/script/length_distribution.R -a {input[0]} -b {input[1]} -o {output[0]} -c {output[1]} -d {output[2]} -e {output[3]} -f {output[4]} -g {output[5]}"
rule gene_expression:
	input:
		"{OutputDir}/snakemake/result/11_data_analysis/finally_fpkm.txt",
		"{OutputDir}/snakemake/result/10_featurecount/all.txt"
	output:
		"{OutputDir}/snakemake/result/10_featurecount/protein_coding.csv",
		"{OutputDir}/snakemake/result/10_featurecount/annotated_lncRNA.csv",
		"{OutputDir}/snakemake/result/10_featurecount/unannotated_lncRNA.csv",
		"{OutputDir}/snakemake/result/10_featurecount/unannotated_protein.csv",
		"{OutputDir}/snakemake/result/11_data_analysis/gene_FPKM_distrubution.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/gene_boxplot_FPKM.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/gene_boxplot_FPKM.png"	
	conda:
		"env/R_Library.yaml"
	params:
		workDir=config["appDir"]
	shell:
		"Rscript {params.workDir}/script/gene_expression_plot.R -a {input[0]} -b {input[1]} -c {output[0]} -d {output[1]} -e {output[2]} -f {output[3]} -g {output[4]} -j {output[5]} -p {output[6]}"	
rule extract_fa:
	input:
		"{OutputDir}/snakemake/result/11_data_analysis/unannotated_lncRNA_finally.gtf",
		"{OutputDir}/snakemake/result/11_data_analysis/unannotated_protein_finally.gtf"
	output:
		"{OutputDir}/snakemake/result/11_data_analysis/unannotated_lncRNA_finally.fa",
		"{OutputDir}/snakemake/result/11_data_analysis/unannotated_protein_finally.fa"
	conda:
		"env/InfoScan.yaml"
	shell:
		"""
		gffread -w {output[0]} -g {FA} {input[0]} && \
		gffread -w {output[1]} -g {FA} {input[1]}
		"""
rule conserved_score:
	input:
		"{OutputDir}/snakemake/result/11_data_analysis/unannotated_lncRNA_finally.gtf",
		"{OutputDir}/snakemake/result/11_data_analysis/unannotated_protein_finally.gtf",
		"{OutputDir}/snakemake/result/11_data_analysis/annotated_lncRNA_finally.gtf",
		"{OutputDir}/snakemake/result/11_data_analysis/protein_coding_finally.gtf",
		"{OutputDir}/snakemake/result/10_featurecount/all.txt"
	output:
		"{OutputDir}/snakemake/result/11_data_analysis/conservation/all_transcript.bed",
		"{OutputDir}/snakemake/result/11_data_analysis/conservation/ConservedScore.txt",
		"{OutputDir}/snakemake/result/11_data_analysis/conservation/gene_conservation_density.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/conservation/gene_conservation_ecdf.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/conservation/gene_conservation_ecdf.png"
	conda:
		"env/InfoScan.yaml"
	message:
		"calculate gene conservation"
	params:
		workDir=config["appDir"]
	shell:
		"""
		cat {input[0]} {input[1]} {input[2]} {input[3]} > {OutputDir}/snakemake/result/11_data_analysis/conservation/all_gene.gtf && \
        python {params.workDir}/script/GtfToBed.py {OutputDir}/snakemake/result/11_data_analysis/conservation/all_gene.gtf {output[0]} && \
        {params.workDir}/script/bigWigAverageOverBed {phastcons} {output[0]}  {output[1]} && \
		Rscript {params.workDir}/script/conservation_1.R -a {output[1]} -b {input[4]} -c {output[2]} -d {output[3]} -e {OutputDir}/snakemake/result/11_data_analysis/conservation/all_gene.gtf -p {output[4]}
		"""
rule transcripts_seq_features:
	input:
		"{OutputDir}/snakemake/result/11_data_analysis/unannotated_lncRNA_finally.fa",
		"{OutputDir}/snakemake/result/11_data_analysis/unannotated_protein_finally.fa",
		lncRNA = "{OutputDir}/snakemake/result/11_data_analysis/annotated_lncRNA_finally.gtf",
		mRNA = "{OutputDir}/snakemake/result/11_data_analysis/protein_coding_finally.gtf"
	output:
		"{OutputDir}/snakemake/result/11_data_analysis/Sequence_features/protein.fa",
		"{OutputDir}/snakemake/result/11_data_analysis/Sequence_features/lncRNA.fa",
		"{OutputDir}/snakemake/result/11_data_analysis/Sequence_features/features.txt",
		"{OutputDir}/snakemake/result/11_data_analysis/Sequence_features/Fickett_score.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/Sequence_features/Max_ORF_length.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/Sequence_features/Max_ORF_coverage.pdf",
		"{OutputDir}/snakemake/result/11_data_analysis/Sequence_features/Fickett_score.png",
		"{OutputDir}/snakemake/result/11_data_analysis/Sequence_features/Max_ORF_length.png",
		"{OutputDir}/snakemake/result/11_data_analysis/Sequence_features/Max_ORF_coverage.png"
	params:
		workDir=config["appDir"]
	shell:
		"""
		gffread -w {output[0]} -g {FA} {input[mRNA]} && \
		gffread -w {output[1]} -g {FA} {input[lncRNA]} && \
		python {params.workDir}/script/plot_extract_feature.py -a {output[0]} -b {output[1]} -c {input[0]} -d {input[1]} -e {output[3]} -f {output[4]} -g {output[5]} -o {output[2]} -k {output[6]} -i {output[7]} -j {output[8]}
		"""
rule classCode:
	input:
		"{OutputDir}/snakemake/result/11_data_analysis/unannotated_lncRNA_finally.gtf",
		"{OutputDir}/snakemake/result/11_data_analysis/unannotated_protein_finally.gtf"
	output:
		"{OutputDir}/snakemake/result/11_data_analysis/ClassCode/"
	shell:
		"""
		cat {input[0]} {input[1]} > {OutputDir}/snakemake/result/11_data_analysis/unannotated_transcripts.gtf && \
		gffcompare -r {GENCODE_protein_gtf} -o {OutputDir}/snakemake/result/11_data_analysis/ClassCode/cufcomp {OutputDir}/snakemake/result/11_data_analysis/unannotated_transcripts.gtf && \
		"""
rule generate_html_report:
	input:
		before_fpkm="{OutputDir}/snakemake/result/11_data_analysis/before_fpkm_filter.pdf",
		after_fpkm="{OutputDir}/snakemake/result/11_data_analysis/after_fpkm_filter.pdf",
		boxplot="{OutputDir}/snakemake/result/11_data_analysis/gene_boxplot_FPKM.pdf",
		exon="{OutputDir}/snakemake/result/11_data_analysis/exon_number.pdf",
		length="{OutputDir}/snakemake/result/11_data_analysis/gene-length-distrubution.pdf",
		ecdf="{OutputDir}/snakemake/result/11_data_analysis/conservation/gene_conservation_ecdf.pdf",
		fickett="{OutputDir}/snakemake/result/11_data_analysis/Sequence_features/Fickett_score.pdf",
		orf_len="{OutputDir}/snakemake/result/11_data_analysis/Sequence_features/Max_ORF_length.pdf",
		orf_cov="{OutputDir}/snakemake/result/11_data_analysis/Sequence_features/Max_ORF_coverage.pdf"
	output:
		html="{OutputDir}/snakemake/result/NovelScan_report_pdf.html"
	params:
		workDir=config["appDir"]
	script:
		"{params.workDir}/script/generate_report.py"

rule generate_html_report_png:
	input:
		before_fpkm="{OutputDir}/snakemake/result/11_data_analysis/before_fpkm_filter.png",
		after_fpkm="{OutputDir}/snakemake/result/11_data_analysis/after_fpkm_filter.png",
		boxplot="{OutputDir}/snakemake/result/11_data_analysis/gene_boxplot_FPKM.png",
		exon="{OutputDir}/snakemake/result/11_data_analysis/exon_number.png",
		length="{OutputDir}/snakemake/result/11_data_analysis/gene-length-distrubution.png",
		ecdf="{OutputDir}/snakemake/result/11_data_analysis/conservation/gene_conservation_ecdf.png",
		fickett="{OutputDir}/snakemake/result/11_data_analysis/Sequence_features/Fickett_score.png",
		orf_len="{OutputDir}/snakemake/result/11_data_analysis/Sequence_features/Max_ORF_length.png",
		orf_cov="{OutputDir}/snakemake/result/11_data_analysis/Sequence_features/Max_ORF_coverage.png"
	output:
		html="{OutputDir}/snakemake/result/NovelScan_report.html"
	params:
		workDir=config["appDir"]
	script:
		"{params.workDir}/script/generate_report_png.py"

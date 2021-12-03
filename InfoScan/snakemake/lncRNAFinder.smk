configfile: "snakemake/config/lncRNAFinder_config.yml"
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
pfam_1 = config["pfam"]
len_filter = config["length_filter"] 
fpkm_threshold = config["fpkm_filter"]
exon_filter=config["exon_filter"]
(SAMPLES,READS,) = glob_wildcards("snakemake/result/2.fastp_output/{sample}_clean_{read}.fq.gz")
READS=["1","2"]
rule all:
	input:
		expand("snakemake/result/6_protein_coding/protein_coding.gtf",sample=SAMPLES),
		expand("snakemake/result/7_annotated_lncRNA/annotated_lncRNA.gtf",sample=SAMPLES),
		expand("snakemake/result/8_unannotated_protein_coding/overlap/overlap.txt",sample=SAMPLES),
		expand("snakemake/result/8_unannotated_protein_coding/unannotated_protein.gtf",sample=SAMPLES),
		expand("snakemake/result/8_unannotated_protein_coding/unannotated_protein.fa",sample=SAMPLES),
		expand("snakemake/result/9_unannotated_lncRNA/overlap/overlap.txt",sample=SAMPLES),
		expand("snakemake/result/9_unannotated_lncRNA/unannotated_lncRNA.gtf",sample=SAMPLES),
		expand("snakemake/result/9_unannotated_lncRNA/unannotated_lncRNA.fa",sample=SAMPLES),
		expand("snakemake/result/10_featurecount/unannotated_transcript.gtf"),
		expand("snakemake/result/10_featurecount/fpkm_matrix.txt"),
		expand("snakemake/result/11_data_analysis/finally_fpkm.txt"),
		expand("snakemake/result/11_data_analysis/gene_FPKM_distrubution.pdf"),
		expand("snakemake/result/11_data_analysis/gene-length-distrubution.pdf"),
		expand("snakemake/result/11_data_analysis/exon_number.pdf"),
		expand("snakemake/result/11_data_analysis/unannotated_lncRNA_finally.fa"),
		expand("snakemake/result/11_data_analysis/sclncRNAfinder.pdf")

rule protein_coding:
	input:
		"snakemake/result/5_stringtie_merge/stringtie_merged.gtf",
		"snakemake/result/5_stringtie_merge/stringtie_merged_transcript.gtf"
	output:
		"snakemake/result/6_protein_coding/protein_coding_1.gtf",
		"snakemake/result/6_protein_coding/protein_name_1.txt",
		"snakemake/result/6_protein_coding/protein_coding.gtf"
	message:
		"identity annotated protein coding gene"
	shell:
		"""
		bedtools intersect -a {GENCODE_protein_gtf} -b {input[1]} -f 0.9 > {output[0]} && \
		snakemake/script/extract_transcript_name.sh {output[0]} {output[1]} && \
		Rscript snakemake/script/extract_gtf.R -i {input[0]} -I {output[1]} -o {output[2]}
		"""
rule annotated_lncRNA:
	input:
		"snakemake/result/5_stringtie_merge/stringtie_merged.gtf",
		"snakemake/result/5_stringtie_merge/stringtie_merged_transcript.gtf"
	output:
		"snakemake/result/7_annotated_lncRNA/annotated_lncRNA_1.gtf",
		"snakemake/result/7_annotated_lncRNA/annotated_lncRNA_name_1.txt",
		"snakemake/result/7_annotated_lncRNA/annotated_lncRNA.gtf"
	message:
		"identity annotated lncRNA"
	shell:
		"""
		bedtools intersect -a {GENCODE_lncRNA_gene_gtf} -b {input[1]} -f 0.9 > {output[0]} && \
		snakemake/script/extract_transcript_name.sh {output[0]} {output[1]} && \
		Rscript snakemake/script/extract_gtf.R -i {input[0]} -I {output[1]} -o {output[2]}
		"""
rule unannotated_lncRNA:
	input:
		"snakemake/result/5_stringtie_merge/stringtie_merged.gtf",
		"snakemake/result/5_stringtie_merge/stringtie_merged_transcript.gtf"
	output:
		"snakemake/result/9_unannotated_lncRNA/cufcomp.stringtie_merged.gtf.tmap",
		"snakemake/result/9_unannotated_lncRNA/lncRNA_1.gtf",
#		"snakemake/result/9_unannotated_lncRNA/lncRNA_1.fa",
		rnasamba="snakemake/result/9_unannotated_lncRNA/rnasamba/rnasamba.txt",
		CPAT="snakemake/result/9_unannotated_lncRNA/CPAT/CPAT.txt",
		lncfinder="snakemake/result/9_unannotated_lncRNA/lncfinder/lncfinder.txt",
		veen="snakemake/result/9_unannotated_lncRNA/overlap/veen-pdf.pdf",
		overlap="snakemake/result/9_unannotated_lncRNA/overlap/overlap.txt",
		lncRNA="snakemake/result/9_unannotated_lncRNA/lncRNA_2.gtf"
	message:
		"identity unannotated lncRNA"
	shell:
		"""
		cuffcompare -r {GTF} -o snakemake/result/9_unannotated_lncRNA/cufcomp {input[0]} && \
		mv snakemake/result/5_stringtie_merge/cufcomp.stringtie_merged.gtf.tmap  {output[0]} && \
		snakemake/script/length_filter.sh {len_filter} {output[0]} snakemake/result/9_unannotated_lncRNA/filter.length_1 snakemake/result/9_unannotated_lncRNA/filter.length_1 snakemake/result/9_unannotated_lncRNA/filter.class_2 snakemake/result/9_unannotated_lncRNA/filter.class_2 snakemake/result/9_unannotated_lncRNA/filter.txt && \
		Rscript snakemake/script/extract_gtf_2.R -i snakemake/result/9_unannotated_lncRNA/cufcomp.combined.gtf -I snakemake/result/9_unannotated_lncRNA/filter.txt -o {output[1]} && \
		gffread  {output[1]} -o- > snakemake/result/9_unannotated_lncRNA/lncRNA_1.gff3 && \
		gffread snakemake/result/9_unannotated_lncRNA/lncRNA_1.gff3 --force-exons -T -o {output[lncRNA]} && \
		gffread -w snakemake/result/9_unannotated_lncRNA/lncRNA_1.fa -g {FA} {output[lncRNA]} && \
		
		rnasamba classify snakemake/result/9_unannotated_lncRNA/rnasamba/result.txt  snakemake/result/9_unannotated_lncRNA/lncRNA_1.fa {model} && \
		cat snakemake/result/9_unannotated_lncRNA/rnasamba/result.txt | grep 'noncoding' | cut -f 1 > snakemake/result/9_unannotated_lncRNA/rnasamba/noncoding-transcript.txt && \
		mv snakemake/result/9_unannotated_lncRNA/rnasamba/noncoding-transcript.txt {output[rnasamba]} && \
		
		Rscript snakemake/script/lncFinder.R -i snakemake/result/9_unannotated_lncRNA/lncRNA_1.fa -o snakemake/result/9_unannotated_lncRNA/lncfinder/result.txt -s {specie} && \
		cat snakemake/result/9_unannotated_lncRNA/lncfinder/result.txt | grep 'NonCoding'| cut -f 1 > snakemake/result/9_unannotated_lncRNA/lncfinder/noncoding-transcript.txt && \
		mv snakemake/result/9_unannotated_lncRNA/lncfinder/noncoding-transcript.txt {output[lncfinder]}

		cpat.py -g  snakemake/result/9_unannotated_lncRNA/lncRNA_1.fa --antisense -d CPAT-3.0.0/dat/{CPATlogitModel} -x CPAT-3.0.0/dat/{CPATHexamer} -o snakemake/result/9_unannotated_lncRNA/CPAT/CPAToutput && \
		Rscript snakemake/script/cpat_noncoding.R -i snakemake/result/9_unannotated_lncRNA/CPAT/CPAToutput.ORF_prob.tsv -p {CPATcutoff} -I snakemake/result/9_unannotated_lncRNA/CPAT/CPAToutput.no_ORF.txt -o {output[CPAT]} && \
								
		Rscript snakemake/script/overlap.R -a {output[rnasamba]} -b {output[lncfinder]} -c {output[CPAT]} -v {output[veen]} -o {output[overlap]} && \
		sed -i 's/^[ ]*//g'  {output[overlap]}
		"""
rule unannotated_lncRNA_2:
	input:
		"snakemake/result/9_unannotated_lncRNA/lncRNA_2.gtf",
		"snakemake/result/9_unannotated_lncRNA/overlap/overlap.txt"
	output:	
		"snakemake/result/9_unannotated_lncRNA/unannotated_lncRNA_3.gtf"
	shell:
		"Rscript snakemake/script/extract_gtf.R -i {input[0]} -I {input[1]} -o {output[0]}"

rule unannotated_lncRNA_3:
	input:
		"snakemake/result/9_unannotated_lncRNA/unannotated_lncRNA_3.gtf"
	output:
		"snakemake/result/9_unannotated_lncRNA/unannotated_lncRNA_3.fa"
	shell:
		"gffread -w {output} -g {FA} {input}"
#使用Transeq 将转录本序列翻译为6个可能的蛋白序列
rule unannotated_lncRNA_4:
	input:
		"snakemake/result/9_unannotated_lncRNA/unannotated_lncRNA_3.fa",
		"snakemake/result/9_unannotated_lncRNA/overlap/overlap.txt",
		"snakemake/result/9_unannotated_lncRNA/unannotated_lncRNA_3.gtf"
	output:
		"snakemake/result/9_unannotated_lncRNA/unannotated_lncRNA.fa",
		"snakemake/result/9_unannotated_lncRNA/unannotated_lncRNA.gtf",
		"snakemake/result/9_unannotated_lncRNA/pcl.fa",
		"snakemake/result/9_unannotated_lncRNA/pcl.gtf",
		"snakemake/result/9_unannotated_lncRNA/pfam/protein.fa"
	shell:
		"""
		transeq {input[0]} {output[4]} -frame=6 &&\
		hmmscan  -o snakemake/result/9_unannotated_lncRNA/pfam/out.txt --tblout snakemake/result/9_unannotated_lncRNA/pfam/tblout.txt --domtblout snakemake/result/9_unannotated_lncRNA/domout.txt --noali -E 1e-5 {pfam_1} {output[4]} && \
		grep -v '^#' snakemake/result/9_unannotated_lncRNA/pfam/tblout.txt | awk '($5< 1e-5){{print $3}}'| cut -b 1-14 | sort | uniq > snakemake/result/9_unannotated_lncRNA/pfam/coding.id && \
		grep -v -f snakemake/result/9_unannotated_lncRNA/pfam/coding.id {input[1]} > snakemake/result/9_unannotated_lncRNA/pfam/filter_transcript_ID && \
		grep -Ff snakemake/result/9_unannotated_lncRNA/pfam/filter_transcript_ID {input[2]} > snakemake/result/9_unannotated_lncRNA/tmp.gtf && \
		gffread --force-exons snakemake/result/9_unannotated_lncRNA/tmp.gtf -T > {output[1]} && \
		gffread -w {output[0]} -g {FA} {output[1]} && \
		grep -Ff snakemake/result/9_unannotated_lncRNA/pfam/coding.id {input[2]} > {output[3]} && \
		gffread -w {output[2]} -g {FA} {output[3]}
		"""
rule unannotated_protein:
	input:
		"snakemake/result/5_stringtie_merge/stringtie_merged.gtf",
		"snakemake/result/5_stringtie_merge/stringtie_merged_transcript.gtf",
		"snakemake/result/9_unannotated_lncRNA/cufcomp.stringtie_merged.gtf.tmap"
	output:
		"snakemake/result/8_unannotated_protein_coding/cufcomp.stringtie_merged.gtf.tmap",
		"snakemake/result/8_unannotated_protein_coding/filter_1.gtf",
		"snakemake/result/8_unannotated_protein_coding/filter_2.gtf",
		"snakemake/result/8_unannotated_protein_coding/filter_2.fa",
		rnasamba="snakemake/result/8_unannotated_protein_coding/rnasamba/rnasamba.txt",
		lncfinder="snakemake/result/8_unannotated_protein_coding/lncfinder/lncfinder.txt",
		CPAT="snakemake/result/8_unannotated_protein_coding/CPAT/CPAT.txt",
		veen="snakemake/result/8_unannotated_protein_coding/overlap/veen-pdf.pdf",
		overlap="snakemake/result/8_unannotated_protein_coding/overlap/overlap.txt",
		unannotated_protein_1="snakemake/result/8_unannotated_protein_coding/unannotated_protein_1.gtf",
		unannotated_protein_1_fa="snakemake/result/8_unannotated_protein_coding/unannotated_protein_1.fa"
	message:
		"identity unannotated protein"
	shell:
		"""
		cp {input[2]}  {output[0]} && \
		snakemake/script/protein_filter.sh  {output[0]} snakemake/result/8_unannotated_protein_coding/filter1.txt snakemake/result/8_unannotated_protein_coding/filter1.txt snakemake/result/8_unannotated_protein_coding/filter.txt && \
		Rscript snakemake/script/extract_gtf_2.R -i snakemake/result/9_unannotated_lncRNA/cufcomp.combined.gtf -I snakemake/result/8_unannotated_protein_coding/filter.txt -o {output[1]} && \
		gffread  {output[1]} -o- > snakemake/result/8_unannotated_protein_coding/filter_1.gff3 && \
		gffread snakemake/result/8_unannotated_protein_coding/filter_1.gff3 --force-exons -T -o {output[2]} && \
		gffread -w {output[3]} -g {FA} {output[2]} && \

		rnasamba classify snakemake/result/8_unannotated_protein_coding/rnasamba/result.txt  {output[3]} {model} && \
		cat snakemake/result/8_unannotated_protein_coding/rnasamba/result.txt|cut -f 1 > snakemake/result/8_unannotated_protein_coding/rnasamba/rnasambaid.txt && \
		cat snakemake/result/8_unannotated_protein_coding/rnasamba/result.txt | grep 'noncoding' | cut -f 1 > snakemake/result/8_unannotated_protein_coding/rnasamba/noncoding-transcript.txt && \
		grep -v -f snakemake/result/8_unannotated_protein_coding/rnasamba/noncoding-transcript.txt snakemake/result/8_unannotated_protein_coding/rnasamba/rnasambaid.txt > snakemake/result/8_unannotated_protein_coding/rnasamba/coding-transcript.txt && \
		mv snakemake/result/8_unannotated_protein_coding/rnasamba/coding-transcript.txt {output[rnasamba]} && \
		sed -i '1d' {output[rnasamba]} && \

		Rscript snakemake/script/lncFinder.R -i {output[3]} -o snakemake/result/8_unannotated_protein_coding/lncfinder/result.txt -s {specie} && \
		cat snakemake/result/8_unannotated_protein_coding/lncfinder/result.txt |cut -f 1 > snakemake/result/8_unannotated_protein_coding/lncfinder/lncfinderid.txt && \
		cat snakemake/result/8_unannotated_protein_coding/lncfinder/result.txt | grep 'NonCoding'| cut -f 1 > snakemake/result/8_unannotated_protein_coding/lncfinder/noncoding-transcript.txt && \
		grep -v -f snakemake/result/8_unannotated_protein_coding/lncfinder/noncoding-transcript.txt snakemake/result/8_unannotated_protein_coding/lncfinder/lncfinderid.txt > snakemake/result/8_unannotated_protein_coding/lncfinder/coding-transcript.txt && \
		mv snakemake/result/8_unannotated_protein_coding/lncfinder/coding-transcript.txt {output[lncfinder]} && \
		
		cpat.py -g  {output[3]} --antisense -d CPAT-3.0.0/dat/{CPATlogitModel} -x CPAT-3.0.0/dat/{CPATHexamer} -o snakemake/result/8_unannotated_protein_coding/CPAT/CPAToutput && \
		Rscript snakemake/script/cpat_coding.R -i snakemake/result/8_unannotated_protein_coding/CPAT/CPAToutput.ORF_prob.tsv -p {CPATcutoff} -o {output[CPAT]}
								
		Rscript snakemake/script/overlap_coding.R -a {output[rnasamba]} -b {output[lncfinder]} -c {output[CPAT]} -v {output[veen]} -o {output[overlap]} && \
		sed -i 's/^[ ]*//g'  {output[overlap]} && \

		Rscript snakemake/script/extract_gtf.R -i {output[2]} -I {output[overlap]} -o {output[unannotated_protein_1]} && \
		gffread -w {output[unannotated_protein_1_fa]} -g {FA} {output[unannotated_protein_1]} 
		"""
rule unannotated_protein_1:
	input:
		"snakemake/result/8_unannotated_protein_coding/unannotated_protein_1.fa",
		"snakemake/result/8_unannotated_protein_coding/overlap/overlap.txt",
		"snakemake/result/8_unannotated_protein_coding/unannotated_protein_1.gtf"
	output:
		"snakemake/result/8_unannotated_protein_coding/unannotated_protein.fa",
		"snakemake/result/8_unannotated_protein_coding/unannotated_protein.gtf",
		"snakemake/result/8_unannotated_protein_coding/pfam/protein.fa"
	shell:
		"""
		transeq {input[0]} {output[2]} -frame=6 &&\
		hmmscan  -o snakemake/result/8_unannotated_protein_coding/pfam/out.txt --tblout snakemake/result/8_unannotated_protein_coding/pfam/tblout.txt --domtblout domout.txt --noali -E 1e-5 {pfam_1} {output[2]} && \
		grep -v '^#' snakemake/result/8_unannotated_protein_coding/pfam/tblout.txt | awk '($5< 1e-5){{print $3}}'| cut -b 1-14 | sort | uniq > snakemake/result/8_unannotated_protein_coding/pfam/coding.id && \
		grep -Ff snakemake/result/8_unannotated_protein_coding/pfam/coding.id {input[2]} > snakemake/result/8_unannotated_protein_coding/tmp.gtf && \
		gffread --force-exons snakemake/result/8_unannotated_protein_coding/tmp.gtf -T > snakemake/result/8_unannotated_protein_coding/unannotated_protein.gtf && \
		gffread -w {output[0]} -g {FA} {output[1]}
		"""
rule merge_gtf:
	input:
		"snakemake/result/9_unannotated_lncRNA/unannotated_lncRNA.gtf",
		"snakemake/result/8_unannotated_protein_coding/unannotated_protein.gtf",
		"snakemake/result/7_annotated_lncRNA/annotated_lncRNA.gtf",
		"snakemake/result/6_protein_coding/protein_coding.gtf",
		"snakemake/result/9_unannotated_lncRNA/pcl.gtf"
	output:
		"snakemake/result/10_featurecount/all.gtf",
		"snakemake/result/10_featurecount/all.txt",
		"snakemake/result/10_featurecount/unannotated_transcript.gtf",
		"snakemake/result/10_featurecount/annotated_transcript.gtf"
	shell:
		"""
		Rscript snakemake/script/merge_gtf.R -a {input[3]} -b {input[2]} -c {input[0]} -d {input[1]} -e {ERCC} -o {output[1]} && \
		cat {input[0]} {input[1]} {input[2]} {input[3]} {input[4]} {ERCC} > {output[0]} && \
		cat {input[0]} {input[1]} {input[4]} > {output[2]} && \
		cat {input[2]} {input[3]} {ERCC} > {output[3]}
		"""
rule featurecounts_unannotated:
	input:
		bam_file=expand("snakemake/result/3_HISAT2_aligned/{sample}.bam",sample=SAMPLES),
		gtf="snakemake/result/10_featurecount/unannotated_transcript.gtf"
	output:
		"snakemake/result/10_featurecount/unannotated_count.txt",
		"snakemake/result/10_featurecount/unannotated_count.txt.summary"
	shell:
		"featureCounts -Q 10 -t transcript -g transcript_id -f -p -T 16 -a {input[gtf]} -o {output[0]} snakemake/result/3_HISAT2_aligned/*.bam"
rule featurecounts_annotated:
	input:
		bam_file=expand("snakemake/result/3_HISAT2_aligned/{sample}.bam",sample=SAMPLES),
		gtf="snakemake/result/10_featurecount/annotated_transcript.gtf"
		
	output:
		"snakemake/result/10_featurecount/annotated_count.txt",
		"snakemake/result/10_featurecount/annotated_count.txt.summary"
	message:
		"featurecounts to quantity"
	shell:
		"featureCounts -Q 10 -g gene_name -p -T 16 -a {input[gtf]} -o {output[0]} snakemake/result/3_HISAT2_aligned/*.bam" 
rule filter_fpkm:
	input:
		"snakemake/result/10_featurecount/annotated_count.txt",
		"snakemake/result/10_featurecount/unannotated_count.txt",
		"snakemake/result/10_featurecount/all.txt"
	output:
		"snakemake/result/10_featurecount/fpkm_matrix.txt"
	message:
		"filter expression matrix by fpkm "
	shell:
		"Rscript snakemake/script/filter_by_fpkm.R -i {input[0]} -I {input[1]} -q {fpkm_threshold} -b {input[2]} -o {output[0]}"
rule filter_exon_number:
	input:
		"snakemake/result/10_featurecount/fpkm_matrix.txt",
		"snakemake/result/10_featurecount/all.gtf"
	output:
		"snakemake/result/11_data_analysis/exon_number.pdf",
		"snakemake/result/11_data_analysis/protein_coding_finally.gtf",
		"snakemake/result/11_data_analysis/annotated_lncRNA_finally.gtf",
		"snakemake/result/11_data_analysis/unannotated_lncRNA_finally.gtf",
		"snakemake/result/11_data_analysis/unannotated_protein_finally.gtf",
		"snakemake/result/11_data_analysis/finally_fpkm.txt"
	shell:
		"Rscript snakemake/script/exon_number.R -a {input[0]} -c {input[1]} -o {output[0]} -d {output[1]} -e {output[2]} -f {output[3]} -g {output[4]} -q {exon_filter}  -O {output[5]}"
rule length_distribution:
	input:
		"snakemake/result/11_data_analysis/finally_fpkm.txt",
		"snakemake/result/10_featurecount/all.txt"
	output:
		"snakemake/result/11_data_analysis/gene-length-distrubution.pdf",
		"snakemake/result/11_data_analysis/before_fpkm_filter.pdf",
		"snakemake/result/11_data_analysis/after_fpkm_filter.pdf"
	shell:
		"Rscript snakemake/script/length_distribution.R -a {input[0]} -b {input[1]} -o {output[0]} -c {output[1]} -d {output[2]}"
rule gene_expression:
	input:
		"snakemake/result/11_data_analysis/finally_fpkm.txt",
		"snakemake/result/10_featurecount/all.txt"
	output:
		"snakemake/result/10_featurecount/protein_coding.csv",
		"snakemake/result/10_featurecount/annotated_lncRNA.csv",
		"snakemake/result/10_featurecount/unannotated_lncRNA.csv",
		"snakemake/result/10_featurecount/unannotated_protein.csv",
		"snakemake/result/11_data_analysis/gene_FPKM_distrubution.pdf",
		"snakemake/result/11_data_analysis/gene_boxplot_FPKM.pdf"	
	shell:
		"Rscript snakemake/script/gene_expression_plot.R -a {input[0]} -b {input[1]} -c {output[0]} -d {output[1]} -e {output[2]} -f {output[3]} -g {output[4]} -j {output[5]}"
	
rule extract_fa:
	input:
		"snakemake/result/11_data_analysis/unannotated_lncRNA_finally.gtf",
		"snakemake/result/11_data_analysis/unannotated_protein_finally.gtf"
	output:
		"snakemake/result/11_data_analysis/unannotated_lncRNA_finally.fa",
		"snakemake/result/11_data_analysis/unannotated_protein_finally.fa"
	shell:
		"""
		gffread -w {output[0]} -g {FA} {input[0]} && \
		gffread -w {output[1]} -g {FA} {input[1]}
		"""
rule merge_pdf:
	input:
		"snakemake/result/8_unannotated_protein_coding/overlap/veen-pdf.pdf",
		"snakemake/result/9_unannotated_lncRNA/overlap/veen-pdf.pdf",
		"snakemake/result/11_data_analysis/gene_boxplot_FPKM.pdf",
		"snakemake/result/11_data_analysis/gene-length-distrubution.pdf",
		"snakemake/result/11_data_analysis/exon_number.pdf",
		"snakemake/result/11_data_analysis/before_fpkm_filter.pdf",
		"snakemake/result/11_data_analysis/after_fpkm_filter.pdf"
	output:
		"snakemake/result/11_data_analysis/sclncRNAfinder.pdf"
	shell:
		"""
		mutool merge -o {output} {input[0]} {input[1]} {input[2]} {input[3]} {input[4]} {input[5]} {input[6]} && \
		mupdf-x11 {output}
		"""

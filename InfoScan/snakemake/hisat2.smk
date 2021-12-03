configfile: "snakemake/config/hisat2_config.yml"
GTF = config["GTFfile"]

(SAMPLES,READS,) = glob_wildcards("snakemake/result/2.fastp_output/{sample}_clean_{read}.fq.gz")
READS=["1","2"]

rule all:
	input:
		expand("snakemake/result/3_HISAT2_aligned/{sample}.sam",sample=SAMPLES),
		expand("snakemake/result/3_HISAT2_aligned/{sample}.bam",sample=SAMPLES),
		expand("snakemake/result/3_HISAT2_aligned/{sample}.bam.bai",sample=SAMPLES),
		expand("snakemake/result/4_stringtie/{sample}.gtf",sample=SAMPLES),
		expand("snakemake/result/5_stringtie_merge/stringtie_merged.gtf",sample=SAMPLES)
rule HISAT2:
	input:
        	"snakemake/result/2.fastp_output/{sample}_clean_1.fq.gz",
        	"snakemake/result/2.fastp_output/{sample}_clean_2.fq.gz"   	
	output:
			"snakemake/result/3_HISAT2_aligned/{sample}.sam"
	message:
		"Begin aligned"
	params: 
		thread = config["thread"],
		index = config["GenomeIndex"],
		params1 = config["params_fast"]+config["params_sensitive"]+config["params_verysensitive"]+" --min-intronlen "+config["min_intron"]+" --max-intronlen "+config["max_intron"]+" "+config["params_qc"]
	log:
		"snakemake/logs/HISAT2/{sample}_HISAT2.log"
	script:
		"script/hisat2.py"
rule samtools:
	input:
		"snakemake/result/3_HISAT2_aligned/{sample}.sam"
	output:
		"snakemake/result/3_HISAT2_aligned/{sample}.bam",
		"snakemake/result/3_HISAT2_aligned/{sample}.bam.bai"
	shell:
		"samtools view -b {input} | samtools sort -O bam -@ 12 -o {output[0]} && samtools index -@ 12 -b {output[0]} {output[1]}" 

rule stringtie:
	input:
		"snakemake/result/3_HISAT2_aligned/{sample}.bam"   	
	output:
		"snakemake/result/4_stringtie/{sample}.gtf",
		"snakemake/result/4_stringtie/{sample}.tab"
	message:
		"Begin stringtie transcript assembly"
	shell:
		"stringtie -B -p 12 -G {GTF} -o {output[0]} -A {output[1]} -i {input}"
rule merge:
	input:
		expand("snakemake/result/3_HISAT2_aligned/{sample}.bam",sample=SAMPLES),
		expand("snakemake/result/4_stringtie/{sample}.gtf",sample=SAMPLES)
	output:
		"snakemake/result/5_stringtie_merge/gtf.txt",
		"snakemake/result/5_stringtie_merge/stringtie_merged.gtf",
		"snakemake/result/5_stringtie_merge/stringtie_merged_transcript.gtf"
	shell:
		"""
		realpath  snakemake/result/4_stringtie/*.gtf >> {output[0]} && \
		stringtie --merge -p 16 -G {GTF} -o {output[1]} {output[0]} && \
		grep -w transcript {output[1]} > {output[2]}
		"""

configfile: "snakemake/config/fastp_config.yml"
data_dir=config["data_dir"]

(SAMPLES,READS,) = glob_wildcards("snakemake/data/{sample}_{read}.fastq.gz")
READS=["1","2"]

rule all:
	input:
		expand("snakemake/result/2.fastp_output/{sample}_clean_1.fq.gz",sample=SAMPLES),
		expand("snakemake/result/2.fastp_output/{sample}_clean_2.fq.gz",sample=SAMPLES)

rule fastp: 
	input:
		"snakemake/data/{sample}_1.fastq.gz",
		"snakemake/data/{sample}_2.fastq.gz"
	output: 
		"snakemake/result/2.fastp_output/{sample}_clean_1.fq.gz",
		"snakemake/result/2.fastp_output/{sample}_clean_2.fq.gz",
		html="snakemake/result/2.fastp_output/{sample}.html",
		json="snakemake/result/2.fastp_output/{sample}.json"
	message:
		"Begin Quality Control"  
	params:
		adapters='--adapter_sequence ' + config["adapter_read1"] + ' --adapter_sequence_r2 ' + config["adapter_read2"],
		quality_phred = config["quality_phred"],
		thread = config["thread"],
		adapter_read1 = config["adapter_read1"],
		adapter_read2 = config["adapter_read2"]
	script:
		"script/fastp.py"
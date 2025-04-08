#configfile: "snakemake/config/fastp_config.yml"
data_dir=config["data_dir"]
OutputDir=config["outputdir"]
workDir=config["appDir"]
sample=config["sample"]
#(WORKDIR,SAMPLES,) = glob_wildcards("{workDir}/data/{sample}_1.fastq.gz")
READS=["1","2"]

rule all:
	input:
#		expand("{OutputDir}/snakemake/result/2.fastp_output/",OutputDir=OutputDir),
		expand("{OutputDir}/snakemake/result/2.fastp_output/{sample}_clean_2.fq.gz",sample=sample,OutputDir=OutputDir),
		expand("{OutputDir}/snakemake/sample_config.yaml",OutputDir=OutputDir)
#		expand("snakemake/sample.yaml")
rule sample_name:
	input:
		expand("{OutputDir}/data/{sample}_1.fastq.gz",sample=sample,OutputDir=OutputDir)
	output:
		"{OutputDir}/snakemake/sample_config.yaml"
	params:
		Outputdir=config["outputdir"]
	run:
		with open(output[0], 'w') as out:
			sample_list= []
			for i in input:
				sample1 = i.split('/')[-1]
				sample2 = sample1.split('_')[0]
				sample_list.append(sample2)
			out.write('"'+'sample'+'"'+':'+' '+str(sample_list))
			out.write('\n')
			out.write('"'+'outputdir'+'"'+':'+' '+params.Outputdir)
rule sample_name2:
	input:
		expand("snakemake/data/{sample}_1.fastq.gz",sample=sample)
	output:
		"snakemake/sample.yaml"
	params:
		Outputdir=config["outputdir"]
	run:
		with open(output[0], 'w') as out:
			sample_list= []
			for i in input:
				sample1 = i.split('/')[-1]
				sample2 = sample1.split('_')[0]
				sample_list.append(sample2)
			out.write('"'+'sample'+'"'+':'+' '+str(sample_list))
			out.write('\n')
			out.write('"'+'outputdir'+'"'+':'+' '+params.Outputdir)
rule fastp: 
	input:
		"{OutputDir}/data/{sample}_1.fastq.gz",
		"{OutputDir}/data/{sample}_2.fastq.gz"
	output: 
		"{OutputDir}/snakemake/result/2.fastp_output/{sample}_clean_1.fq.gz",
		"{OutputDir}/snakemake/result/2.fastp_output/{sample}_clean_2.fq.gz",
		html="{OutputDir}/snakemake/result/2.fastp_output/{sample}.html"
	message:
		"Begin Quality Control"
#	conda:
#		"env/InfoScan.yaml"
	params:
		adapters='--adapter_sequence ' + config["adapter_read1"] + ' --adapter_sequence_r2 ' + config["adapter_read2"],
		quality_phred = config["quality_phred"],
		thread = config["thread"],
		adapter_read1 = config["adapter_read1"],
		adapter_read2 = config["adapter_read2"],
		workDir=config["appDir"]  # 新增此行
	script:
		"{params.workDir}/script/fastp.py"

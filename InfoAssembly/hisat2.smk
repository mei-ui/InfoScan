OutputDir=config["outputDir"]
SAMPLES=config["sample"]
GTF = config["GTFfile"]
workDir=config["appDir"]

#(SAMPLES,READS,) = glob_wildcards("{OutputDir}/snakemake/result/2.fastp_output/{sample}_clean_{read}.fq.gz")
READS=["1","2"]

rule all:
	input:
		expand("{OutputDir}/snakemake/result/3_HISAT2_aligned/{sample}.bam",sample=SAMPLES,OutputDir=OutputDir),
		expand("{OutputDir}/snakemake/result/3_HISAT2_aligned/{sample}.bam.bai",sample=SAMPLES,OutputDir=OutputDir),
		expand("{OutputDir}/snakemake/result/4_stringtie/{sample}.gtf",sample=SAMPLES,OutputDir=OutputDir),
		expand("{OutputDir}/snakemake/result/5_stringtie_merge/stringtie_merged.gtf",sample=SAMPLES,OutputDir=OutputDir)

rule HISAT2:
	input:
		"{OutputDir}/snakemake/result/2.fastp_output/{sample}_clean_1.fq.gz",
		"{OutputDir}/snakemake/result/2.fastp_output/{sample}_clean_2.fq.gz"
	output:
		temp("{OutputDir}/snakemake/result/3_HISAT2_aligned/{sample}.sam")
	message:
		"Begin aligned"
	conda:
		"env/samtool.yaml"
	params: 
		thread = config["thread"],
		index = config["GenomeIndex"],
		workDir=config["appDir"],
		params1 = config["params_fast"]+config["params_sensitive"]+config["params_verysensitive"]+" --min-intronlen "+config["min_intron"]+" --max-intronlen "+config["max_intron"]+" "+config["params_qc"]
	log:
		"{OutputDir}/snakemake/logs/HISAT2/{sample}_HISAT2.log"
	script:
		"{params.workDir}/script/hisat2.py"
rule samtools:
	input:
		"{OutputDir}/snakemake/result/3_HISAT2_aligned/{sample}.sam"
	output:
		"{OutputDir}/snakemake/result/3_HISAT2_aligned/{sample}.bam",
		"{OutputDir}/snakemake/result/3_HISAT2_aligned/{sample}.bam.bai"
	conda:
		"env/samtool.yaml"
	shell:
		"samtools view -b {input} | samtools sort -O bam -@ 12 -o {output[0]} && samtools index -@ 12 -b {output[0]} {output[1]}" 

rule stringtie:
	input:
		"{OutputDir}/snakemake/result/3_HISAT2_aligned/{sample}.bam"   	
	output:
		"{OutputDir}/snakemake/result/4_stringtie/{sample}.gtf"
	conda:
		"env/samtool.yaml"	
	message:
		"Begin stringtie transcript assembly"
	shell:
		"stringtie -B -p 12 -G {GTF} -o {output} -i {input}"
rule merge:
    input:
        expand("{OutputDir}/snakemake/result/3_HISAT2_aligned/{sample}.bam", sample=SAMPLES, OutputDir=OutputDir),
        expand("{OutputDir}/snakemake/result/4_stringtie/{sample}.gtf", sample=SAMPLES, OutputDir=OutputDir)
    output:
        "{OutputDir}/snakemake/result/5_stringtie_merge/gtf.txt",
        "{OutputDir}/snakemake/result/5_stringtie_merge/stringtie_merged.gtf",
        "{OutputDir}/snakemake/result/5_stringtie_merge/stringtie_merged_transcript.gtf"
    conda:
        "env/samtool.yaml"
    params:
        workDir=config["appDir"],
        OutputDir=lambda wildcards: wildcards.OutputDir
    shell:
        """
        ulimit -n 4096 && \
        {params.workDir}/script/realpath/realpath {params.OutputDir}/snakemake/result/4_stringtie/*.gtf > {output[0]} && \
        
        # 判断行数是否大于4000
        line_count=$(wc -l < {output[0]})
        if [ "$line_count" -gt 3000 ]; then
            split -l 3000 {output[0]} {params.OutputDir}/snakemake/result/4_stringtie/split_
            FILES_TO_MERGE={params.OutputDir}/snakemake/result/4_stringtie/split_*
        else
            # 直接使用原文件，无需拆分
            FILES_TO_MERGE={output[0]}
        fi

        # 处理所有需要合并的文件
        MERGED_FILES=""
        for file in $FILES_TO_MERGE; do
            output_merged="merged_$(basename "$file").gtf"
            stringtie --merge -p 16 -G {GTF} -o "$output_merged" "$file"
            MERGED_FILES="$MERGED_FILES $output_merged"
        done

        # 合并最终结果
        echo $MERGED_FILES | xargs {params.workDir}/script/realpath/realpath > {params.OutputDir}/snakemake/result/5_stringtie_merge/gtf_2.txt && \
        stringtie --merge -p 16 -G {GTF} -o {output[1]} $(cat {params.OutputDir}/snakemake/result/5_stringtie_merge/gtf_2.txt) && \
        grep -w transcript {output[1]} > {output[2]}
        """

#configfile: "{OutputDir}/snakemake/config/stringtie_config.yml"
GTF = config["GTFfile"]
SAMPLES = config["sample"]
OutputDir = config["outputdir"]
#workDir = config["appDir"]
#(SAMPLES,) = glob_wildcards("{OutputDir}/snakemake/result/3_HISAT2_aligned/{sample}.bam")
READS=["1","2"]
rule all:
	input:
		expand("{OutputDir}/snakemake/result/3_HISAT2_aligned/{sample}.bam.bai",sample=SAMPLES,OutputDir=OutputDir),
		expand("{OutputDir}/snakemake/result/4_stringtie/{sample}.gtf",sample=SAMPLES,OutputDir=OutputDir),
		expand("{OutputDir}/snakemake/result/5_stringtie_merge/stringtie_merged.gtf",OutputDir=OutputDir)
rule samtools:
	input:
		"{OutputDir}/snakemake/result/3_HISAT2_aligned/{sample}.bam"
	output:
		"{OutputDir}/snakemake/result/3_HISAT2_aligned/{sample}.bam.bai"
	conda:
		"env/samtool.yaml"
	shell:
		"samtools index -@ 12 -b {input} {output}" 
rule stringtie:
	input:
		"{OutputDir}/snakemake/result/3_HISAT2_aligned/{sample}.bam",
		"{OutputDir}/snakemake/result/3_HISAT2_aligned/{sample}.bam.bai"
	output:
		"{OutputDir}/snakemake/result/4_stringtie/{sample}.gtf"
	message:
		"Begin stringtie transcript assembly"
	conda:
		"env/samtool.yaml"
	shell:
		"stringtie -p 12 -G {GTF} -o {output[0]} -i {input[0]}"
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
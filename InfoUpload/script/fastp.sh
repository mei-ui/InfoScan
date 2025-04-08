#!/bin/bash
#snakemake -F -s $1 --config workdir=$2

Help() {
echo
"bash fastp.sh -i SRR1663493.1_1 -t SRR1663494.1_1 -c fastp_config.yml" 

}

                                   
usage() {                                      # Function: Print a help message.
  echo "Usage: $0 [ -h help] [ -i read1 ] [ -t read2 ] [ -c fastp_config.yml ]" 1>&2
}
exit_abnormal() {                              # Function: Exit with error.
  usage
  exit 1
}


while getopts ":h:i:c:" options; do              
                                                                                           
  case "${options}" in
    h|:)
      usage
      Help
      exit 0
      ;;
    i)                                         
      workDir=${OPTARG}
      if ! [[ -n $$workDir ]] ; then                      
        echo "You didn't set the read1 sample"                                
      fi                          
      ;;
    c)                                         
      outDir=${OPTARG}
      if ! [[ -n $outDir ]] ; then                      
        echo "You didn't set the output dir"                                
      fi                          
      ;;
    \?) # incorrect option
      echo "Error: -${OPTARG} Invalid option"
      exit_abnormal
      ;;
  esac
done

shift $(($OPTIND - 1))
export PATH=~/miniconda3/bin:$PATH
echo -e "activate conda envs"
conda_activate="source activate InfoScan"
eval $conda_activate

echo -e "==Starting a fastp job...==:\n"

# #fastqc
# echo -e "fastqc:${fastp_array[23]}"
# if [[ "${fastp_array[23]}" =~ "true"  ]] ;
# then
#     echo -e "Getting fastqc report...\n"
#     snakemake -s snakemake/Snakefile --cores 8 snakemake/qc/read1/"${read1}".html
#     snakemake -s snakemake/Snakefile --cores 8 snakemake/qc/read2/"${read2}".html
# fi

#cutadapt-se
#echo -e "ad_remove_read1:${fastp_array[33]}"
#if [[ "${fastp_array[33]}" =~ "true"  ]] ;
#then
#chmod 755 $workDir/*
#which snakemake
#export PATH=/Users/meisq/miniconda3/envs/InfoScan/bin/snakemake:$PATH
#chmod 755 /Users/meisq/miniconda3/envs/InfoScan/bin/snakemake
#chmod 755 /Users/meisq/miniconda3/envs/InfoScan/lib/python3.6/site-packages/snakemake/persistence.py
# 从 YAML 配置中提取 thread 参数（假设已通过工具解析成变量）
thread=$(grep 'thread:' $outDir/snakemake/config/upload_fastq_config.yml | awk '{print $2}' | tr -d '"')  # 去除引号

# 验证是否为整数
if [[ "$thread" =~ ^[0-9]+$ ]]; then
    echo "Valid thread value: $thread"
else
    echo "Error: 'thread' must be an integer." >&2
    exit 1
fi

# 后续使用时强制转为整数（可选）
thread=$((thread))
echo -e "Getting reads trimmed fastq...\n"
library=`grep -Eo "paired-end" $outDir/snakemake/config/upload_fastq_config.yml`
cat $outDir/snakemake/sample.yaml $outDir/snakemake/config/upload_fastq_config.yml > $outDir/snakemake/config/fastp_config.yaml && 
if [ $library = "paired-end" ]; \
  then snakemake -s $workDir/fastp.smk --configfile $outDir/snakemake/config/fastp_config.yaml --unlock;snakemake -s $workDir/fastp.smk --cores $thread --configfile $outDir/snakemake/config/fastp_config.yaml; \
  else snakemake -s $workDir/fastp2.smk --cores 4 --use-conda;fi
#fi

#emit done signal
echo -e "==All done!=="

#!/bin/bash
#snakemake -F -s $1 --config workdir=$2

export TERM=xterm
source ~/miniconda3/etc/profile.d/conda.sh 

Help() {
echo
"bash hisat2.sh -i SRR1663493.1_1 -t SRR1663494.1_1 -c hisat2_config.yml" 

}

                                   
usage() {                                      # Function: Print a help message.
  echo "Usage: $0 [ -h help] [ -i log file ] [ -t workDir ] [ -c hisat2_config.yml ]  [ -o outDir ]" 1>&2
}
exit_abnormal() {                              # Function: Exit with error.
  usage
  exit 1
}


while getopts ":h:i:t:c:o:" options; do              
                                                                                           
  case "${options}" in
    h|:)
      usage
      Help
      exit 0
      ;;
    c)                                         
      config=${OPTARG}
      if ! [[ -n $config ]] ; then                      
        echo "You didn't set the config yml file"                                
      fi                          
      ;;
    t)                                         
      workDir=${OPTARG}
      if ! [[ -n $workDir ]] ; then                      
        echo "You didn't set the workDir"                                
      fi                          
      ;;
    o)                                         
      outDir=${OPTARG}
      if ! [[ -n $outDir ]] ; then                      
        echo "You didn't set the outDir"                                
      fi                          
      ;;
    i)                                         
      log=${OPTARG}
      if ! [[ -n $log ]] ; then                      
        echo "You didn't set the log file"                                
      fi                          
      ;;
    \?) # incorrect option
      echo "Error: -${OPTARG} Invalid option"
      exit_abnormal
      ;;
  esac
done

shift $(($OPTIND - 1))

hisat2_config=$(cat $config) 
hisat2_array=($(echo $hisat2_config | tr ":" "\n"))


echo -e "==Starting a hisat2 job...==:\n" > $log

echo -e "activate conda envs"
conda activate InfoScan

# 从 YAML 配置中提取 thread 参数（假设已通过工具解析成变量）
thread=$(grep 'thread:' $config | awk '{print $2}' | tr -d '"')  # 去除引号

# 验证是否为整数
if [[ "$thread" =~ ^[0-9]+$ ]]; then
    echo "Valid thread value: $thread"
else
    echo "Error: 'thread' must be an integer." >&2
    exit 1
fi

# 后续使用时强制转为整数（可选）
thread=$((thread))

cat $outDir/snakemake/sample.yaml $config > $outDir/snakemake/config/InfoAssembly_config.yaml && 
cd $workDir/../MacOS
snakemake -s $workDir/hisat2.smk --configfile $outDir/snakemake/config/InfoAssembly_config.yaml --unlock  && \
echo -e "Getting transcript alignment...\n"
snakemake -s $workDir/hisat2.smk  --cores $thread --configfile $outDir/snakemake/config/InfoAssembly_config.yaml --rerun-incomplete >> $log 2>&1

#emit done signal
echo -e "==All done!==" >> $log
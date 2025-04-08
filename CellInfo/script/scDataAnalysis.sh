#!/bin/bash
#snakemake -F -s $1 --config workdir=$2

Help() {
echo
"bash scRNAdataAnalysis.sh -i log -c DataAnalysis_config.yml" 

}

                                   
usage() {                                      # Function: Print a help message.
  echo "Usage: $0 [ -h help] [ -i log ] [ -c DataAnalysis_config.yml ]" 1>&2
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
    c)                                         
      config=${OPTARG}
      if ! [[ -n $config ]] ; then                      
        echo "You didn't set the config yml file"                                
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

DataAnalysis_config=$(cat $config) 
scRNAdataAnalysis_array=($(echo $DataAnalysis_config | tr ":" "\n"))

# for i in "${scRNAdataAnalysis_array[@]}"
# do
#    echo $i
# done

#echo -e "clear snakemake history"
#rm -rf .snakemake
echo -e "activate conda envs"
conda_activate="source activate InfoScan3"
eval $conda_activate

echo -e "==Starting a scRNAdataAnalysis job...==:\n" > $log

# #fastqc
# echo -e "fastqc:${scRNAdataAnalysis_array[23]}"
# if [[ "${scRNAdataAnalysis_array[23]}" =~ "true"  ]] ;
# then
#     echo -e "Getting fastqc report...\n"
#     snakemake -s snakemake/Snakefile --cores 8 snakemake/qc/read1/"${read1}".html
#     snakemake -s snakemake/Snakefile --cores 8 snakemake/qc/read2/"${read2}".html
# fi
echo -e "Getting lncRNA Find...\n"
snakemake -s snakemake/scRNAdataAnalysis.smk --unlock
snakemake -s snakemake/scRNAdataAnalysis.smk --cores 1 --use-conda --configfile $config >> $log 2>&1
#emit done signal
echo -e "==All done!==" >> $log
#!/bin/bash
#snakemake -F -s $1 --config workdir=$2

Help() {
echo
"bash scRNAdataAnalysis.sh -i log -c FunctionAnalysis_config.yml" 

}

                                   
usage() {                                      # Function: Print a help message.
  echo "Usage: $0 [ -h help] [ -i log ] [ -c FunctionAnalysis_config.yml ]" 1>&2
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

FunctionAnalysis_config=$(cat $config) 
scRNAdataAnalysis_array=($(echo $FunctionAnalysis_config | tr ":" "\n"))

# for i in "${scRNAdataAnalysis_array[@]}"
# do
#    echo $i
# done

#echo -e "clear snakemake history"
#rm -rf .snakemake
#echo -e "activate conda envs"
##eval $conda_activate

echo -e "==Starting a FuncScan job...==:\n" > $log
outdir=${config%/snakemake/config*}
cat $outDir/snakemake/config/NovelScan_config.yml $config > $outDir/snakemake/config/FuncScan.yml && 
snakemake -s snakemake/FunctionAnalysis.smk --cores 1 -F --use-conda --conda-frontend conda --configfile $outDir/snakemake/config/FuncScan.yml>> $log 2>&1
#emit done signal
echo -e "==All done!==" >> $log
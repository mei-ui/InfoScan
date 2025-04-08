#!/bin/bash
#snakemake -F -s $1 --config workdir=$2
export TERM=xterm
source ~/miniconda3/etc/profile.d/conda.sh 
Help() {
echo
"bash stringtie.sh -i SRR1663493.1_1 -t SRR1663494.1_1 -c stringtie_config.yml" 

}

                                   
usage() {                                      # Function: Print a help message.
  echo "Usage: $0 [ -h help] [ -i work_path ] [ -c stringtie_config.yml ]" 1>&2
}
exit_abnormal() {                              # Function: Exit with error.
  usage
  exit 1
}


while getopts ":h:i:c:t:" options; do
                                                                                           
  case "${options}" in
    h|:)
      usage
      Help
      exit 0
      ;;
    i)                                         
      workDir=${OPTARG}
      if ! [[ -n $workDir ]] ; then                      
        echo "You didn't set the work path"
      fi                          
      ;;
    c)                                         
      config=${OPTARG}
      if ! [[ -n $config ]] ; then                      
        echo "You didn't set the config yml file"                                
      fi                          
      ;;
    t)                                         
      outDir=${OPTARG}
      if ! [[ -n $outDir ]] ; then                      
        echo "You didn't set the outdir file"                                
      fi                          
      ;;
    \?) # incorrect option
      echo "Error: -${OPTARG} Invalid option"
      exit_abnormal
      ;;
  esac
done

shift $(($OPTIND - 1))

echo -e "==Starting a stringtie job...==:\n"
echo -e "conda activate InfoScan"
conda activate InfoScan

echo -e "Getting stringtie...\n"
cat $outDir/snakemake/sample_bam.yaml $config > $outDir/snakemake/config/Stringrie_config.yaml 
cd $workDir/../MacOS
snakemake -s $workDir/stringtie.smk --configfile $outDir/snakemake/config/Stringrie_config.yaml --unlock
echo -e "Getting transcript alignment...\n"
snakemake -s $workDir/stringtie.smk --cores 1 --configfile $outDir/snakemake/config/Stringrie_config.yaml  --rerun-incomplete
#fi

#emit done signal
echo -e "==All done!=="

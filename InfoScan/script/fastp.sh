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


while getopts ":h:i:t:c:" options; do              
                                                                                           
  case "${options}" in
    h|:)
      usage
      Help
      exit 0
      ;;
    i)                                         
      read1=${OPTARG}
      if ! [[ -n $read1 ]] ; then                      
        echo "You didn't set the read1 sample"                                
      fi                          
      ;;
    t)                                         
      read2=${OPTARG}
      if ! [[ -n $read2 ]] ; then                      
        echo "You didn't set the read2 sample"                                
      fi                          
      ;;
    c)                                         
      config=${OPTARG}
      if ! [[ -n $config ]] ; then                      
        echo "You didn't set the config yml file"                                
      fi                          
      ;;
    \?) # incorrect option
      echo "Error: -${OPTARG} Invalid option"
      exit_abnormal
      ;;
  esac
done

shift $(($OPTIND - 1))

fastp_config=$(cat $config) 
fastp_array=($(echo $fastp_config | tr ":" "\n"))

# for i in "${fastp_array[@]}"
# do
#    echo $i
# done

#echo -e "clear snakemake history"
#rm -rf .snakemake

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
echo -e "Getting reads trimmed fastq...\n"
snakemake -s snakemake/fastp.smk --cores 8 
#fi

#emit done signal
echo -e "==All done!=="
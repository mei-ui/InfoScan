#!/bin/bash
#snakemake -F -s $1 --config workdir=$2

Help() {
echo
"bash SpatialScan.sh -i log -c SpatialScan_config.yml" 

}

                                   
usage() {                                      # Function: Print a help message.
  echo "Usage: $0 [ -h help] [ -i log ] [ -c SpatialScan_config.yml ]" 1>&2
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

SpatialScan_config=$(cat $config) 
SpatialScan_array=($(echo $SpatialScan_config | tr ":" "\n"))

echo -e "==Starting a SpatialScan job...==:\n" > $log

snakemake -s snakemake/SpatialScan.smk --cores 1 -F --use-conda --conda-frontend conda --configfile $config>> $log 2>&1
#emit done signal
echo -e "==All done!==" >> $log
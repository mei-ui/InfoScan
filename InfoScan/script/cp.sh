#! /bin/bash

if [[ $# -ne 3 ]]
then
    echo 'Usage: ./'$0 ' sample_path work_path log_file'
    exit 1
fi

sample_path=$1
work_path=$2
log_file=$3


echo -e "clear fsatp.sh history"
echo -e "Your fastp job is starting, please wait..." > $log_file

#.fastq.gz
if [[  -d "$work_path/data" ]]; then
    rm -rf $work_path/data
fi

if [[ ! -d "$work_path/data" ]]; then
    mkdir $work_path/data
fi



if [[ "$sample_path" != "" ]] ;
then
    echo -e "$(tput setaf 2)read1: unzipping and redirecting fastq file...\n$(tput sgr 0)"
    cp -r $sample_path/*  $work_path/data &
fi


wait


echo -e "Starting snakemake workflow..."
#snakemake
bash snakemake/script/fastp.sh -i $work_path -c fastp_config.yml >> $log_file 2>&1 &

wait

echo -e "Finished snakemake workflow!" >> $log_file
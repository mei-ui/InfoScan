#! /bin/bash
export TERM=xterm
source ~/miniconda3/etc/profile.d/conda.sh 

if [[ $# -ne 4 ]]
then
    echo 'Usage: ./'$0 ' sample_path work_path out_path log_file'
    exit 1
fi

sample_path=$1
work_path=$2
out_path=$3
log_file=$4


echo -e "clear stringtie.sh history"
echo -e "Your stringtie job is starting, please wait..." > $log_file

#.bam
if [[ ! -d "$out_path" ]]; then
    mkdir $out_path
fi

if [[  -d "$out_path/snakemake/result/3_HISAT2_aligned" ]]; then
    rm -rf $out_path/snakemake/result/3_HISAT2_aligned/*
fi
if [[ ! -d "$out_path/snakemake/result/3_HISAT2_aligned" ]]; then
    mkdir -p $out_path/result
fi


if [[ "$sample_path" != "" ]] ;
then
    echo -e "copying bam files to $out_path/snakemake/result/3_HISAT2_aligned"
    cp -r $sample_path/*  $out_path/snakemake/result/3_HISAT2_aligned &
fi

wait

echo "Generating sample manifest..."
find "$out_path/snakemake/result/3_HISAT2_aligned" -name '*.bam'  -exec basename {} \; | \
    sed "s/.bam//" | \
    awk 'BEGIN{print "sample: ["} {printf "\047%s\047,\n", $0} END{print "]"}' | \
    sed '$ s/,$//' > "$out_path/snakemake/sample_bam.yaml"

echo -e "conda activate InfoScan"
conda activate InfoScan

echo -e "Starting snakemake workflow..."
#snakemake
bash $work_path/script/stringtie.sh -i $work_path -c $out_path/snakemake/config/upload_bam_config.yml -t $out_path >> $log_file 2>&1 &
wait
echo -e "Finished snakemake workflow!" >> $log_file

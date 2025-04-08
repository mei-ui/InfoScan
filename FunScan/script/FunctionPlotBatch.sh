#!/bin/bash
#snakemake -F -s $1 --config workdir=$2
config=$1
log=$2

echo -e "==Starting a FunctionPlotBatch job...==:\n"  > $log
outdir=${config%/snakemake/config*}
cat $outDir/snakemake/config/NovelScan_config.yaml $config > $outDir/snakemake/config/FuncScan2.yaml && 
snakemake -s snakemake/FunctionPlotBatch.smk --cores 1 --use-conda --conda-frontend conda --configfile $outDir/snakemake/config/FuncScan2.yaml >> $log 2>&1
#emit done signal
echo -e "==All done!==" >> $log
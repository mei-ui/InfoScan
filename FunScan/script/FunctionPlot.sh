#!/bin/bash
log=$2
config=$1
#echo -e "clear snakemake history"
#rm -rf .snakemake
echo -e "activate conda envs"
conda_activate="source activate InfoScan3"
eval $conda_activate

echo -e "==Starting a Function analysis job...==:\n" > $log

echo -e "Getting plot...\n" >> $log
snakemake -s snakemake/FunctionPlot.smk --core 8  -F --configfile $config >> $log 2>&1
#emit done signal
echo -e "==All done!==" >> $log
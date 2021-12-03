from snakemake.shell import shell
import re


n = len(snakemake.input)
assert n == 2, "Input must contain 2 (paired-end) elements."

quality_phred = snakemake.params.get("quality_phred", "")
adapters = snakemake.params.get("adapters", "")
thread = snakemake.params.get("thread", "")
adapter_read1 = snakemake.params.get("adapter_read1", "")
adapter_read2 = snakemake.params.get("adapter_read2", "")


assert (
    adapters != ""
), "No options provided to cutadapt. Please use 'params: adapters=' or 'params: extra='."

pattern = re.compile('[agctnAGCTN]+')
match1 = bool(re.match(pattern,adapter_read1))
match2 = bool(re.match(pattern,adapter_read2))

if match1 and match2:        
    print('start running fastp process...')
    shell(
    "fastp"
    " {snakemake.params.adapters}"
    " -q {snakemake.params.quality_phred}"
    " --thread {snakemake.params.thread}"
    " -h {snakemake.output.html}"
    " -j {snakemake.output.json}"
    " -i {snakemake.input[0]}"
    " -I {snakemake.input[1]}"
    " -o {snakemake.output[0]}"
    " -O {snakemake.output[1]}")   
else:
    print('no adapters are provided, auto-detection for adapter...')
    shell(
    "fastp"
    " -q {snakemake.params.quality_phred}"
    " --thread {snakemake.params.thread}"
    " -h {snakemake.output.html}"
    " -j {snakemake.output.json}"
    " -i {snakemake.input[0]}"
    " -I {snakemake.input[1]}"
    " -o {snakemake.output[0]}"
    " -O {snakemake.output[1]}" )   
  


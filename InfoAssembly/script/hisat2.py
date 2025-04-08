from snakemake.shell import shell
import re

n = len(snakemake.input)
if n==2:
	assert n == 2, "Input must contain 2 (paired-end) elements."
	
	thread = snakemake.params.get("thread", "")
	index = snakemake.params.get("index", "")
	params1 = snakemake.params.get("params","")
	
	
	print('start running hisat2 process...')
	shell(
	"hisat2"
	" -p {snakemake.params.thread}"
	" -x {snakemake.params.index}"
	" -1 {snakemake.input[0]}"
	" -2 {snakemake.input[1]}"
	" --dta"
	" --new-summary"
	" -S {snakemake.output[0]}"
	" {snakemake.params.params1} 2> {snakemake.log}")
if n == 1 :
	assert n == 1, "Input must contain 1 (single-end) elements."

	thread = snakemake.params.get("thread", "")
	index = snakemake.params.get("index", "")
	params1 = snakemake.params.get("params","")
	
	
	print('start running hisat2 process...')
	shell(
	"hisat2"
	" -p {snakemake.params.thread}"
	" -x {snakemake.params.index}"
	" -U {snakemake.input}"
	" --dta"
	" --new-summary"
	" -S {snakemake.output[0]}"
	" {snakemake.params.params1} 2> {snakemake.log}")
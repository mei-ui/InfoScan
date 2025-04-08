#!/bin/bash
set -eo pipefail

export TERM=xterm
source ~/miniconda3/etc/profile.d/conda.sh 

if [[ $# -ne 4 ]]; then
    echo "Usage: $0 sample_path work_path out_path log_file"
    exit 1
fi

sample_path="$1"
work_path="$2"
out_path="$3"
log_file="$4"

mkdir -p "$out_path/snakemake/logs"
rm -f "$out_path/snakemake/logs/upload_fastq.log"
touch "$out_path/snakemake/logs/upload_fastq.log"

echo -e "clear fsatp.sh history"
echo -e "Your fastp job is starting, please wait..." > "$log_file"

if [[ -d "$out_path/data" ]]; then
    echo "Removing old data directory..."
    rm -rf "$out_path/data"
fi
mkdir -p "$out_path/data"

if [[ -n "$sample_path" ]]; then
    echo -e "$(tput setaf 2)read1: unzipping and redirecting fastq file...$(tput sgr 0)"
    if ! cp -r "$sample_path"/* "$out_path/data/"; then
        echo "Error: Failed to copy sample files!" >&2
        exit 1
    fi
fi

shopt -s nullglob

rename_files() {
    local pattern=$1
    local from_suffix=$2
    local to_suffix=$3
    for i in "$out_path/data/"*$pattern; do
        file=$(basename "$i")
        sample="${file%%_*}"
        if [[ "$file" == "${sample}_${from_suffix}" ]]; then
            new_path="$out_path/data/${sample}_${to_suffix}"
            if [[ -e "$new_path" ]]; then
                echo "Notice: Target file '$new_path' exists. Skipping."
                continue
            fi
            mv "$i" "$new_path"
        fi
    done
}

rename_patterns=(
    "_R1.fastq.gz:R1.fastq.gz:_1.fastq.gz"
    "_R1.fq.gz:R1.fq.gz:_1.fastq.gz"
    "_r1.fq.gz:r1.fq.gz:_1.fastq.gz"
    "_r1.fastq.gz:r1.fastq.gz:_1.fastq.gz"
    "_1.fq.gz:1.fq.gz:_1.fastq.gz"
    "_R1.fastq:R1.fastq:_1.fastq.gz"
    "_R1.fq:R1.fq:_1.fastq.gz"
    "_r1.fq:r1.fq:_1.fastq.gz"
    "_r1.fastq:r1.fastq:_1.fastq.gz"
    "_1.fastq:1.fastq:_1.fastq.gz"
    "_1.fq:1.fq:_1.fastq.gz"
)

for pattern in "${rename_patterns[@]}"; do
    IFS=":" read -r match_pattern from to <<< "$pattern"
    rename_files "$match_pattern" "$from" "$to"
done

for ext in fastq fq; do
    for i in "$out_path/data/"*."$ext"; do
        if [[ -f "$i" ]]; then
            echo "Compressing $i..."
            pigz -p 8 "$i"
        fi
    done
done

echo "Generating sample manifest..."
find "$out_path/data" -name '*_1.fastq.gz' ! -name '._*_1.fastq.gz' -exec basename {} \; | \
    sed "s/_1.fastq.gz//" | \
    awk 'BEGIN{print "sample: ["} {printf "\047%s\047,\n", $0} END{print "]"}' | \
    sed '$ s/,$//' > "$out_path/snakemake/sample.yaml"

echo -e "activate conda envs"
conda activate InfoScan

if ! command -v stringtie >/dev/null 2>&1; then
    echo "Error: stringtie not found in PATH!" >&2
    exit 1
fi

echo -e "Starting snakemake workflow..."
if ! bash "$work_path/script/fastp.sh" -i "$work_path" -c "$out_path" >> "$log_file" 2>&1; then
    echo "Error: Snakemake workflow failed!" >&2
    exit 1
fi

echo -e "Finished snakemake workflow!" >> "$log_file"
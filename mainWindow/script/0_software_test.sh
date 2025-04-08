#!/bin/bash
set -eo pipefail

Help() {
    echo "Usage: $0 [ -h help ] [ -y YAML ] [ -d WORKDIR ]"
    echo "Example:"
    echo "  bash $0 -y env.yaml -d /path/to/InfoScan.app/Contents/Resources"
}

init_pkgmgr() {
    declare -a micromamba_paths=(
        "${HOME}/micromamba/micromamba"
        "/usr/local/bin/micromamba" 
        "${HOME}/.local/bin/micromamba"
        "/opt/micromamba/bin/micromamba"
    )

    for path in "${micromamba_paths[@]}"; do
        if [ -x "$path" ]; then
            echo "INFO: Found Micromamba at $path"
            eval "$("$path" shell hook -s bash)"
            export MICROMAMBA_ROOT_PREFIX="${path%/*}"
            return 0
        fi
    done

    declare -a conda_paths=(
        "${HOME}/anaconda3"
        "${HOME}/miniconda3"
        "/opt/anaconda3"
        "/opt/miniconda3"
        "/usr/local/anaconda3"
    )

    for path in "${conda_paths[@]}"; do
        local conda_sh="${path}/etc/profile.d/conda.sh"
        if [ -f "$conda_sh" ]; then
            echo "INFO: Found Conda at $path"
            source "$conda_sh"
            return 0
        fi
    done

    echo "ERROR: Neither Micromamba nor Conda found!"
    echo "Checked Micromamba paths:"
    printf '  - %s\n' "${micromamba_paths[@]}"
    echo "Checked Conda paths:"
    printf '  - %s\n' "${conda_paths[@]}"
    exit 1
}

check_env_exists() {
    local env_name="$1"
    if command -v micromamba &>/dev/null; then
        micromamba env list --json | grep -q "\"name\": \"$env_name\""
    else
        conda env list | grep -qw "$env_name"
    fi
}

get_env_name() {
    local yaml_file="$1"
    name_line=$(grep -E '^[[:space:]]*name:' "$yaml_file" | head -n1)
    [ -z "$name_line" ] && { echo "ERROR: Missing 'name' field in YAML"; exit 1; }
    env_name=$(echo "$name_line" | awk -F: '{print $2}' | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//' -e "s/^['\"]//" -e "s/['\"]$//")
    [ -z "$env_name" ] && { echo "ERROR: Empty name value in YAML"; exit 1; }
    echo "$env_name"
}

validate_workdir() {
    local dir="$1"
    [ -d "$dir/snakemake/script" ] || {
        echo "ERROR: Invalid workdir structure, missing subdirs: $dir/snakemake/{script,env}"
        exit 1
    }
}

main() {
    init_pkgmgr

    local yaml="" workDir="" env_name=""
    while getopts ":h:y:d:" opt; do
        case $opt in
            h) Help; exit 0 ;;
            y) yaml=$(realpath  "$OPTARG") ;;
            d) workDir=$(realpath "$OPTARG") ;;
            \?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
            :) echo "Option -$OPTARG requires argument" >&2; exit 1 ;;
        esac
    done

    [ -z "$workDir" ] && { echo "ERROR: -d parameter required" >&2; exit 1; }
    validate_workdir "$workDir"

    if [ -n "$yaml" ]; then
        env_name=$(get_env_name "$yaml") || exit 1
    else
        echo "WARN: Using default env name 'InfoScan'"
        env_name="InfoScan"
    fi

    if check_env_exists "$env_name"; then
        echo "[$(date +%T)] Activating existing environment: $env_name"
    else
        [ -n "$yaml" ] || { echo "ERROR: New env requires YAML file"; exit 1; }
        echo "[$(date +%T)] Creating environment from $yaml..."
        if command -v micromamba &>/dev/null; then
            micromamba env create -f "$yaml" -y
        else
            conda env create -f "$yaml"
        fi
    fi

    local retries=3
    for ((i=1; i<=retries; i++)); do
        if { command -v micromamba &>/dev/null && micromamba activate "$env_name"; } || conda activate "$env_name"; then
            break
        elif [ $i -eq $retries ]; then
            echo "ERROR: Failed to activate $env_name after $retries attempts"
            exit 1
        else
            echo "[Retry $i/$retries] Attempting environment repair..."
            if command -v micromamba &>/dev/null; then
                micromamba update -f "$yaml" --prune -y
            else
                conda env update -f "$yaml" --prune
            fi
            sleep $((i * 2))
        fi
    done

    echo "[$(date +%T)] Setting executable permissions..."
    find "$workDir/snakemake/script" -type f -exec chmod -v a+x {} +

    # 修改后的依赖声明
    deps=(
        "stringtie:bioconda"
        "shyaml:conda-forge"
    )

    for entry in "${deps[@]}"; do
        dep=$(echo "$entry" | cut -d: -f1)
        channel=$(echo "$entry" | cut -d: -f2)
        if ! command -v "$dep" >/dev/null; then
            echo "[$(date +%T)] Installing $dep from $channel..."
            if command -v micromamba &>/dev/null; then
                micromamba install -y -c "$channel" "$dep"
            else
                conda install -y -c "$channel" "$dep"
            fi
        fi
    done

    echo "[$(date +%T)] Environment validation passed!"
    echo "Active environment info:"
    if command -v micromamba &>/dev/null; then
        micromamba info
    else
        conda info
    fi
}

trap 'echo "[ERROR] Failed at line $LINENO: $BASH_COMMAND"; exit 1' ERR

main "$@"
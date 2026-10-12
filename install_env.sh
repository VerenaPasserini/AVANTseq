#!/bin/bash
# Create (or update) the "avantseq" conda environment with all tools needed by AVANTseq.
# Usage:  bash install_env.sh
set -e
cd "$(dirname "$0")"
CONDA=$(command -v mamba || command -v conda) || { echo "conda or mamba not found"; exit 1; }

if $CONDA env list | grep -qE "^avantseq\s"; then
    $CONDA env update -n avantseq -f envs/avantseq.yaml --prune
else
    $CONDA env create -f envs/avantseq.yaml
fi

echo; echo "Done. Activate with:  conda activate avantseq"

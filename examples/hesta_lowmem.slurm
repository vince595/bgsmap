#!/usr/bin/env bash
#SBATCH --job-name=bgsmap-hesta
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --output=bgsmap-hesta-%j.log

# Submit from an environment where bgsmap 1.2.0 is already on PATH.
# sbatch examples/hesta_lowmem.slurm GSS H5AD BURDEN CACHE_DIR NEW_RESULT_DIR
set -euo pipefail

if [[ $# -ne 5 ]]; then
    echo "Usage: sbatch $0 GSS.feather SPATIAL.h5ad BURDEN.tsv.gz CACHE_DIR NEW_RESULT_DIR" >&2
    exit 2
fi
gss=$1
h5ad=$2
burden=$3
cache=$4
out=$5

if [[ -e "$out" ]]; then
    echo "Choose a new output directory: $out" >&2
    exit 2
fi

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cpus=${SLURM_CPUS_PER_TASK:-8}
workers=$((cpus > 1 ? cpus - 1 : 1))

bgsmap cache-gss --gss "$gss" --out "$cache" --projection-size 4096
bgsmap-lowmem run \
    --gss "$gss" --h5ad "$h5ad" --annotation substructure \
    --burden "$burden" --gene-col gene_symbol --p-col Pvalue_Burden \
    --score-mode logp --duplicate-policy bonferroni --permutations 0 --seed 42 \
    --threads "$workers" --block-size 512 --gss-cache "$cache" \
    --skip-plots --out "$out"

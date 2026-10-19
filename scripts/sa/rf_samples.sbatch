#!/bin/bash
#SBATCH --job-name=vfb-rf
#SBATCH --account=def-nicoleli
#SBATCH --time=0-03:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000M
#SBATCH --gpus-per-node=h100:1
#SBATCH --output=logs/rf_%A_%a.out
#SBATCH --error=logs/rf_%A_%a.err
# One array task runs CHUNK samples in parallel on one GPU node.
# Submit from the repo root, e.g. for 512 samples in chunks of 16:
#   sbatch --array=0-31 scripts/sa/rf_samples.sbatch
set -euo pipefail
module load StdEnv/2023 gcc/12.3 python/3.11 cuda/12.2
export OMP_NUM_THREADS=1
CHUNK=${CHUNK:-16}
START=$(( SLURM_ARRAY_TASK_ID * CHUNK + 1 ))
END=$(( START + CHUNK - 1 ))
python3 scripts/sa/rf_run_chunk.py --start "${START}" --end "${END}" --jobs "${SLURM_CPUS_PER_TASK}"

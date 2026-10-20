#!/bin/bash
#SBATCH --job-name=vfb-verify
#SBATCH --account=def-nicoleli
#SBATCH --time=0-04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --gpus-per-node=h100:1
#SBATCH --mem=32000M
#SBATCH --output=logs/verify_%j.out
#SBATCH --error=logs/verify_%j.err
# Runs scripts/verify_rules.py as one job and archives it like --keep.
# Submit from the repo root, e.g.
#   sbatch --export=ALL,KEEP_NAME=verify-migration-2x,KEEP_WHY="migration rerun at 2x v",VERIFY_ARGS="--only migration2x --seeds 5 --jobs 8" scripts/verify.sbatch
set -euo pipefail
source scripts/job_defaults.env
: "${IVDBM_RESULTS_DIR:=${HOME}/projects/def-nicoleli/results/vfb-abm-results}"
module load StdEnv/2023 python/3.11 cuda
export OMP_NUM_THREADS=1

COMMIT="$(git describe --always --dirty)"
OUT="output/verify_${SLURM_JOB_ID}"
mkdir -p "${OUT}" logs
echo "commit ${COMMIT}  args: ${VERIFY_ARGS:-}"

python3 scripts/verify_rules.py ${VERIFY_ARGS:-} \
    --jobs "${SLURM_CPUS_PER_TASK}" --out "${OUT}" | tee "${OUT}/verdict.txt"

if [[ -n "${KEEP_NAME:-}" ]]; then
  RUN_FOLDER="${IVDBM_RESULTS_DIR}/$(date +%F)_${SLURM_JOB_ID}_${KEEP_NAME}"
  mkdir -p "${RUN_FOLDER}/summary"
  cp -r "${OUT}" "${RUN_FOLDER}/raw"
  cp "${OUT}/verdict.txt" "${RUN_FOLDER}/summary/"
  cat > "${RUN_FOLDER}/manifest.txt" <<EOF
date: $(date '+%Y-%m-%d %H:%M')
cluster: ${CC_CLUSTER:-$(hostname)}
job_id: ${SLURM_JOB_ID}
commit: ${COMMIT}
script: scripts/verify_rules.py ${VERIFY_ARGS:-}
why: ${KEEP_WHY:-}
EOF
  echo "Archived to ${RUN_FOLDER}"
fi

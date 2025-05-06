#!/bin/bash

#SBATCH --job-name=ldsc_array
#SBATCH --partition=interruptible_cpu
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4000
#SBATCH --cpus-per-task=1
#SBATCH --time=0-02:00
#SBATCH --array=0-433

# Navigate to the working directory
cd /cephfs/volumes/hpc_data_usr/k20051211/8e486e0b-f00f-4534-8dcc-d66a07c10fd5/scratch_tmp/s_ldsc/ldsc/ || { echo "Working directory not found!"; exit 1; }

module load anaconda3/2021.05-gcc-13.2.0
source activate /cephfs/volumes/hpc_data_usr/k20051211/8e486e0b-f00f-4534-8dcc-d66a07c10fd5/scratch_tmp/s_ldsc/ldsc/ldsc_bioconda_env

# Directories (unchanged)
SUMSTATS_DIR="/cephfs/volumes/hpc_data_usr/k20051211/8e486e0b-f00f-4534-8dcc-d66a07c10fd5/scratch_tmp/s_ldsc/ldsc/sumstats_final"
REF_LD_DIR="/cephfs/volumes/hpc_data_usr/k20051211/8e486e0b-f00f-4534-8dcc-d66a07c10fd5/scratch_tmp/s_ldsc/ldsc/sldsc_ref/1000G_EUR_Phase3_baseline"
WEIGHTS_DIR="/cephfs/volumes/hpc_data_usr/k20051211/8e486e0b-f00f-4534-8dcc-d66a07c10fd5/scratch_tmp/s_ldsc/ldsc/sldsc_ref/1000G_Phase3_weights_hm3_no_MHC"
FRQ_DIR="/cephfs/volumes/hpc_data_usr/k20051211/8e486e0b-f00f-4534-8dcc-d66a07c10fd5/scratch_tmp/s_ldsc/ldsc/sldsc_ref/1000G_Phase3_frq"
SUPERCELL_DIR="/cephfs/volumes/hpc_data_usr/k20051211/8e486e0b-f00f-4534-8dcc-d66a07c10fd5/scratch_tmp/s_ldsc/ldsc/superclusters"
OUT_DIR="/cephfs/volumes/hpc_data_usr/k20051211/8e486e0b-f00f-4534-8dcc-d66a07c10fd5/scratch_tmp/s_ldsc/ldsc/results_updated"

mkdir -p "$OUT_DIR"

COMBOS_FILE="/cephfs/volumes/hpc_data_usr/k20051211/8e486e0b-f00f-4534-8dcc-d66a07c10fd5/scratch_tmp/s_ldsc/ldsc/combos.txt"

# Make sure combos.txt exists and has enough lines for the array
TOTAL_LINES=$(wc -l < "$COMBOS_FILE")
if (( SLURM_ARRAY_TASK_ID >= TOTAL_LINES )); then
  echo "Invalid SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID. Max index: $((TOTAL_LINES - 1))"
  exit 1
fi

# Get the specific combination for this array job
TASK_COMBO=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$COMBOS_FILE")
IFS=',' read -r TRAIT SUMSTATS_FILE CELL_NAME CELL_TYPE <<<"$TASK_COMBO"

OUT_PREFIX="${OUT_DIR}/${TRAIT}_${CELL_NAME}"

ldsc.py \
  --h2 "$SUMSTATS_FILE" \
  --ref-ld-chr "${REF_LD_DIR}/baseline.,${CELL_TYPE}/baseline." \
  --w-ld-chr "${WEIGHTS_DIR}/weights.hm3_noMHC." \
  --overlap-annot \
  --frqfile-chr "${FRQ_DIR}/1000G.EUR.QC." \
  --out "$OUT_PREFIX"

echo "Task $SLURM_ARRAY_TASK_ID for $TRAIT and $CELL_NAME completed."

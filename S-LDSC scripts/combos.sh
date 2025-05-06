#!/bin/bash

SUMSTATS_DIR="/cephfs/volumes/hpc_data_usr/k20051211/8e486e0b-f00f-4534-8dcc-d66a07c10fd5/scratch_tmp/s_ldsc/ldsc/sumstats_final"
SUPERCELL_DIR="/cephfs/volumes/hpc_data_usr/k20051211/8e486e0b-f00f-4534-8dcc-d66a07c10fd5/scratch_tmp/s_ldsc/ldsc/superclusters"
COMBOS_FILE="/cephfs/volumes/hpc_data_usr/k20051211/8e486e0b-f00f-4534-8dcc-d66a07c10fd5/scratch_tmp/s_ldsc/ldsc/combos.txt"

# Overwrite existing combos file
> "$COMBOS_FILE"

for SUMSTATS_FILE in "$SUMSTATS_DIR"/*.sumstats.gz; do
  if [[ -f "$SUMSTATS_FILE" ]]; then
    TRAIT=$(basename "$SUMSTATS_FILE" .sumstats.gz)
    for CELL_TYPE in "$SUPERCELL_DIR"/*; do
      if [[ -d "$CELL_TYPE" ]]; then
        CELL_NAME=$(basename "$CELL_TYPE")
        echo "$TRAIT,$SUMSTATS_FILE,$CELL_NAME,$CELL_TYPE" >> "$COMBOS_FILE"
      fi
    done
  fi
done

echo "Generated combos.txt with $(wc -l < "$COMBOS_FILE") lines."

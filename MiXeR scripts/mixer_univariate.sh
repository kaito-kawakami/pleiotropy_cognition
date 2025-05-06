#!/bin/bash
#SBATCH --job-name=mixer_bivar
#SBATCH --output=mixer_bivar-%A_%a.txt
#SBATCH --error=mixer_bivar-%A_%a.txt
#SBATCH --time=08:00:00
#SBATCH --partition=interruptible_cpu
#SBATCH --cpus-per-task=25
#SBATCH --mem-per-cpu=32G
#SBATCH --array=1-20

# Purge modules and load Singularity
module purge

# Set paths for reference data and Singularity image.
# Here COMORMENT is assumed to be the directory holding the reference data.
export COMORMENT=/scratch_tmp/users/k20051211/network_final/mixer
export SINGULARITY_BIND="$COMORMENT/mixer/reference:/REF:ro"
export SIF=$COMORMENT/mixer/singularity
export MIXER_COMMON_ARGS="--ld-file /REF/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld --bim-file /REF/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim --threads 25"
export REP="rep${SLURM_ARRAY_TASK_ID}"
export EXTRACT="--extract /REF/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.$REP.snps"

# Use the Singularity image's Python. --home ensures your current directory is available in the container.
export PYTHON="singularity exec --home $PWD:/home $SIF/mixer.sif python"

# Define paths to your GWAS summary statistics for each trait.
export MEMORY="/REF/sumstats/memory_qc_noMHC.csv.gz"
export RT="/REF/sumstats/rt_qc_noMHC.csv.gz"
export PROSPECTIVE="/REF/sumstats/prospective_qc_noMHC.csv.gz"
export SYMBOL="/REF/sumstats/symbol_qc_noMHC.csv.gz"
export NUMERIC="/REF/sumstats/numeric_qc_noMHC.csv.gz"
export TMTB="/REF/sumstats/tmtb_qc_noMHC.csv.gz"
export VNR="/REF/sumstats/vnr_qc_noMHC.csv.gz"

# Define output directory and ensure it exists.
export output="/scratch_tmp/users/k20051211/network_final/mixer/mixer/mixer_output"
mkdir -p $output

# Loop over the traits and run fit and test analyses for each.
for trait in MEMORY RT PROSPECTIVE SYMBOL NUMERIC TMTB VNR; do
    # Using indirect expansion to retrieve the file path stored in the variable named by $trait.
    file=${!trait}
    
    # Define a prefix for output files based on the trait name.
    out_prefix=$output/${trait}

    echo "Running fit1 for trait: $trait"
    $PYTHON /tools/mixer/precimed/mixer.py fit1 $MIXER_COMMON_ARGS $EXTRACT --trait1-file $file --out ${out_prefix}.fit.${REP}.json

    echo "Running test1 for trait: $trait"
    $PYTHON /tools/mixer/precimed/mixer.py test1 $MIXER_COMMON_ARGS --trait1-file $file --load-params ${out_prefix}.fit.${REP}.json.json --out ${out_prefix}.test.${REP}
done

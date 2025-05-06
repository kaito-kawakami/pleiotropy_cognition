#!/bin/bash
# This script submits separate jobs for each phenotype pair
source settings.sh

# set the walltime & node (assumes slurm, adapt as necessary)
WT=2:00
NODE=normal 

# iterate over all unique phenotype pairs
N_PAIRS=$(wc -l pheno.pairs.txt | awk '{print $1}')

for I in $(seq 1 $N_PAIRS); do
	# extract phenotype IDs
	P1=$(awk 'NR=='$I' {print $1}' pheno.pairs.txt)
	P2=$(awk 'NR=='$I' {print $2}' pheno.pairs.txt)

	echo "$P1 $P2"

	# check if output files for these phenotypes already exist
	if [[ -f $OUT/$P1.$P2.univ ]]; then
		echo "*  output files already exist"
	else
		# if not, SUBMIT JOB! (note: this is written for slurm, adapt the line below as necessary)
		sbatch -t $WT:00 -J $P1.$P2.lava -o slurm.$P1-$P2.%A.out lava_rg.job $SCRIPTS/settings.sh $P1 $P2
	fi
done

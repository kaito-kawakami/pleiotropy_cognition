# this 'settings' file contains relevant directories and file names used for the analyses
# you need to at least adapt the MAIN path 

# key directories
MAIN=/scratch_tmp/users/k20051211/network_final 	# base directory ( ADAPT !!! )
SCRIPTS=$MAIN/lava		# scripts dir
OUT=$MAIN/lava/results	# results dir 
DATA=/scratch_tmp/users/k20051211/network_final


# names of key input files
LOCFILE="/scratch_tmp/users/k20051211/network_final/lava/block_partition/lava-partitioning/ldblock_output_size5000_max25.blocks"
INFOFILE=$MAIN/input_info_file.tsv
OVERLAP=$MAIN/sample_overlap.tsv

# reference data dir & prefix
REFDAT=$MAIN/g1000_eur

#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 4:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f3_pancancer_mutation_ByRate_slurm/slurm_U87_REST_percentile_T_ExtendMerge.out

time python find_overlap_keep_info_NOT_sep_strand_lastColMarked_keep_bfile_info.py \
-a ../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//U87/U87_REST_percentile_T_ExtendMerge.merge.bed \
-b data/pancancer_mutation.bed \
-s hg38 -p f3_pancancer_mutation_ByRate/U87/U87_REST_percentile_T_ExtendMerge_mutatioinCount.bed

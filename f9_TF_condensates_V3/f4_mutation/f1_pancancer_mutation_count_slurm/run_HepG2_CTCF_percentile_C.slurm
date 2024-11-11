#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 4:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f1_pancancer_mutation_count_slurm/slurm_HepG2_CTCF_percentile_C.out

time bedtools intersect \
-a ../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//HepG2/HepG2_CTCF_percentile_C.merge.bed \
-b data/pancancer_mutation.bed \
-c > f1_pancancer_mutation_count/HepG2/HepG2_CTCF_percentile_C_mutatioinCount.bed

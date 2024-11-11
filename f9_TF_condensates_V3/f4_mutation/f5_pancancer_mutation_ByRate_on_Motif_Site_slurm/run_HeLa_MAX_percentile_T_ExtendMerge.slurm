#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 4:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f5_pancancer_mutation_ByRate_on_Motif_Site_slurm/slurm_HeLa_MAX_percentile_T_ExtendMerge.out

time bedtools intersect \
-a /nv/vol190/zanglab/shared/Motif/sites/hg38_fimo_jarspar/results/MAX_MA0058.3.bed \
-b ../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged//HeLa/HeLa_MAX_percentile_T_ExtendMerge.merge.bed \
-wa -u > f5_pancancer_mutation_ByRate_on_Motif_Site/HeLa/HeLa_MAX_percentile_T_ExtendMerge_Motif.overlapped.bed

time bedtools merge \
-i f5_pancancer_mutation_ByRate_on_Motif_Site/HeLa/HeLa_MAX_percentile_T_ExtendMerge_Motif.overlapped.bed \
> f5_pancancer_mutation_ByRate_on_Motif_Site/HeLa/HeLa_MAX_percentile_T_ExtendMerge_Motif.overlapped.merge.bed

time python find_overlap_keep_info_NOT_sep_strand_lastColMarked_keep_bfile_info.py \
-a f5_pancancer_mutation_ByRate_on_Motif_Site/HeLa/HeLa_MAX_percentile_T_ExtendMerge_Motif.overlapped.merge.bed \
-b data/pancancer_mutation.bed \
-s hg38 -p f5_pancancer_mutation_ByRate_on_Motif_Site/HeLa/HeLa_MAX_percentile_T_ExtendMerge_mutatioinCount.bed

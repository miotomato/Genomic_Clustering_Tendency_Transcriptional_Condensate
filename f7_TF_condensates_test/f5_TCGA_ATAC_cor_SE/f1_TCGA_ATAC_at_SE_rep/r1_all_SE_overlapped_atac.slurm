#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 24:00:00
#SBATCH -p gpu
#SBATCH -A cphg_cz3d
#SBATCH -o py2_append.out


## === constitutive
python find_overlap_keep_info_NOT_sep_strand_revised.py \
-a ../data_TCGA/mynorm_TCGA-ATAC_PanCan_Log2_QuantileNorm_Counts_plus5.txt \
-b ../data/SE_hg38/all_hg38_SE.bed \
-s hg38 -p mynorm_TCGA-ATAC_PanCan_all_SE_overlapped.bed -q nonoverlapped.bed -e2 0


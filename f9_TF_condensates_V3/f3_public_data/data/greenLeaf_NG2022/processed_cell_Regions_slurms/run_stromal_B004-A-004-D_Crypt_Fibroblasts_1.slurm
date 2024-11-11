#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o processed_cell_Regions_slurms/slurm_stromal_B004-A-004-D_Crypt_Fibroblasts_1.out

zcat GSE201336_RAW/GSM6058856_B004-A-004-D_20200817_fragments.tsv.gz \
GSE201336_RAW/GSM6058855_B004-A-004-D_20200715_fragments.tsv.gz | grep -f processed_cell_ID_per_sample_per_cellType/stromal/Crypt_Fibroblasts_1/B004-A-004-D.txt > \
processed_cell_Regions_per_sample_per_cellType/stromal/Crypt_Fibroblasts_1/B004-A-004-D.bed

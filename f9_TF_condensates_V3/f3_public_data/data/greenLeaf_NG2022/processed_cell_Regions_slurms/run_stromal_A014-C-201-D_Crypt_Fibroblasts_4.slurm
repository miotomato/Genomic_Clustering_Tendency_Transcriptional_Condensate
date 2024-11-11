#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o processed_cell_Regions_slurms/slurm_stromal_A014-C-201-D_Crypt_Fibroblasts_4.out

zcat GSE201336_RAW/GSM6058823_A014-C-201-D_20200702_fragments.tsv.gz | grep -f processed_cell_ID_per_sample_per_cellType/stromal/Crypt_Fibroblasts_4/A014-C-201-D.txt > \
processed_cell_Regions_per_sample_per_cellType/stromal/Crypt_Fibroblasts_4/A014-C-201-D.bed

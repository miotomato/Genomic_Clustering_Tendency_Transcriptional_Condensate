#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o processed_cell_Regions_slurms/slurm_stromal_F091-D_Villus_Fibroblasts_WNT5B+.out

zcat GSE201336_RAW/GSM6058849_F091-D_20200817_fragments.tsv.gz | grep -f processed_cell_ID_per_sample_per_cellType/stromal/Villus_Fibroblasts_WNT5B+/F091-D.txt > \
processed_cell_Regions_per_sample_per_cellType/stromal/Villus_Fibroblasts_WNT5B+/F091-D.bed

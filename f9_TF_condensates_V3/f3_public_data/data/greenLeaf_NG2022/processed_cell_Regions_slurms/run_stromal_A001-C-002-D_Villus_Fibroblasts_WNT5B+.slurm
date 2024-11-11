#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o processed_cell_Regions_slurms/slurm_stromal_A001-C-002-D_Villus_Fibroblasts_WNT5B+.out

zcat GSE201336_RAW/GSM6058770_A001-C-002-D_20200702_fragments.tsv.gz | grep -f processed_cell_ID_per_sample_per_cellType/stromal/Villus_Fibroblasts_WNT5B+/A001-C-002-D.txt > \
processed_cell_Regions_per_sample_per_cellType/stromal/Villus_Fibroblasts_WNT5B+/A001-C-002-D.bed

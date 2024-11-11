#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o processed_cell_Regions_slurms/slurm_stromal_A002-C-201-D_Endothelial.out

zcat GSE201336_RAW/GSM6058801_A002-C-201-D_08042020_fragments.tsv.gz | grep -f processed_cell_ID_per_sample_per_cellType/stromal/Endothelial/A002-C-201-D.txt > \
processed_cell_Regions_per_sample_per_cellType/stromal/Endothelial/A002-C-201-D.bed

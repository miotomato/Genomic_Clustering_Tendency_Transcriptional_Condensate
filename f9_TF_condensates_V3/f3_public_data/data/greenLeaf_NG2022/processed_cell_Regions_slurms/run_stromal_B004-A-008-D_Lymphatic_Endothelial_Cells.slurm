#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o processed_cell_Regions_slurms/slurm_stromal_B004-A-008-D_Lymphatic_Endothelial_Cells.out

zcat GSE201336_RAW/GSM6058857_B004-A-008-D_20200817_fragments.tsv.gz | grep -f processed_cell_ID_per_sample_per_cellType/stromal/Lymphatic_Endothelial_Cells/B004-A-008-D.txt > \
processed_cell_Regions_per_sample_per_cellType/stromal/Lymphatic_Endothelial_Cells/B004-A-008-D.bed

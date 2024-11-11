#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o processed_cell_Regions_slurms/slurm_stromal_A001-C-023-D-R2_Endothelial.out

zcat GSE201336_RAW/GSM6058773_A001-C-023-D_20200214_fragments.tsv.gz \
GSE201336_RAW/GSM6058774_A001-C-023-D_20200715_fragments.tsv.gz | grep -f processed_cell_ID_per_sample_per_cellType/stromal/Endothelial/A001-C-023-D-R2.txt > \
processed_cell_Regions_per_sample_per_cellType/stromal/Endothelial/A001-C-023-D-R2.bed

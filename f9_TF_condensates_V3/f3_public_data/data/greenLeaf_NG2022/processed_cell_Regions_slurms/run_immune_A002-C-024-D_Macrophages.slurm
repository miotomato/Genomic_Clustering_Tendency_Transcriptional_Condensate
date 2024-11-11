#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o processed_cell_Regions_slurms/slurm_immune_A002-C-024-D_Macrophages.out

zcat GSE201336_RAW/GSM6058792_A002-C-024-D_20200715_fragments.tsv.gz | grep -f processed_cell_ID_per_sample_per_cellType/immune/Macrophages/A002-C-024-D.txt > \
processed_cell_Regions_per_sample_per_cellType/immune/Macrophages/A002-C-024-D.bed

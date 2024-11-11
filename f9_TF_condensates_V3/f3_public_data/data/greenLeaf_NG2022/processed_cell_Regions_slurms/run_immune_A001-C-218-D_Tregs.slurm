#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o processed_cell_Regions_slurms/slurm_immune_A001-C-218-D_Tregs.out

zcat GSE201336_RAW/GSM6058785_A001-C-218-D_20200214_fragments.tsv.gz | grep -f processed_cell_ID_per_sample_per_cellType/immune/Tregs/A001-C-218-D.txt > \
processed_cell_Regions_per_sample_per_cellType/immune/Tregs/A001-C-218-D.bed

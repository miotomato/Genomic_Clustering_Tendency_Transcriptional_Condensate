#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o processed_cell_Regions_slurms/slurm_immune_B004-A-204-D_CD8+.out

zcat GSE201336_RAW/GSM6058858_B004-A-204-D_20200702_fragments.tsv.gz | grep -f processed_cell_ID_per_sample_per_cellType/immune/CD8+/B004-A-204-D.txt > \
processed_cell_Regions_per_sample_per_cellType/immune/CD8+/B004-A-204-D.bed

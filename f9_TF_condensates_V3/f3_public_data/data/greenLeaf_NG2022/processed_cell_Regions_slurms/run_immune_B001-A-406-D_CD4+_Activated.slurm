#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o processed_cell_Regions_slurms/slurm_immune_B001-A-406-D_CD4+_Activated.out

zcat GSE201336_RAW/GSM6058853_B001-A-406-D_fragments.tsv.gz | grep -f processed_cell_ID_per_sample_per_cellType/immune/CD4+_Activated/B001-A-406-D.txt > \
processed_cell_Regions_per_sample_per_cellType/immune/CD4+_Activated/B001-A-406-D.bed

#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o processed_cell_Regions_slurms/slurm_immune_A002-C-016-D_CD4+_Activated.out

zcat GSE201336_RAW/GSM6058790_A002-C-016-D_20200715_fragments.tsv.gz | grep -f processed_cell_ID_per_sample_per_cellType/immune/CD4+_Activated/A002-C-016-D.txt > \
processed_cell_Regions_per_sample_per_cellType/immune/CD4+_Activated/A002-C-016-D.bed

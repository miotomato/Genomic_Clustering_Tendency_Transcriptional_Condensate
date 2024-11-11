#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o processed_cell_Regions_slurms/slurm_immune_A014-C-043-D_CD8+.out

zcat GSE201336_RAW/GSM6058816_A014-C-043-D_08042020_fragments.tsv.gz | grep -f processed_cell_ID_per_sample_per_cellType/immune/CD8+/A014-C-043-D.txt > \
processed_cell_Regions_per_sample_per_cellType/immune/CD8+/A014-C-043-D.bed

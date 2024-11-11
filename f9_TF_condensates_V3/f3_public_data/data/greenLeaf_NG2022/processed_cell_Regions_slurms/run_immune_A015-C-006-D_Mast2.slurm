#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o processed_cell_Regions_slurms/slurm_immune_A015-C-006-D_Mast2.out

zcat GSE201336_RAW/GSM6058827_A015-C-006-D_20200811_fragments.tsv.gz | grep -f processed_cell_ID_per_sample_per_cellType/immune/Mast2/A015-C-006-D.txt > \
processed_cell_Regions_per_sample_per_cellType/immune/Mast2/A015-C-006-D.bed

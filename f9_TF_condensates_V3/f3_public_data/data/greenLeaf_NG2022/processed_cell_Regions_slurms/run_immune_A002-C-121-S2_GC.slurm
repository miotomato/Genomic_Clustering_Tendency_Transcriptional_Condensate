#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o processed_cell_Regions_slurms/slurm_immune_A002-C-121-S2_GC.out

zcat GSE201336_RAW/GSM6058800_A002-C-121-S2_20200310_fragments.tsv.gz | grep -f processed_cell_ID_per_sample_per_cellType/immune/GC/A002-C-121-S2.txt > \
processed_cell_Regions_per_sample_per_cellType/immune/GC/A002-C-121-S2.bed

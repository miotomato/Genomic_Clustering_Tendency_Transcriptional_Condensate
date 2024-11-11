#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o processed_cell_Regions_slurms/slurm_stromal_CRC-2-15564-D_Crypt_Fibroblasts_1.out

zcat GSE201336_RAW/GSM6058843_CRC-2-15564-D_20200917_fragments.tsv.gz | grep -f processed_cell_ID_per_sample_per_cellType/stromal/Crypt_Fibroblasts_1/CRC-2-15564-D.txt > \
processed_cell_Regions_per_sample_per_cellType/stromal/Crypt_Fibroblasts_1/CRC-2-15564-D.bed

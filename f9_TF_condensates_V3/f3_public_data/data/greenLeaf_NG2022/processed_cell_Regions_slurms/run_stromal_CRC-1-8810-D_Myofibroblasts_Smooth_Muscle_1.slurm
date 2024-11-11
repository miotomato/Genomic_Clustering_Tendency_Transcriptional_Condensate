#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o processed_cell_Regions_slurms/slurm_stromal_CRC-1-8810-D_Myofibroblasts_Smooth_Muscle_1.out

zcat GSE201336_RAW/GSM6058842_CRC-1-8810-D_20200917_fragments.tsv.gz | grep -f processed_cell_ID_per_sample_per_cellType/stromal/Myofibroblasts_Smooth_Muscle_1/CRC-1-8810-D.txt > \
processed_cell_Regions_per_sample_per_cellType/stromal/Myofibroblasts_Smooth_Muscle_1/CRC-1-8810-D.bed

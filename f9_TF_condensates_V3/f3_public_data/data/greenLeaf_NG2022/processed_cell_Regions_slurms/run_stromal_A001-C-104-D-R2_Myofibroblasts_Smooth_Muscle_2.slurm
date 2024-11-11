#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o processed_cell_Regions_slurms/slurm_stromal_A001-C-104-D-R2_Myofibroblasts_Smooth_Muscle_2.out

zcat GSE201336_RAW/GSM6058775_A001-C-104-D_20200214_fragments.tsv.gz \
GSE201336_RAW/GSM6058776_A001-C-104-D_20200811_fragments.tsv.gz | grep -f processed_cell_ID_per_sample_per_cellType/stromal/Myofibroblasts_Smooth_Muscle_2/A001-C-104-D-R2.txt > \
processed_cell_Regions_per_sample_per_cellType/stromal/Myofibroblasts_Smooth_Muscle_2/A001-C-104-D-R2.bed

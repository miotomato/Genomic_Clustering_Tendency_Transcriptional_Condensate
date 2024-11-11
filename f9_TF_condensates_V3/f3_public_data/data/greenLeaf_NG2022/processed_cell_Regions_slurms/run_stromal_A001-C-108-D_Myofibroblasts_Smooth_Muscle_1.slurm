#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o processed_cell_Regions_slurms/slurm_stromal_A001-C-108-D_Myofibroblasts_Smooth_Muscle_1.out

zcat GSE201336_RAW/GSM6058777_A001-C-108-D_20200811_fragments.tsv.gz | grep -f processed_cell_ID_per_sample_per_cellType/stromal/Myofibroblasts_Smooth_Muscle_1/A001-C-108-D.txt > \
processed_cell_Regions_per_sample_per_cellType/stromal/Myofibroblasts_Smooth_Muscle_1/A001-C-108-D.bed

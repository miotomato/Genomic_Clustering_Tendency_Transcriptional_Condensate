#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o processed_cell_Regions_slurms/slurm_immune_A002-C-010-D-R1_Mast2.out

zcat GSE201336_RAW/GSM6058787_A002-C-010-D_20200310_fragments.tsv.gz \
GSE201336_RAW/GSM6058788_A002-C-010-D_20200702_fragments.tsv.gz | grep -f processed_cell_ID_per_sample_per_cellType/immune/Mast2/A002-C-010-D-R1.txt > \
processed_cell_Regions_per_sample_per_cellType/immune/Mast2/A002-C-010-D-R1.bed

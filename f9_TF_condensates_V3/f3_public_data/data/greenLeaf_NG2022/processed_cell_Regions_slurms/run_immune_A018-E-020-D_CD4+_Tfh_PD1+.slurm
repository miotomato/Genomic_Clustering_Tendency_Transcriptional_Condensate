#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o processed_cell_Regions_slurms/slurm_immune_A018-E-020-D_CD4+_Tfh_PD1+.out

zcat GSE201336_RAW/GSM6058840_A018-E-020-D_20200817_fragments.tsv.gz | grep -f processed_cell_ID_per_sample_per_cellType/immune/CD4+_Tfh_PD1+/A018-E-020-D.txt > \
processed_cell_Regions_per_sample_per_cellType/immune/CD4+_Tfh_PD1+/A018-E-020-D.bed

#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o processed_cell_Regions_slurms/slurm_immune_B001-A-302-D_Exhausted_T_cells.out

zcat GSE201336_RAW/GSM6058851_B001-A-302-D_fragments.tsv.gz | grep -f processed_cell_ID_per_sample_per_cellType/immune/Exhausted_T_cells/B001-A-302-D.txt > \
processed_cell_Regions_per_sample_per_cellType/immune/Exhausted_T_cells/B001-A-302-D.bed

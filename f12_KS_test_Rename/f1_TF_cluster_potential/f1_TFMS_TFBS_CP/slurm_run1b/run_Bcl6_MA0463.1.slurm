#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 12:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_run1b/slurm_Bcl6_MA0463.1.out

time python run1b_TFMS_enrich_at_SE.py -i CP_TFMS_vs_random/_csv/Bcl6_MA0463.1_Exponential_Pvalue.csv


#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 12:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_run1b/slurm_Gata1_MA0035.3.out

time python run1b_TFMS_enrich_at_SE.py -i CP_TFMS_vs_random/_csv/Gata1_MA0035.3_Exponential_Pvalue.csv


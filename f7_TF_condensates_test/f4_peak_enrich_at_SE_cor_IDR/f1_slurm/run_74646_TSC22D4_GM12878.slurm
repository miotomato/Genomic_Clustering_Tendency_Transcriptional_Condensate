#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f1_slurm/slurm_74646_TSC22D4_GM12878.out

time python py1_TFBS_enrichAt_SE.py -i 74646 -f TSC22D4 -c GM12878


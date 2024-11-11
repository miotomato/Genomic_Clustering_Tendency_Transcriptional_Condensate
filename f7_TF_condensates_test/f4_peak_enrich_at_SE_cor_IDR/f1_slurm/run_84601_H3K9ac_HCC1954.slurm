#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f1_slurm/slurm_84601_H3K9ac_HCC1954.out

time python py1_TFBS_enrichAt_SE.py -i 84601 -f H3K9ac -c HCC1954


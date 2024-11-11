#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f1_slurm/slurm_81141_H3K18ac_HCT-116.out

time python py1_TFBS_enrichAt_SE.py -i 81141 -f H3K18ac -c HCT-116


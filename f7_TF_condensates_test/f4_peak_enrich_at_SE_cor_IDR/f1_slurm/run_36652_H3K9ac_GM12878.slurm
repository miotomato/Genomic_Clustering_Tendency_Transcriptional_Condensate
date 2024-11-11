#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f1_slurm/slurm_36652_H3K9ac_GM12878.out

time python py1_TFBS_enrichAt_SE.py -i 36652 -f H3K9ac -c GM12878


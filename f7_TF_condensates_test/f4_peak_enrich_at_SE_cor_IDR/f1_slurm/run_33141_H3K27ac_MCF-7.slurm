#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f1_slurm/slurm_33141_H3K27ac_MCF-7.out

time python py1_TFBS_enrichAt_SE.py -i 33141 -f H3K27ac -c MCF-7


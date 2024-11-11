#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f1_slurm/slurm_33176_H3K27ac_H2171.out

time python py1_TFBS_enrichAt_SE.py -i 33176 -f H3K27ac -c H2171


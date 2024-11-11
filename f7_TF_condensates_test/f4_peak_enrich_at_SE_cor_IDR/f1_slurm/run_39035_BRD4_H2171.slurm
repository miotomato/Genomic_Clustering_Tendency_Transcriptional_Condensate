#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f1_slurm/slurm_39035_BRD4_H2171.out

time python py1_TFBS_enrichAt_SE.py -i 39035 -f BRD4 -c H2171


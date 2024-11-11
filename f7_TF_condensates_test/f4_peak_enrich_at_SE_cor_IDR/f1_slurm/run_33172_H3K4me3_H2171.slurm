#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f1_slurm/slurm_33172_H3K4me3_H2171.out

time python py1_TFBS_enrichAt_SE.py -i 33172 -f H3K4me3 -c H2171


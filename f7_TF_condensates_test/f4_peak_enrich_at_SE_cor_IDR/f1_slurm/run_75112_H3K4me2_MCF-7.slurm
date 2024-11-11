#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f1_slurm/slurm_75112_H3K4me2_MCF-7.out

time python py1_TFBS_enrichAt_SE.py -i 75112 -f H3K4me2 -c MCF-7


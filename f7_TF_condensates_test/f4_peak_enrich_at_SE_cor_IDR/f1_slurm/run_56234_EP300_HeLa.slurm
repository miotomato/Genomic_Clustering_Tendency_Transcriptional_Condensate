#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f1_slurm/slurm_56234_EP300_HeLa.out

time python py1_TFBS_enrichAt_SE.py -i 56234 -f EP300 -c HeLa


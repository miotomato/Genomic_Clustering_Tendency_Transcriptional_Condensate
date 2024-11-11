#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f1_slurm/slurm_74153_NR3C1_MCF-7.out

time python py1_TFBS_enrichAt_SE.py -i 74153 -f NR3C1 -c MCF-7


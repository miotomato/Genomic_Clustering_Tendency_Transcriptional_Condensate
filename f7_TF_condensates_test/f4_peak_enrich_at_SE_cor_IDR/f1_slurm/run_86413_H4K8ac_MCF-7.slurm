#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f1_slurm/slurm_86413_H4K8ac_MCF-7.out

time python py1_TFBS_enrichAt_SE.py -i 86413 -f H4K8ac -c MCF-7


#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f1_slurm/slurm_62610_H2AFZ_MCF-7.out

time python py1_TFBS_enrichAt_SE.py -i 62610 -f H2AFZ -c MCF-7


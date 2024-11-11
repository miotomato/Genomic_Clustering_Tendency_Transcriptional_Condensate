#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f1_slurm/slurm_5786_H2B_HeLa.out

time python py1_TFBS_enrichAt_SE.py -i 5786 -f H2B -c HeLa


#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f1_slurm/slurm_54234_RUNX1_Jurkat.out

time python py1_TFBS_enrichAt_SE.py -i 54234 -f RUNX1 -c Jurkat


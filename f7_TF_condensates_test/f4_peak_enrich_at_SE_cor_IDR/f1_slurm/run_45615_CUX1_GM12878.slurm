#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f1_slurm/slurm_45615_CUX1_GM12878.out

time python py1_TFBS_enrichAt_SE.py -i 45615 -f CUX1 -c GM12878


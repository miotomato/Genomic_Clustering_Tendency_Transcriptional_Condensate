#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f1_slurm/slurm_64382_SUZ12_GM12878.out

time python py1_TFBS_enrichAt_SE.py -i 64382 -f SUZ12 -c GM12878


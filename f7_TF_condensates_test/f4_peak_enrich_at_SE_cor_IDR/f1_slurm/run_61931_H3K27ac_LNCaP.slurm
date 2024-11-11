#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f1_slurm/slurm_61931_H3K27ac_LNCaP.out

time python py1_TFBS_enrichAt_SE.py -i 61931 -f H3K27ac -c LNCaP


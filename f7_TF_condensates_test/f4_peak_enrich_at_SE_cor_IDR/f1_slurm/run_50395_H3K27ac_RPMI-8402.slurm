#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f1_slurm/slurm_50395_H3K27ac_RPMI-8402.out

time python py1_TFBS_enrichAt_SE.py -i 50395 -f H3K27ac -c RPMI-8402


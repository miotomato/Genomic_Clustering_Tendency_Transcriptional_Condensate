#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f1_slurm/slurm_64081_NFE2L1_K562.out

time python py1_TFBS_enrichAt_SE.py -i 64081 -f NFE2L1 -c K562


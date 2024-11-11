#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f1_slurm/slurm_48994_AFF4_HeLa.out

time python py1_TFBS_enrichAt_SE.py -i 48994 -f AFF4 -c HeLa


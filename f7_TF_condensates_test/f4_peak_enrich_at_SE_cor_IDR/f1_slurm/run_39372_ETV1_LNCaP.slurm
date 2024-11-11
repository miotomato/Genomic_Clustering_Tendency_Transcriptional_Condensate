#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f1_slurm/slurm_39372_ETV1_LNCaP.out

time python py1_TFBS_enrichAt_SE.py -i 39372 -f ETV1 -c LNCaP


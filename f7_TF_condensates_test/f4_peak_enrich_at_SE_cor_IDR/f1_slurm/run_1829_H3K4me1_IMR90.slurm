#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f1_slurm/slurm_1829_H3K4me1_IMR90.out

time python py1_TFBS_enrichAt_SE.py -i 1829 -f H3K4me1 -c IMR90


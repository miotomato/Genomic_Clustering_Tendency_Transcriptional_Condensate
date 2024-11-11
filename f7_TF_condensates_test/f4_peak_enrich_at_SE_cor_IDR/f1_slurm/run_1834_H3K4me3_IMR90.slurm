#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f1_slurm/slurm_1834_H3K4me3_IMR90.out

time python py1_TFBS_enrichAt_SE.py -i 1834 -f H3K4me3 -c IMR90


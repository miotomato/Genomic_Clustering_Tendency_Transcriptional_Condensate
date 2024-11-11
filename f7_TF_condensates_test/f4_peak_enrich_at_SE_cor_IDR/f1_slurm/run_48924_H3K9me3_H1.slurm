#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f1_slurm/slurm_48924_H3K9me3_H1.out

time python py1_TFBS_enrichAt_SE.py -i 48924 -f H3K9me3 -c H1


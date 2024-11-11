#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f1_slurm/slurm_6421_H3K36me3_H1.out

time python py1_TFBS_enrichAt_SE.py -i 6421 -f H3K36me3 -c H1


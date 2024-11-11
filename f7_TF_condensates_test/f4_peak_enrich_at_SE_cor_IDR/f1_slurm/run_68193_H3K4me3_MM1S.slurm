#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f1_slurm/slurm_68193_H3K4me3_MM1S.out

time python py1_TFBS_enrichAt_SE.py -i 68193 -f H3K4me3 -c MM1S


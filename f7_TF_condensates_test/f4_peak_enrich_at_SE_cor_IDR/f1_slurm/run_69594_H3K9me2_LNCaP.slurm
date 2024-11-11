#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f1_slurm/slurm_69594_H3K9me2_LNCaP.out

time python py1_TFBS_enrichAt_SE.py -i 69594 -f H3K9me2 -c LNCaP


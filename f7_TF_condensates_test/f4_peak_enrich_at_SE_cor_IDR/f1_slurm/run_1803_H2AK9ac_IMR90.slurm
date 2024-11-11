#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f1_slurm/slurm_1803_H2AK9ac_IMR90.out

time python py1_TFBS_enrichAt_SE.py -i 1803 -f H2AK9ac -c IMR90


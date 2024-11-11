#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f1_slurm/slurm_9196_GATA3_RPMI-8402.out

time python py1_TFBS_enrichAt_SE.py -i 9196 -f GATA3 -c RPMI-8402


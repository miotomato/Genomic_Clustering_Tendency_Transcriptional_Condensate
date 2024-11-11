#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f1_slurm/slurm_9194_TCF7L1_RPMI-8402.out

time python py1_TFBS_enrichAt_SE.py -i 9194 -f TCF7L1 -c RPMI-8402


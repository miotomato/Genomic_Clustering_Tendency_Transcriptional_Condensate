#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f1_slurm/slurm_43061_TCF7L2_LNCaP.out

time python py1_TFBS_enrichAt_SE.py -i 43061 -f TCF7L2 -c LNCaP


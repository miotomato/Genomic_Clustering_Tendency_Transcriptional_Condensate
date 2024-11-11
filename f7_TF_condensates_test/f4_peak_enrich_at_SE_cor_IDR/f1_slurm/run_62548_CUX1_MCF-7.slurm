#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f1_slurm/slurm_62548_CUX1_MCF-7.out

time python py1_TFBS_enrichAt_SE.py -i 62548 -f CUX1 -c MCF-7


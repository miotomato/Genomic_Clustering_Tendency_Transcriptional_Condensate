#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=6000
#SBATCH -t 12:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o run_py1_slurms/slurm_Astrocytes.out

time python py1_cluster_cor_SE_fisher_exact.py -i Astrocytes
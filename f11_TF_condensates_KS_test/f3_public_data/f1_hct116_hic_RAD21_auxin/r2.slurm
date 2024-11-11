#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=80000
#SBATCH -t 48:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o out_r2.log

# Run program


time python r2_DCI_overlap_TFBS_clustered.py


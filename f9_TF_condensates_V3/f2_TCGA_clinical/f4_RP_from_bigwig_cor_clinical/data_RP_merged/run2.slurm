#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=60000
#SBATCH -t 48:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o out2.log

# Run program
time python r2_avg_PR_per_sample.py



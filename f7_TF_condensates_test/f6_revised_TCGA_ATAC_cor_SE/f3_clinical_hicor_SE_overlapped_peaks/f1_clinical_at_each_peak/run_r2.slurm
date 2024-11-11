#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 48:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o out_run_r2.log

# Run program
time python r2_caseID_openness_each_SE_vs_clinical.py


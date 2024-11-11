#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=80000
#SBATCH -t 20:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o out_run2.log

time python run2_CI_intersect_TFBS.py


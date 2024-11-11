#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 122:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o out_r3.log


time python r3_clinical_per_region.py 
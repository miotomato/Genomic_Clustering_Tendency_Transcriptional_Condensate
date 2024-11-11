#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=80000
#SBATCH -t 40:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o out_run1.log

time python run1_get_chromatin_interaction.py


#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=80000
#SBATCH -t 10:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o out_run3.log

time python run3_TFBS_CI_plot.py


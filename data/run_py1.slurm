#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=300000
#SBATCH -t 12:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o out_py1.log2

# Run program

# time python py1_cistrome_merged_peaks.py
time python py1_cistrome_merged_peaks_append.py

#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_run2/slurm_SK-N-SH_YY1_46342.out

module load gcc/7.1.0
module load bedtools/2.26.0

time python run2_TFBS_vs_TFMS_sample.py -i 46342


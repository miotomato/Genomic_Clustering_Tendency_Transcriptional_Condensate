#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 12:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_run3/slurm_50989_CTCF_H1.out

module load gcc/7.1.0
module load bedtools/2.26.0

time python run3_TFBS_with_motif_vs_TFMS_sample.py -i 50989


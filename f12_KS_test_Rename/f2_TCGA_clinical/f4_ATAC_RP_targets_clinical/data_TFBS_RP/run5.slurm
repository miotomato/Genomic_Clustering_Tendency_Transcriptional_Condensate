#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=60000
#SBATCH -t 48:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o out5.log

# Run program
time python r5_gene_RP_by_ATAC_overlap_TFBS_per_patient.py



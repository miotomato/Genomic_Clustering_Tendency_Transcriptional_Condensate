#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=60000
#SBATCH -t 48:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o out4.log

# Run program
time python r4_clinical_per_gene_by_TFBS_RP.py



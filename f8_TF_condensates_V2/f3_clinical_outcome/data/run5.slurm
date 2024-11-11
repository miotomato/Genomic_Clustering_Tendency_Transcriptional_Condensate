#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 48:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o out_run5.slurm.log

# Run program

time python r5_gene_RP_by_ATAC_overlap_SE_overlap_TFMS_per_patient.py > out_run5.log



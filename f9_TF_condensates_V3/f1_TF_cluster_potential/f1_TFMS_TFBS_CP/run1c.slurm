#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 122:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o out_run1c.log

time python run1c_TFMS_enrich_at_each_CT_SE.py




#!/bin/bash
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --mem=60gb
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_step2_manu.log


export PATH="/nv/vol190/zanglab/zw5j/env/anaconda2/bin:$PATH"
export PATH="/nv/vol190/zanglab/zw5j/env/hicpro/installation/HiC-Pro_2.10.0/bin/utils/:$PATH"
module load gcc/7.1.0


# Run ICE Normalization ...
/nv/vol190/zanglab/zw5j/env/hicpro/installation/HiC-Pro_2.10.0/scripts/ice_norm.sh -c /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/f11_TF_condensates_KS_test/f3_public_data/data_processed/public_hic_wt/config_mboi.txt >> hicpro.log 



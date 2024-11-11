#!/bin/bash
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --mem=60gb
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o slurm_step2_manu.log

# Merge multiple files from the same sample ...
# /nv/vol190/zanglab/zw5j/env/hicpro/installation/HiC-Pro_2.10.0/scripts/merge_valid_interactions.sh -c /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/f11_TF_condensates_KS_test/f3_public_data/data_processed/GSE104333_hct116_hic/confit_test_latest_Mboi.txt >> hicpro.log


export PATH="/nv/vol190/zanglab/zw5j/env/anaconda2/bin:$PATH"
export PATH="/nv/vol190/zanglab/zw5j/env/hicpro/installation/HiC-Pro_2.10.0/bin/utils/:$PATH"
module load gcc/7.1.0

# Merge stat files per sample ...
/nv/vol190/zanglab/zw5j/env/hicpro/installation/HiC-Pro_2.10.0/scripts/merge_stats.sh -c /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/f11_TF_condensates_KS_test/f3_public_data/data_processed/GSE104333_hct116_hic/confit_test_latest_Mboi.txt >> hicpro.log

# Generate binned matrix files ...
/nv/vol190/zanglab/zw5j/env/hicpro/installation/HiC-Pro_2.10.0/scripts/build_raw_maps.sh -c /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/f11_TF_condensates_KS_test/f3_public_data/data_processed/GSE104333_hct116_hic/confit_test_latest_Mboi.txt

# Run quality checks for all samples ...
/nv/vol190/zanglab/zw5j/env/hicpro/installation/HiC-Pro_2.10.0/scripts/make_plots.sh -c /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/f11_TF_condensates_KS_test/f3_public_data/data_processed/GSE104333_hct116_hic/confit_test_latest_Mboi.txt -p "all" >> hicpro.log

# Run ICE Normalization ...
/nv/vol190/zanglab/zw5j/env/hicpro/installation/HiC-Pro_2.10.0/scripts/ice_norm.sh -c /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/f11_TF_condensates_KS_test/f3_public_data/data_processed/GSE104333_hct116_hic/confit_test_latest_Mboi.txt >> hicpro.log 

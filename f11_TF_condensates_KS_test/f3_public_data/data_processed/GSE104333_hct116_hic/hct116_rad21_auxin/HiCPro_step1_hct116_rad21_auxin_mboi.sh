#!/bin/bash
#SBATCH -n 8
#SBATCH -t 72:00:00
#SBATCH --mem-per-cpu=20gb
#SBATCH -p parallel
#SBATCH -A cphg_cz3d

#SBATCH --mail-user=zhenjia@virginia.edu
#SBATCH --mail-type=end
#SBATCH --job-name=HiCpro_s1_hct116_rad21_auxin_mboi
#SBATCH --export=ALL
#SBATCH --array=1-125

FASTQFILE=$SLURM_SUBMIT_DIR/inputfiles_hct116_rad21_auxin_mboi.txt; export FASTQFILE
make --file /nv/vol190/zanglab/zw5j/env/hicpro/installation/HiC-Pro_2.10.0/scripts/Makefile CONFIG_FILE=/nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/f11_TF_condensates_KS_test/f3_public_data/data_processed/GSE104333_hct116_hic/confit_test_latest_Mboi.txt CONFIG_SYS=/nv/vol190/zanglab/zw5j/env/hicpro/installation/HiC-Pro_2.10.0/config-system.txt all_sub 2>&1

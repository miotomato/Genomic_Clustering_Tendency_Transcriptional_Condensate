#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=80000
#SBATCH -t 48:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o out_bart3d.log

# Run program


outdir=f1_bart3d
# mkdir ${outdir}

project_dir=/nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/f11_TF_condensates_KS_test
hicpro_results=${project_dir}//f3_public_data/data_processed/GSE104333_hct116_hic/hct116_rad21_auxin/hic_results/matrix

treat_arr=(HCT116_RAD21_6hr_auxin_rep1 HCT116_RAD21_6hr_auxin_rep2)
ctrl_arr=(HCT116_no_auxin_rep1 HCT116_no_auxin_rep2)


for genomicDistance in 100000 200000 500000
do
    time bart3d -f hicpro -s hg38 \
           -t ${hicpro_results}/${treat_arr[0]}/raw/5000/${treat_arr[0]}_5000.matrix \
           -c ${hicpro_results}/${ctrl_arr[0]}/raw/5000/${ctrl_arr[0]}_5000.matrix \
           --bedFileHicpro ${hicpro_results}/${treat_arr[0]}/raw/5000/${treat_arr[0]}_5000_abs.bed \
           --genomicDistance ${genomicDistance} \
           --outdir f1_bart3d \
           --outFileName RAD21_6hr_auxin_over_NT_rep1_dis$(($genomicDistance/1000))k 


    time bart3d -f hicpro -s hg38 \
           -t ${hicpro_results}/${treat_arr[1]}/raw/5000/${treat_arr[1]}_5000.matrix \
           -c ${hicpro_results}/${ctrl_arr[1]}/raw/5000/${ctrl_arr[1]}_5000.matrix \
           --bedFileHicpro ${hicpro_results}/${treat_arr[0]}/raw/5000/${treat_arr[0]}_5000_abs.bed \
           --genomicDistance ${genomicDistance} \
           --outdir f1_bart3d \
           --outFileName RAD21_6hr_auxin_over_NT_rep2_dis$(($genomicDistance/1000))k 


    time bart3d -f hicpro -s hg38 \
           -t ${hicpro_results}/${treat_arr[0]}/raw/5000/${treat_arr[0]}_5000.matrix,${hicpro_results}/${treat_arr[1]}/raw/5000/${treat_arr[1]}_5000.matrix \
           -c ${hicpro_results}/${ctrl_arr[0]}/raw/5000/${ctrl_arr[0]}_5000.matrix,${hicpro_results}/${ctrl_arr[1]}/raw/5000/${ctrl_arr[1]}_5000.matrix \
           --bedFileHicpro ${hicpro_results}/${treat_arr[0]}/raw/5000/${treat_arr[0]}_5000_abs.bed \
           --genomicDistance ${genomicDistance} \
           --outdir f1_bart3d \
           --outFileName RAD21_6hr_auxin_over_NT_all_reps_dis$(($genomicDistance/1000))k 
done



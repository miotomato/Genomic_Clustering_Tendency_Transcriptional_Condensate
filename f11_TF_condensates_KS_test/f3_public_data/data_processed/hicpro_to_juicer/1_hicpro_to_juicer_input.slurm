#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=80000
#SBATCH -t 48:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o hicpro_to_juicer.out

# source /nv/vol190/zanglab/zw5j/env/hicpro/run_hic_pro_doc/run_hic_pro_prepare.sh
# export PATH="/nv/vol190/zanglab/zw5j/env/hicpro/installation/HiC-Pro_2.10.0/bin/utils/:$PATH"
export PATH="/nv/vol190/zanglab/zw5j/env/anaconda2/bin:$PATH"
export PATH="/nv/vol190/zanglab/zw5j/env/hicpro/installation/HiC-Pro_2.10.0/bin/utils/:$PATH"

outdir="juicer_hic"
mkdir $outdir

genome_size="/nv/vol190/zanglab/zw5j/data//Genome/ucsc/hg38/hg38.chrom.sizes"
juicer_tool="/nv/vol190/zanglab/zw5j/env/HiC_tools/Juicebox/juicer_tools.1.8.9_jcuda.0.8.jar"

res_mboi="/nv/vol190/zanglab/zw5j/env/hicpro/annotation/hg38_mboi.bed"
res_dpnii="/nv/vol190/zanglab/zw5j/env/hicpro/annotation/hg38_dpnii.bed"
res_bglii="/nv/vol190/zanglab/zw5j/env/hicpro/annotation/hg38_bglii.bed"
res_hindiii="/nv/vol190/zanglab/zw5j/env/hicpro/annotation/hg38_hindiii.bed"


filedir="/nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/f11_TF_condensates_KS_test/f3_public_data/data_processed/public_hic_wt"
hicpro2juicebox.sh -i ${filedir}/bglii/hic_results/data/LNCaP/LNCaP_allValidPairs -g ${genome_size} -r ${res_bglii} -o $outdir -j ${juicer_tool}
hicpro2juicebox.sh -i ${filedir}/dpnii/hic_results/data/H1/H1_allValidPairs -g ${genome_size} -r ${res_dpnii} -o $outdir -j ${juicer_tool}
hicpro2juicebox.sh -i ${filedir}/hindiii/hic_results/data/HepG2/HepG2_allValidPairs -g ${genome_size} -r ${res_hindiii} -o $outdir -j ${juicer_tool}
hicpro2juicebox.sh -i ${filedir}/hindiii/hic_results/data/MCF7/MCF7_allValidPairs -g ${genome_size} -r ${res_hindiii} -o $outdir -j ${juicer_tool}
hicpro2juicebox.sh -i ${filedir}/hindiii/hic_results/data/Panc1/Panc1_allValidPairs -g ${genome_size} -r ${res_hindiii} -o $outdir -j ${juicer_tool}
hicpro2juicebox.sh -i ${filedir}/mboi/hic_results/data/GM12878/GM12878_allValidPairs -g ${genome_size} -r ${res_mboi} -o $outdir -j ${juicer_tool}
hicpro2juicebox.sh -i ${filedir}/mboi/hic_results/data/HeLa/HeLa_allValidPairs -g ${genome_size} -r ${res_mboi} -o $outdir -j ${juicer_tool}
hicpro2juicebox.sh -i ${filedir}/mboi/hic_results/data/IMR90/IMR90_allValidPairs -g ${genome_size} -r ${res_mboi} -o $outdir -j ${juicer_tool}
hicpro2juicebox.sh -i ${filedir}/mboi/hic_results/data/K562/K562_allValidPairs -g ${genome_size} -r ${res_mboi} -o $outdir -j ${juicer_tool}

allValidPairs_file="/nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/f11_TF_condensates_KS_test/f3_public_data/data_processed/GSE104333_hct116_hic/hct116_rad21_auxin/hic_results/data/HCT116_no_auxin_rep1/HCT116_no_auxin_rep1_allValidPairs"
hicpro2juicebox.sh -i ${allValidPairs_file} -g ${genome_size} -r ${res_mboi} -o $outdir -j ${juicer_tool}


# jurkat_matrix_dir=/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/11_HiC_analysis/f1_preprocess/HiC_Pro/hicPro_raw_fq/panos_hic/hic_results/matrix 
# jurkat_hic_file=/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/9_get_CTCF_signals_plus_appended_data/f8_hic_TADs/f0_hic_data_hg38/view_hic/Jurkat_inter_30.hic 

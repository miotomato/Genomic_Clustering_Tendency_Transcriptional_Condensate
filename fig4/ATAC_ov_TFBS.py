import os
import numpy as np
import pandas as pd

# == Output directory
outdir = 'f1_ATAC_overlap_TFBS_caseID'
os.makedirs(outdir, exist_ok=True)

# == Project directory
project_dir = '/standard/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR'

# == Read matched names between TCGA ATAC and SE data
name_match_file = f"{project_dir}/f7_TF_condensates_test/f6_revised_TCGA_ATAC_cor_SE/data/TCGA/TCGA-ATAC_SE_cancerType_match.xlsx"
name_match = pd.read_excel(name_match_file, index_col=0)
name_match = name_match.dropna()

# == Path to the Python script for overlap calculation
pyfile = "find_overlap_keep_info_NOT_sep_strand_revised.py"

# == TCGA ATAC-seq normalized counts file
tcga_file = f"{project_dir}/f7_TF_condensates_test/f6_revised_TCGA_ATAC_cor_SE/data/TCGA/mynorm_TCGA-ATAC_PanCan_Log2_QuantileNorm_Counts_plus5.caseID.avg.txt"

# == Directory containing clustered TFBS data
c_tfbs_dir = f"{project_dir}/f15_revision/f1_TF_cluster_potential/f3_clustered_TFBS/f5_atac_overlap_coBinding_TFBS"

# == Get top 3 transcription factors (TFs) based on TFBS CP rank and z-score
selected_factors = {}
tfbs_cp_dir = f"{project_dir}/f15_revision/f1_TF_cluster_potential/f2_cor_CP_SE_AICAP/f9_per_CT_TFBS_CP_cor_zscore_CP_with_motif_SE/TFBS_CP"

for ct in ['MCF-7', 'HCT-116']:
    # Read the TFBS CP ranking data
    df = pd.read_csv(f"{tfbs_cp_dir}/_CP_TFBS_nonBlackList_vs_TFMS_{ct}.csv", index_col=0)
    # Select the top 3 TFs based on TFBS CP rank and average rank (z-score)
    selected_factors[f"{ct} top_TFBSCP"] = df['TFBS CP rank'].sort_values().iloc[:3].index
    selected_factors[f"{ct} top_zscored_TFBSCP"] = df['avg rank'].sort_values().iloc[:3].index

# == Process overlap for each cancer type and TFBS category
for cancertype in ['BRCA', 'COAD']:
    # Get the cell type (SE) associated with the cancer type
    ct = name_match.loc[cancertype].SE

    # Iterate over treatment flags and factor types
    for treat_flag in ['percentile_T_ExtendMerge']:
        for factorType in ['top_zscored_TFBSCP']:
            # Create output subdirectory for the current cell type and factor type
            subdir = f"{ct}_{factorType}"
            os.makedirs(f"{outdir}/{subdir}", exist_ok=True)

            # Get the selected TFs for the current cell type and factor type
            factors = selected_factors[f"{ct} {factorType}"]

            # Define the x-tick labels for different combinations of TFs
            xticklabels = [
                'Union',
                '-'.join([factors[0]]),
                '-'.join([factors[1]]),
                '-'.join([factors[2]]),
                '-'.join(factors[:2]),
                '-'.join(factors[:3]),
            ]

            # Iterate over each TF combination and perform overlap analysis
            for label in xticklabels:
                tfbs_file = f"{c_tfbs_dir}/{subdir}/{treat_flag}_{ct}_{label}.bed"
                overlap_file = f"{outdir}/{subdir}/ATAC_overlap_{treat_flag}_{ct}_{label}.bed"

                # Construct and execute the command for overlap analysis
                command = f"python {pyfile} -a {tcga_file} -b {tfbs_file} -s hg38 -p {overlap_file} -e2 0"
                print(command)
                os.system(command)

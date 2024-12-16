import os
import numpy as np
import pandas as pd

# == Input and Output Directories
indir = 'ATAC_overlap_TFBS_per_patient_sig'
outdir = 'ATAC_overlap_TFBS_per_patient_sig_RP'
os.makedirs(outdir, exist_ok=True)

# == Load Filtered TCGA ATAC Clustered Samples
filtered_file = "data/TCGA/TCGA-ATAC_clustered_samples.xlsx"
filtered_df = pd.read_excel(filtered_file, index_col=0)

# == Read Matched Names Between TCGA ATAC and SE Data
name_match_file = "data/TCGA-ATAC_SE_cancerType_match.xlsx"
name_match = pd.read_excel(name_match_file, index_col=0)
name_match = name_match.dropna()

# == Get Top 3 Factors Based on TFBS CP Rank and Z-score
selected_factors = {}
tfbs_cp_dir = "TFBS_CP"

for ct in ['MCF-7', 'HCT-116']:
    df = pd.read_csv(f"{tfbs_cp_dir}/_CP_TFBS_nonBlackList_vs_TFMS_{ct}.csv", index_col=0)
    selected_factors[f"{ct} top_TFBSCP"] = df['TFBS CP rank'].sort_values().iloc[:3].index
    selected_factors[f"{ct} top_zscored_TFBSCP"] = df['avg rank'].sort_values().iloc[:3].index

# == Process Overlap for Each Cancer Type and TFBS Category
for cancertype in ['BRCA', 'COAD']:
    # Get the cell type (SE) associated with the cancer type
    ct = name_match.loc[cancertype].SE
    # Get filtered patient IDs for the current cancer type
    filtered_ids = filtered_df[filtered_df.cohort == cancertype].case_id

    # Iterate over treatment flags and factor types
    for treat_flag in ['percentile_T_ExtendMerge']:
        for factorType in ['top_zscored_TFBSCP']:
            # Create output subdirectory for the current cell type and factor type
            subdir = f"{ct}_{factorType}"
            os.makedirs(f"{outdir}/{subdir}", exist_ok=True)

            # Get the selected TFs for the current cell type and factor type
            factors = selected_factors[f"{ct} {factorType}"]

            # Define x-tick labels for different combinations of TFs
            xticklabels = [
                'ALL',
                'Union',
                '-'.join([factors[0]]),
                '-'.join([factors[1]]),
                '-'.join([factors[2]]),
                '-'.join(factors[:2]),
                '-'.join(factors[:3]),
            ]

            # Iterate over each TF combination and patient ID
            for label in xticklabels:
                for patient in filtered_ids:
                    # Define the input BED file and output file paths
                    bed_file = f"{indir}/{subdir}/ATAC_overlap_{treat_flag}_{ct}_{label}_{patient}.bed"
                    out_file = f"{outdir}/{subdir}/ATAC_overlap_{treat_flag}_{ct}_{label}_{patient}_RP.tsv"

                    # Check if the input BED file exists before executing the command
                    if os.path.isfile(bed_file):
                        commandLine = f"python2 get-regulatory-potential-on-genes_peak_level.py -b {bed_file} -s hg38 -g hg38_unique_geneSymbol.ucsc -o {out_file}"
                        os.system(commandLine)
                        print(commandLine)

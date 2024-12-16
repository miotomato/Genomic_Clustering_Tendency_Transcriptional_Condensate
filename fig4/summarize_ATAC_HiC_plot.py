
import os
import sys
import argparse
import glob
import re
import bisect
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import gridspec
import scipy.optimize

# Set default plot styles
sns.set(font_scale=1)
sns.set_style("whitegrid", {'axes.grid': False})
matplotlib.rcParams['font.size'] = 12
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
sns.set_style("ticks")

# Define project directory (modify this path according to your project location)
project_dir = 'f12_KS_test_Rename'

# Define paths to different datasets
diff_atac_dir = '{}/f2_TCGA_clinical/f1_diff_ATAC/f9b_diff_ATAC_overlap_TFBS_clustered_data_figs'.format(project_dir)
atac_rp_dir = '{}/f2_TCGA_clinical/f2_RP_from_bigwig/f4_avg_RP_per_sample_across_patients_figs'.format(project_dir)
hic_dir = '{}/f3_public_data/f2_TFBS_CI/f3_CI_figs'.format(project_dir)
rank_dir = '../../f1_TF_cluster_potential/f2_cor_CP_SE_AICAP/f9_per_CT_TFBS_CP_cor_zscore_CP_with_motif_SE/TFBS_CP/'

# Load TFBS data
tfbs_file = '{}/data_merged_SE_overlapped.csv'.format(project_dir)
tfbs_df = pd.read_csv(tfbs_file, index_col=0)

# Load name match data for SE and cancer types
name_match = pd.read_excel('{}/TCGA/TCGA-ATAC_SE_cancerType_match.xlsx'.format(project_dir), index_col=0)
name_match = name_match.dropna()

# Output directory
outdir = 'f2_ATAC_HiC_compr_TFBS_heatmap'
os.makedirs(outdir, exist_ok=True)

# Define transcription factors to exclude
removed_tfs = ['HSF1', 'T', 'NR2C2']

# List of cancer types to analyze
cancertypes = ['BRCA', 'CESC', 'COAD', 'LIHC', 'PRAD']

# Treatment flags for different processing modes
treat_flags = ['percentile_T', 'percentile_T_ExtendMerge']

# Genomic distances in kilobases
genomic_dis_kbs = [20, 50, 100, 200, 500]

# Example: Restricting analysis to a subset of genomic distances
genomic_dis_kbs = [50]

# Main analysis loop
for genomic_dis in genomic_dis_kbs:
    for treat_flag in treat_flags[1:]:  # Loop through specified treatment flags
        for cancertype in cancertypes:
            # Get SE identifier for the current cancer type
            ct = name_match.loc[cancertype, 'SE']

            # Load differential ATAC-seq data
            atac_diff_file = '{}/{}_by_{}.csv'.format(diff_atac_dir, cancertype, treat_flag)
            atac_diff = pd.read_csv(atac_diff_file, index_col=0)

            # Load ATAC-seq RP data
            atac_rp_file = '{}/{}_halflife_10000_by_{}.merge.csv'.format(atac_rp_dir, cancertype, treat_flag)
            atac_rp = pd.read_csv(atac_rp_file, index_col=0)

            # Load Hi-C data
            hic_file = '{}/{}_{}_CI_{}KB.csv'.format(hic_dir, ct, treat_flag, genomic_dis)
            hic_df = pd.read_csv(hic_file, index_col=0)

            # Load rank data for TFBS
            rank_file = '{}/_CP_TFBS_nonBlackList_vs_TFMS_{}.csv'.format(rank_dir, ct)
            rank_df = pd.read_csv(rank_file, index_col=0)

            # Filter transcription factors that are not in the exclusion list
            kept_factors = [i for i in rank_df.index if '{} {}'.format(ct, i) in tfbs_df.index]

            # Perform analysis for each kept factor
            for kept_factor in kept_factors:
                tfbs_index = '{} {}'.format(ct, kept_factor)

                # Get TFBS data for the current factor
                se_overlap_count = tfbs_df.loc[tfbs_index, '# {} on SE'.format(treat_flag)]

                # Visualization: Plot the TFBS overlap count
                plt.figure(figsize=(8, 6))
                sns.barplot(x=[kept_factor], y=[se_overlap_count])
                plt.title('TFBS Overlap Count for {} in {}'.format(kept_factor, cancertype))
                plt.ylabel('Overlap Count')
                plt.xlabel('Transcription Factor')
                plt.tight_layout()

                # Save the plot to the output directory
                plot_filename = '{}/{}_{}_{}_{}KB.png'.format(outdir, cancertype, kept_factor, treat_flag, genomic_dis)
                plt.savefig(plot_filename)
                plt.close()

                print('Plot saved:', plot_filename)

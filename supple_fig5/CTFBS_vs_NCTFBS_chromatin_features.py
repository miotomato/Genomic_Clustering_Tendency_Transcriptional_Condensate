import os
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import gridspec
import matplotlib

# Set plot parameters
matplotlib.rcParams['font.size'] = 12
sns.set(font_scale=1)
sns.set_style("whitegrid", {'axes.grid': False})
sns.set_style("ticks")
matplotlib.rcParams["font.sans-serif"] = ["Arial"]

# Paths
project_dir = 'f12_KS_test_Rename'
diff_atac_dir = f'{project_dir}/f2_TCGA_clinical/f1_diff_ATAC/f9b_diff_ATAC_overlap_TFBS_clustered_data_figs'
atac_rp_dir = f'{project_dir}/f2_TCGA_clinical/f2_RP_from_bigwig/f4_avg_RP_per_sample_across_patients_figs'
hic_dir = f'{project_dir}/f3_public_data/f2_TFBS_CI/f3_CI_figs'
tfbs_file = f'{project_dir}/f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged/data_merged_SE_overlapped.csv'

# Read data
tfbs_df = pd.read_csv(tfbs_file, index_col=0)
name_match = pd.read_excel(f'{project_dir}/../f9_TF_condensates_V3/data/TCGA/TCGA-ATAC_SE_cancerType_match.xlsx', index_col=0).dropna()

# Output directory
outdir = 'f2_ATAC_HiC_compr_TFBS_heatmap'
os.makedirs(outdir, exist_ok=True)

removed_tfs = ['HSF1', 'T', 'NR2C2']
cancertypes = ['BRCA', 'CESC', 'COAD', 'LIHC', 'PRAD']
treat_flags = ['percentile_T', 'percentile_T_ExtendMerge']
genomic_dis_kbs = [50]

rank_dir = 'TFBS_CP'

for genomic_dis in genomic_dis_kbs:
    for treat_flag in treat_flags[1:]:
        for cancertype in cancertypes:
            ct = name_match.loc[cancertype, 'SE']
            atac_diff = pd.read_csv(f'{diff_atac_dir}/{cancertype}_by_{treat_flag}.csv', index_col=0)
            atac_rp = pd.read_csv(f'{atac_rp_dir}/{cancertype}_halflife_10000_by_{treat_flag}.merge.csv', index_col=0)
            hic_df = pd.read_csv(f'{hic_dir}/_{ct}_{treat_flag}_CI_{genomic_dis}KB.csv', index_col=0)
            rank_df = pd.read_csv(f'{rank_dir}/_CP_TFBS_nonBlackList_vs_TFMS_{ct}.csv', index_col=0)

            kept_factors = [i for i in rank_df.index if f'{ct} {i}' in tfbs_df.index]

            out_df = pd.DataFrame()
            for kept_factor in kept_factors:
                tfbs_index = f'{ct} {kept_factor}'
                a = tfbs_df.loc[tfbs_index, f'# {treat_flag} on SE']
                b = tfbs_df.loc[tfbs_index, f'# {treat_flag}']
                c = tfbs_df.loc[tfbs_index, '# percentile_C on SE']
                d = tfbs_df.loc[tfbs_index, '# percentile_C']
                s, p = stats.fisher_exact([[a, b], [c, d]])
                out_df.loc[kept_factor, 'SE enrich statistics'] = s
                out_df.loc[kept_factor, 'SE enrich pvalue'] = p

            # Replace zeros in p-values for plotting
            out_df.replace(0, 1e-299, inplace=True)

            # Save results
            out_df.to_csv(f'{outdir}/data_{cancertype}_{treat_flag}_{genomic_dis}KB.csv')

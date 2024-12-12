import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import warnings

# == Configure Plotting Parameters
warnings.filterwarnings("ignore")
matplotlib.rcParams['font.size'] = 13
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"

sns.set(font_scale=1)
sns.set_style("whitegrid", {'axes.grid': False})
sns.set_style("ticks")

# == Constants
p_thre = 0.05
flag = 'ATAC-seq peaks'

# == Functions

def add_logrank_info(out_df, df, motif_name, flag):
    """Adds log-rank test information to the output DataFrame."""
    df_p = df[df['log rank p'] < p_thre]
    df_pt = df[(df['treat time'] < df['ctrl time']) & (df['log rank p'] < p_thre)]

    out_df.loc[motif_name, f'{flag} total'] = df.shape[0]
    out_df.loc[motif_name, f'{flag} #P<{p_thre}'] = df_p.shape[0]
    out_df.loc[motif_name, f'{flag} %P<{p_thre}'] = np.round(df_p.shape[0] / df.shape[0], 4)
    out_df.loc[motif_name, f'{flag} #TreatTime<CtrlTime&P<{p_thre}'] = df_pt.shape[0]
    out_df.loc[motif_name, f'{flag} %TreatTime<CtrlTime&P<{p_thre}'] = np.round(df_pt.shape[0] / df.shape[0], 4)
    return out_df


def stack_bar(df, xticklabels, figname, cancertype, flag):
    """Generates and saves a stacked bar plot for log-rank test results."""
    plt.figure(figsize=(3, 2.6))
    positions = np.arange(df.shape[0])

    # Get total and significant counts for 'All'
    t_all = df.loc['All', f'{flag} total']
    p_all = df.loc['All', f'{flag} #P<{p_thre}']

    # Plot bars and add significance markers
    for position in positions:
        index = df.index[position]
        total = df.loc[index, f'{flag} total']
        significant = df.loc[index, f'{flag} #P<{p_thre}']

        # Fisher's exact test for significance
        _, p_value = stats.fisher_exact([[significant, total - significant], [p_all - significant, t_all - p_all - total + significant]])

        if p_value < 0.05:
            plt.text(position, 91 * significant / total, '*', ha='center', va='bottom', color='tab:red', fontsize=22)

        plt.bar(position, 100 * significant / total, width=0.68, lw=0, color='tab:purple', alpha=0.6)

    plt.title(cancertype)
    plt.xticks(positions, xticklabels, rotation=45, ha='right', fontsize=13, color='k')
    plt.ylabel('% TFBS overlapped ATAC-seq \npeaks w/ logrank P<0.05', fontsize=13)
    if cancertype == 'COAD':
        plt.ylim([0, 8])
    plt.savefig(figname, bbox_inches='tight', pad_inches=0.1, dpi=600, transparent=True)
    plt.close()


# == Directories
indir = 'f5_atac_overlap_coBinding_TFBS'
outdir = 'f6_atac_overlap_coBinding_TFBS_figs'
os.makedirs(outdir, exist_ok=True)

# == Data Files
atac_file = '../../../f9_TF_condensates_V3/data/TCGA/tcga_atac.bed'
name_match_file = '../../../f9_TF_condensates_V3/data/TCGA/TCGA-ATAC_SE_cancerType_match.xlsx'
clinical_dir = '../../../f7_TF_condensates_test/f6_revised_TCGA_ATAC_cor_SE/f2_clinical_hicor_atac_peaks/f1_clinical_at_each_peak/f2_caseID_each_peak_vs_clinical/'

# == Load Cancer Type to Cell Type Mapping
name_match = pd.read_excel(name_match_file, index_col=0).dropna()

# == Load Selected Factors for Each Cell Type
selected_factors = {}
tfbs_cp_dir = '../f2_cor_CP_SE_AICAP/f9_per_CT_TFBS_CP_cor_zscore_CP_with_motif_SE/TFBS_CP/'

for ct in ['MCF-7', 'HCT-116']:
    df = pd.read_csv(f"{tfbs_cp_dir}/_CP_TFBS_nonBlackList_vs_TFMS_{ct}.csv", index_col=0)
    selected_factors[f"{ct} top_TFBSCP"] = df['TFBS CP rank'].sort_values().iloc[:3].index
    selected_factors[f"{ct} top_zscored_TFBSCP"] = df['avg rank'].sort_values().iloc[:3].index

# == Process Each Cancer Type
for cancertype in ['BRCA', 'COAD']:
    ct = name_match.loc[cancertype].SE
    clinical_df = pd.read_csv(f"{clinical_dir}/{cancertype}_logrank_info.csv", index_col=0)

    for treat_flag in ['percentile_T_ExtendMerge']:
        for factorType in ['top_TFBSCP', 'top_zscored_TFBSCP']:
            subdir = f"{ct}_{factorType}"
            os.makedirs(f"{outdir}/{subdir}", exist_ok=True)

            factors = selected_factors[f"{ct} {factorType}"]
            xticklabels = [
                'Union',
                '-'.join([factors[0]]),
                '-'.join([factors[1]]),
                '-'.join([factors[2]]),
                '-'.join(factors[:2]),
                '-'.join(factors[:3])
            ]

            # Add log-rank information
            out_df = pd.DataFrame()
            add_logrank_info(out_df, clinical_df, 'All', flag)

            for ii, label in enumerate(xticklabels):
                overlapped_bed = f"{indir}/{subdir}/atac_overlapped_{treat_flag}_{ct}_{label}.bed"
                overlapped_df = pd.read_csv(overlapped_bed, sep='\t', index_col=3, header=None)
                overlapped_df.columns = ['chr', 'start', 'end', 'score', 'annotation']

                # Save merged clinical information
                merged_df = pd.concat([overlapped_df, clinical_df], axis=1).dropna()
                merged_df.to_csv(f"{outdir}/{subdir}/{cancertype}_{treat_flag}_{label}_logRank.csv")

                clinical_df_tf = clinical_df.loc[overlapped_df.index]
                add_logrank_info(out_df, clinical_df_tf, label, flag)

            # Save and plot the log-rank information
            out_df.to_csv(f"{outdir}/{subdir}/{cancertype}_{treat_flag}_clinical.csv")
            figname = f"{outdir}/{subdir}/{cancertype}_{treat_flag}_clinical.pdf"
            stack_bar(out_df, ['All'] + xticklabels, figname, cancertype, flag)

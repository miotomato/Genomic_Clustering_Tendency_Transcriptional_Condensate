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
matplotlib.rcParams['font.size'] = 14
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"

sns.set(font_scale=1)
sns.set_style("whitegrid", {'axes.grid': False})
sns.set_style("ticks")


def mark_pvalue(compr_pos, positions, box_vals):
    """Perform t-test between two groups and mark significance with an asterisk if p < 0.05."""
    s, p = stats.ttest_ind(box_vals[compr_pos[1]], box_vals[compr_pos[0]], nan_policy='omit')
    x1, x2 = positions[compr_pos[0]], positions[compr_pos[1]]
    y = np.percentile(box_vals[compr_pos[0]], 97.5)
    if p < 0.05:
        plt.text(x2, y, '*', ha='center', va='center_baseline', color='r', fontsize=27)
    return s, p


# == Input and Output Directories
indir = 'f3_DCI_overlap_coBinding_TFBS'
outdir = 'f4_DCI_overlap_coBinding_TFBS_figs'
os.makedirs(outdir, exist_ok=True)

# == Load Cancer Type and Cell Type Mapping
name_match_file = '../../../f9_TF_condensates_V3/data/TCGA/TCGA-ATAC_SE_cancerType_match.xlsx'
name_match = pd.read_excel(name_match_file, index_col=0).dropna()

# == Load Selected Factors for Each Cell Type
selected_factors = {}
tfbs_cp_dir = '../../f1_TF_cluster_potential/f2_cor_CP_SE_AICAP/f9_per_CT_TFBS_CP_cor_zscore_CP_with_motif_SE/TFBS_CP/'

for ct in ['MCF-7', 'HCT-116']:
    df = pd.read_csv(f"{tfbs_cp_dir}/_CP_TFBS_nonBlackList_vs_TFMS_{ct}.csv", index_col=0)
    selected_factors[f"{ct} top_TFBSCP"] = df['TFBS CP rank'].sort_values().iloc[:3].index
    selected_factors[f"{ct} top_zscored_TFBSCP"] = df['avg rank'].sort_values().iloc[:3].index

# == DCI Data Directory and Parameters
dci_dir = '../../../f11_TF_condensates_KS_test/f3_public_data/f1_hct116_hic_RAD21_auxin/f1_bart3d/'
genomicDistances = [100000, 200000, 500000]
reps = ['rep1', 'rep2', 'all_reps']

# == Loop Over Cancer Types, Treatment Flags, and Factor Types
for cancertype in ['COAD']:
    ct = name_match.loc[cancertype].SE
    
    for treat_flag in ['percentile_T_ExtendMerge']:
        for factorType in ['top_zscored_TFBSCP']:
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

            # == Box Plot of DCI
            for genomicDis in genomicDistances[:1]: 
                for rep in reps[2:]:
                    box_vals = []
                    out_df = pd.DataFrame()

                    # == Load Genomic DCI Data
                    dci_file = f"{dci_dir}/RAD21_6hr_auxin_over_NT_{rep}_dis{int(genomicDis / 1000)}k_differential_score.bed"
                    dci_df = pd.read_csv(dci_file, sep='\t', header=None)
                    box_vals.append(dci_df[3])

                    # == Load DCI of Co-binding for Each Factor Combination
                    for ii, label in enumerate(xticklabels):
                        overlapped_bed = f"{indir}/{subdir}/DCI_RAD21_KO_{rep}_dis{int(genomicDis / 1000)}k_overlapped_{treat_flag}_{ct}_{label}.bed"
                        overlapped_df = pd.read_csv(overlapped_bed, sep='\t', header=None)
                        box_vals.append(overlapped_df[3])

                    # == Create Box Plot
                    figname = f"{cancertype}_{rep}_{int(genomicDis / 1000)}k_{treat_flag}_DCI"
                    positions = np.arange(len(xticklabels) + 1)

                    plt.figure(figsize=(3, 2.5))
                    plt.boxplot(
                        box_vals, positions=positions, widths=0.6, patch_artist=True,
                        boxprops=dict(color='k', facecolor='w', fill=None, lw=1),
                        medianprops=dict(color='grey'), showfliers=False
                    )

                    # Mark p-values for comparisons
                    for ii in range(1, len(xticklabels) + 1):
                        s, p = mark_pvalue([0, ii], positions, box_vals)
                        out_df.loc[xticklabels[ii - 1], 'ttest-s'] = s
                        out_df.loc[xticklabels[ii - 1], 'ttest-p'] = p

                    # Plot Aesthetics
                    plt.xticks(positions, ['All'] + xticklabels, fontsize=12, rotation=45, ha='right')
                    plt.axhline(y=0, color='k', lw=1.2, ls='--')
                    plt.ylabel('DCI of RAD21 \n Depletion in HCT-116', fontsize=14)
                    plt.title(cancertype)

                    # Save Figure
                    plt.savefig(f"{outdir}/{subdir}/{figname}.pdf", bbox_inches='tight', pad_inches=0.1, dpi=600, transparent=True)
                    plt.close()

                    # Save t-test Results
                    out_df.to_csv(f"{outdir}/{subdir}/{figname}.csv")

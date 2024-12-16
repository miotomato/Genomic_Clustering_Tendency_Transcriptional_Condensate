import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns
from collections import Counter

# == Configure Plotting Parameters
matplotlib.rcParams['font.size'] = 14
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams["mathtext.rm"] = "Arial"

sns.set(font_scale=1.1)
sns.set_style("whitegrid", {'axes.grid': False})
sns.set_style("ticks")

# == Function Definitions

def add_logrank_info(out_df, df, motif_name):
    """Add log-rank test results to the output DataFrame."""
    df_p = df[df['log rank p'] < p_thre]
    df_pt = df[(df['treat time'] < df['ctrl time']) & (df['log rank p'] < p_thre)]

    out_df.loc[motif_name, 'total'] = df.shape[0]
    out_df.loc[motif_name, f'#P<{p_thre}'] = df_p.shape[0]
    out_df.loc[motif_name, f'%P<{p_thre}'] = np.round(df_p.shape[0] / df.shape[0], 2)
    out_df.loc[motif_name, f'#TreatTime<CtrlTime&P<{p_thre}'] = df_pt.shape[0]
    out_df.loc[motif_name, f'%TreatTime<CtrlTime&P<{p_thre}'] = np.round(df_pt.shape[0] / df.shape[0], 2)

    return out_df

# == Input and Output Directories
indir = 'f4_clinical_per_gene_by_TFBS_RP'
outdir = 'f5_TF_targets_clinical_by_gene_TFBS_RP'
os.makedirs(outdir, exist_ok=True)

# == Load Matched Names Between TCGA ATAC and SE Data
name_match = pd.read_excel("data/TCGA-ATAC_SE_cancerType_match.xlsx", index_col=0).dropna()

# == Load Top 3 Factors for Each Cell Type
selected_factors = {}
tfbs_cp_dir = "TFBS_CP"

for ct in ['MCF-7', 'HCT-116']:
    df = pd.read_csv(f"{tfbs_cp_dir}/_CP_TFBS_nonBlackList_vs_TFMS_{ct}.csv", index_col=0)
    selected_factors[f"{ct} top_TFBSCP"] = df['TFBS CP rank'].sort_values().iloc[:3].index
    selected_factors[f"{ct} top_zscored_TFBSCP"] = df['avg rank'].sort_values().iloc[:3].index

# == Process Each Cancer Type and Treatment Flag
p_thres = [0.1, 0.05]
for p_thre in p_thres:
    for cancertype in ['BRCA', 'COAD']:
        ct = name_match.loc[cancertype].SE

        for treat_flag in ['percentile_T_ExtendMerge']:
            for factorType in ['top_zscored_TFBSCP']:
                subdir = f"{ct}_{factorType}"
                os.makedirs(f"{outdir}/{subdir}/fig", exist_ok=True)
                factors = selected_factors[f"{ct} {factorType}"]

                # Define x-tick labels
                xticklabels = [
                    'ALL', 'Union',
                    '-'.join([factors[0]]),
                    '-'.join([factors[1]]),
                    '-'.join([factors[2]]),
                    '-'.join(factors[:2]),
                    '-'.join(factors[:3])
                ]

                # Process each factor combination
                out_df = pd.DataFrame()
                for ii, label in enumerate(xticklabels):
                    basename = f"{treat_flag}_{ct}_{label}"
                    clinical_file = f"{indir}/{subdir}/{basename}_logrank_info.csv"
                    clinical_df = pd.read_csv(clinical_file, index_col=0)

                    out_df = add_logrank_info(out_df, clinical_df, label)

                # Drop 'Union' row
                out_df = out_df.drop(['Union'])

                # Plot results
                check_col = f"#P<{p_thre}"
                plt.figure(figsize=(3, 2.6))
                positions = np.arange(out_df.shape[0])

                t_all = out_df.loc['ALL', 'total']
                p_all = out_df.loc['ALL', check_col]

                for position in positions:
                    total = out_df.iloc[position]['total']
                    a = out_df.iloc[position][check_col]
                    s, p = stats.fisher_exact([[a, total - a], [p_all - a, t_all - p_all - total + a]])

                    out_df.loc[out_df.index[position], f"{check_col} fisher s"] = s
                    out_df.loc[out_df.index[position], f"{check_col} fisher p"] = p

                    if p < 0.05:
                        plt.text(position, 92 * a / total, '*', ha='center', va='bottom', color='tab:red', fontsize=20)

                    plt.bar(position, 100 * a / total, width=0.68, lw=0, color='tab:purple', alpha=0.6)

                # Plot aesthetics
                plt.title(cancertype)
                plt.xticks(positions, out_df.index, rotation=45, ha='right', fontsize=13, color='k')
                plt.ylabel(f'% ATAC-seq target genes \n w/ logrank P<{p_thre}', fontsize=13)

                # Save the plot
                figname = f"{outdir}/{subdir}/{treat_flag}_{ct}_p{p_thre}.pdf"
                plt.savefig(figname, bbox_inches='tight', pad_inches=0.1, dpi=600, transparent=True)
                plt.close()

                # Save the results DataFrame
                out_df.to_csv(f"{outdir}/{subdir}/{treat_flag}_{ct}_p{p_thre}_data.csv")

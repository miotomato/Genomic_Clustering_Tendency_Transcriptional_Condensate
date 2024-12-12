import os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
import warnings

# == Configure Plotting Parameters
warnings.filterwarnings("ignore")
matplotlib.rcParams['font.size'] = 16
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams["mathtext.rm"] = "Arial"

sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid': False})
sns.set_style("ticks", {'ytick.color': 'k', 'axes.edgecolor': 'k'})

# == Chromosomes List
chroms = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']

# == Directories
indir = '../../f12_KS_test_Rename/f4_mutation/f6_pancancer_mutation_data_collection'
outdir = 'f6b_pancancer_mutation_data_collection_fig'
os.makedirs(outdir, exist_ok=True)

# == Genome Control Data
genome_control = {
    'Rate': 13474.4 * 6285 / 3298912062,
    'Count': 76088647 / 3298912062
}

# == Cell Types
cts = ['MCF-7', 'HCT-116', 'HeLa', 'LNCaP', 'U87', 'HepG2']

# == Treatment Flags
treat_flags = ['percentile_T', 'percentile_T_ExtendMerge']

# == Rank Directory
rank_dir = '../f1_TF_cluster_potential/f2_cor_CP_SE_AICAP/f9_per_CT_TFBS_CP_cor_zscore_CP_with_motif_SE/TFBS_CP/'

# == Loop Over Mutation Types, Treatment Flags, and Cell Types
for mutationType in ['Rate']:  # Only processing 'Rate'
    for treat_flag in treat_flags[1:]:  # Only processing 'percentile_T_ExtendMerge'
        for ct in cts:
            # Load Rank Data for Current Cell Type
            rank_df = pd.read_csv(f"{rank_dir}/_CP_TFBS_nonBlackList_vs_TFMS_{ct}.csv", index_col=0)
            
            # Load Mutation Data for Current Cell Type
            mutation_file = f"{indir}/{ct}_murationRate.csv"
            df = pd.read_csv(mutation_file, index_col=0)
            df = df.drop(df.index.intersection(['T']))  # Drop 'T' row if present

            # Filter and Rank TFs Present in Both Mutation and Rank Data
            cp_ranked_tf = [i for i in rank_df.index if i in df.index]
            df = df.loc[cp_ranked_tf]

            # Data Columns for Mutation Rates
            data_col = [
                f"{treat_flag} %mutation{mutationType} per bp",
                f"percentile_C %mutation{mutationType} per bp"
            ]

            # Plotting Bar Graph
            positions = np.arange(df.shape[0])
            plt.figure(figsize=(df.shape[0] / 3.5, 2.6))

            # Plot Bars for C-TFBS and NC-TFBS
            a = plt.bar(positions - 0.2, df[data_col[0]], width=0.3, color='salmon', label='C-TFBS')
            b = plt.bar(positions + 0.2, df[data_col[1]], width=0.3, color='royalblue', label='NC-TFBS')

            # Add Genome Control Line
            plt.axhline(y=genome_control[mutationType], color='k', lw=1.2, ls='--')

            # Configure Plot Labels and Aesthetics
            plt.ylabel('Mutation Rate')
            plt.xlim([-1, df.shape[0]])
            plt.xticks(positions, df.index, rotation=60, ha='right', fontsize=12, color='k')
            plt.legend(fontsize=10, bbox_to_anchor=[1, 1.2], loc="upper right", frameon=False)

            # Save the Figure
            figname = f"{outdir}/{ct}_{treat_flag}_mutation{mutationType}.pdf"
            plt.savefig(figname, bbox_inches='tight', pad_inches=0.1, dpi=600, transparent=True)
            plt.close()

            # Perform T-Test and Print Results
            s, p = stats.ttest_ind(df[data_col[0]], df[data_col[1]])
            print(f"T-test for {ct}, {treat_flag}, {mutationType}: Statistic={s:.4f}, P-value={p:.4e}")

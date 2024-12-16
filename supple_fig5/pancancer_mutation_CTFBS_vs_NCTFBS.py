import os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

# Plot settings
matplotlib.rcParams['font.size'] = 16
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid': False})
sns.set_style("ticks", {'ytick.color': 'k', 'axes.edgecolor': 'k'})
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams["mathtext.rm"] = "Arial"

# Chromosome list
chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
          'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
          'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']

# Input and output directories
indir = 'f6_pancancer_mutation_data_collection'
outdir = 'f6b_pancancer_mutation_data_collection_fig'
os.makedirs(outdir, exist_ok=True)

# Genome control values for normalization/comparison
genome_control = {'Rate': 13474.4 * 6285 / 3298912062,
                  'Count': 76088647 / 3298912062}

# Cell types to analyze
cts = ['MCF-7', 'HCT-116', 'HeLa', 'LNCaP', 'U87', 'HepG2']

# Treatment flags
treat_flags = ['percentile_T', 'percentile_T_ExtendMerge']

# Rank directory for TFBS data
rank_dir = 'TFBS_CP'

# Loop over mutation types, treatment flags, and cell types
for mutationType in ['Rate', 'Count'][:1]:
    for treat_flag in treat_flags[1:]:
        for ct in cts:
            # Load TFBS ranking data
            rank_df = pd.read_csv(f'{rank_dir}/_CP_TFBS_nonBlackList_vs_TFMS_{ct}.csv', index_col=0)

            # Load mutation rate data
            df = pd.read_csv(f'{indir}/{ct}_murationRate.csv', index_col=0)
            df = df.drop(df.index.intersection(['T']))

            # Filter and sort based on TFBS ranking
            cp_ranked_tf = [i for i in rank_df.index if i in df.index]
            df = df.loc[cp_ranked_tf]

            # Define the columns for mutation rates
            data_col = [f'{treat_flag} %mutation{mutationType} per bp',
                        f'percentile_C %mutation{mutationType} per bp']

            # Plot the data
            positions = np.arange(df.shape[0])
            plt.figure(figsize=(df.shape[0] / 3.5, 2.6))
            a = plt.bar(positions - 0.2, df[data_col[0]], width=0.3, color='salmon')
            b = plt.bar(positions + 0.2, df[data_col[1]], width=0.3, color='royalblue')
            plt.xlim([-1, df.shape[0]])
            plt.axhline(y=genome_control[mutationType], color='k', lw=1.2, ls='--')
            plt.ylabel('Mutation Rate')
            plt.legend([a, b], ['C-TFBS', 'NC-TFBS'], fontsize=10, bbox_to_anchor=[1, 1.2],
                       borderaxespad=0.1, labelspacing=0.1, loc="upper right", frameon=False)
            plt.xticks(positions, df.index, rotation=60, ha='right', fontsize=12, color='k')

            # Save the figure
            figname = f'{outdir}/{ct}_{treat_flag}_mutation{mutationType}.pdf'
            plt.savefig(figname, bbox_inches='tight', pad_inches=0.1, dpi=600, transparent=True)
            plt.close()

            # Perform t-test and print results
            s, p = stats.ttest_ind(df[data_col[0]], df[data_col[1]])
            print(f'{ct} - {treat_flag} - {mutationType}: t-statistic = {s:.3f}, p-value = {p:.3e}')

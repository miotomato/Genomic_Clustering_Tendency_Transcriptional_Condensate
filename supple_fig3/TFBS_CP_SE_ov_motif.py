import os
import glob
import re
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sns
import scipy

# Plot settings
matplotlib.rcParams['font.size'] = 14
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
sns.set(font_scale=1.1)
sns.set_style("whitegrid", {'axes.grid': False})
sns.set_style("ticks")

# Suppress warnings
import warnings
warnings.filterwarnings("ignore")

# Configuration dictionaries
group_info = {
    'CP_KStest': {'data_median': ['log10-dis ks_2samp-s signed', 'log10-dis ks_2samp-p signed'],
                  'data_fisher': ['log10-dis ks_2samp-s signed', 'log10-dis ks_2samp-p signed Fisher-TF-Celltype-combined']},
    'CP_Ttest': {'data_median': ['log10-dis T-test-s', 'log10-dis T-test-p'],
                 'data_fisher': ['log10-dis T-test-s', 'log10-dis T-test-p Fisher-TF-Celltype-combined']},
    'CP_RankSum': {'data_median': ['log10-dis Wilcoxon-rank-sum-s', 'log10-dis Wilcoxon-rank-sum-p'],
                   'data_fisher': ['log10-dis Wilcoxon-rank-sum-s', 'log10-dis Wilcoxon-rank-sum-p Fisher-TF-Celltype-combined']},
    'SE': {'data_median': ['fisher_exact_s', 'fisher_exact_p'],
           'data_fisher': ['fisher_exact_s', 'fisher_exact_p Fisher-TF-Celltype-combined']}
}

or_thre = {'CP_KStest': 0.3, 'CP_Ttest': 60, 'CP_RankSum': 60, 'SE': 4}
label_name = {'CP_KStest': 'TFBS CP', 'CP_Ttest': 'TFBS CP', 'CP_RankSum': 'TFBS CP', 'SE': 'log2 OR'}
cmp_dict = {'CP_KStest': plt.cm.PiYG_r, 'CP_Ttest': plt.cm.PiYG_r, 'CP_RankSum': plt.cm.PiYG_r, 'SE': plt.cm.bwr}

# Cutoffs
index_cutoff = 20
column_cutoff = 7

# Paths
SE_files = glob.glob('data/SE_hg38/*.bed')
SEs = [os.path.basename(i).split(".bed")[0] for i in SE_files]
indirs = ['CP_TFBS_nonBlackList_vs_TFMS',
          'CP_TFBS_nonBlackList_overlap_motif_vs_TFMS',
          'CP_TFBS_nonBlackList_NOT_overlap_motif_vs_TFMS']
fisher_P_key = 'data_fisher'
profileTypes = ['with_motif_SE', 'with_motif_only']

# Main loop
for profileType in profileTypes:
    outdir = f'f3_TFBS_CP_heatmap_{profileType}'
    os.makedirs(f'{outdir}/_csv', exist_ok=True)

    for indir in indirs:
        infile = f'{indir}/per_TF_per_Celltype_{indir}.csv'
        df = pd.read_csv(infile, index_col=0).replace(np.inf, np.nan).dropna()
        correlation_data = {}
        correlation_p = {}

        for or_key in ['CP_KStest', 'SE']:
            s_col = group_info[or_key][fisher_P_key][0]
            p_col = group_info[or_key][fisher_P_key][1]
            df_s = pd.DataFrame()
            df_p = pd.DataFrame()
            basename = f'{fisher_P_key}_{indir}_{or_key}'

            for ii in df.index:
                factor, celltype = ii.split('_')[1], ii.split('_')[0]
                if profileType == 'with_motif_SE' and celltype not in SEs:
                    continue
                df_s.loc[factor, celltype] = df.loc[ii, s_col]
                df_p.loc[factor, celltype] = df.loc[ii, p_col]

            df_s.rename(index={'RBPJ': 'NOTCH1'}, inplace=True)
            df_p.rename(index={'RBPJ': 'NOTCH1'}, inplace=True)
            df_s.to_csv(f'{outdir}/_csv/{basename}_statistics.csv')
            df_p.to_csv(f'{outdir}/_csv/{basename}_P.csv')

            # Z-score normalization
            df_s_z = df_s.apply(lambda x: scipy.stats.zscore(x, nan_policy='omit') if len(x.dropna()) > 1 else x, axis=1)
            df_s_z.to_csv(f'{outdir}/_csv/{basename}_statistics_zscored.csv')

            # Filter data
            selected_index = df_s.notnull().sum(axis=1).nlargest(index_cutoff).index[::-1]
            selected_col = df_s.notnull().sum(axis=0).nlargest(column_cutoff).index
            df_s = df_s.loc[selected_index, selected_col]
            df_p = df_p.loc[selected_index, selected_col].replace(0, 1e-299)

            if or_key == 'SE':
                df_s = np.log2(df_s)

            # Save filtered data
            df_s.to_csv(f'{outdir}/{basename}_S.csv')
            df_p.to_csv(f'{outdir}/{basename}_P.csv')
            correlation_data[or_key] = df_s
            correlation_p[or_key] = df_p

            # Plot heatmap
            plt.figure(figsize=(3.6, 7.5))
            gs = gridspec.GridSpec(3, 4, width_ratios=[4, 0.0, 0.0, 0.4], height_ratios=[1, 1, 0.8], wspace=0.15, hspace=0.5)
            pal = cmp_dict[or_key]
            vmax = or_thre[or_key]
            vmin = -vmax
            norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
            color_map = matplotlib.cm.ScalarMappable(norm=norm, cmap=pal)
            scatter_scale = 0.9

            ax = plt.subplot(gs[:, 0])
            for ii, idx in enumerate(df_s.index):
                for jj, col in enumerate(df_s.columns):
                    statistics = df_s.loc[idx, col]
                    pvalue = df_p.loc[idx, col]
                    if not np.isnan(statistics):
                        color = color_map.to_rgba(statistics)
                        size = -np.log10(pvalue)
                        ax.scatter(jj, ii, s=scatter_scale * size, color=color)
                    else:
                        ax.scatter(jj, ii, marker='x', s=scatter_scale * 30, color='k')

            ax.set_yticks(range(len(df_s.index)))
            ax.set_yticklabels(df_s.index, fontsize=13)
            ax.set_xticks(range(len(df_s.columns)))
            ax.set_xticklabels(df_s.columns, rotation=60, fontsize=13)
            ax.set_xlim([-0.6, len(df_s.columns) - 0.4])
            ax.set_ylim([-0.6, len(df_s.index) - 0.4])
            ax.set_title(f'{or_key}\n{fisher_P_key}\n{indir}', fontsize=13, va='bottom')

            # Save heatmap
            plt.savefig(f'{outdir}/{basename}.pdf', bbox_inches='tight', pad_inches=0.02, transparent=True)
            plt.close()

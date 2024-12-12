import os
import glob
import re
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import gridspec

# Suppress warnings for cleaner output
import warnings
warnings.filterwarnings("ignore")

# Set plotting styles
plt.rcParams['font.size'] = 14
sns.set(font_scale=1.0)
sns.set_style("whitegrid", {'axes.grid': False})
sns.set_style("ticks")
plt.rcParams["font.sans-serif"] = ["Arial"]

# Group information for different statistical tests
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

# Cutoff thresholds for statistical tests
or_thre = {'CP_KStest': 0.3, 'CP_Ttest': 60, 'CP_RankSum': 60, 'SE': 4}

# Labels for plots
label_name = {'CP_KStest': 'TFBS CP', 'CP_Ttest': 'TFBS CP', 'CP_RankSum': 'TFBS CP', 'SE': 'log2 OR'}

# Colormap dictionary
cmp_dict = {'CP_KStest': plt.cm.PiYG_r, 'CP_Ttest': plt.cm.PiYG_r, 'CP_RankSum': plt.cm.PiYG_r, 'SE': plt.cm.bwr}

# Parameters for filtering data
index_cutoff = 20
column_cutoff = 7

# SE files for filtering based on cell types
SE_files = glob.glob('../../../f12_KS_test_Rename/data/SE_hg38/*.bed')
SEs = [os.path.basename(i).split(".bed")[0] for i in SE_files]

# Input directories
indirs = [
    'CP_TFBS_nonBlackList_vs_TFMS',
    'CP_TFBS_nonBlackList_overlap_motif_vs_TFMS',
    'CP_TFBS_nonBlackList_NOT_overlap_motif_vs_TFMS'
]

# Profile types for analysis
profileTypes = ['with_motif_SE', 'with_motif_only']

# Process each profile type
for profileType in profileTypes:
    outdir = f'f3_TFBS_CP_heatmap_{profileType}'
    os.makedirs(f'{outdir}/_csv', exist_ok=True)

    for indir in indirs:
        infile = f'../f1_TFMS_TFBS_CP/{indir}/per_TF_per_Celltype_{indir}.csv'
        df = pd.read_csv(infile, index_col=0)
        df = df.replace(np.inf, np.nan).dropna()
        
        # Process CP and SE enrichment data
        for or_key in ['CP_KStest', 'SE']:
            s_col = group_info[or_key]['data_fisher'][0]
            p_col = group_info[or_key]['data_fisher'][1]
            df_s = pd.DataFrame()
            df_p = pd.DataFrame()

            basename = f'{indir}_{or_key}'
            print(f'Processing: {basename}')

            for ii in df.index:
                factor = ii.split('_')[1]
                celltype = ii.split('_')[0]
                if profileType == 'with_motif_SE' and celltype not in SEs:
                    continue
                df_s.loc[factor, celltype] = df.loc[ii, s_col]
                df_p.loc[factor, celltype] = df.loc[ii, p_col]

            # Save processed statistics
            df_s.to_csv(f'{outdir}/_csv/{basename}_statistics.csv')
            df_p.to_csv(f'{outdir}/_csv/{basename}_P.csv')

            # Plot heatmap
            plt.figure(figsize=(3.6, 7.5))
            gs = gridspec.GridSpec(3, 4, width_ratios=[4, 0.0, 0.0, 0.4], height_ratios=[1, 1, 0.8], wspace=0.15, hspace=0.5)
            pal = cmp_dict[or_key]
            vmax = or_thre[or_key]
            vmin = -vmax
            norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
            color_map = matplotlib.cm.ScalarMappable(norm=norm, cmap=pal)

            ax = plt.subplot(gs[:, 0])
            for ii in range(len(df_s.index)):
                for jj in range(len(df_s.columns)):
                    statistics = df_s.loc[df_s.index[ii], df_s.columns[jj]]
                    pvalue = df_p.loc[df_s.index[ii], df_s.columns[jj]]
                    if not np.isnan(statistics):
                        color = color_map.to_rgba(statistics)
                        size = -1 * np.log10(pvalue)
                        ax.scatter(jj, ii, s=size, color=color)
                    else:
                        ax.scatter(jj, ii, marker='x', s=30, color='k')

            ax.set_yticks(np.arange(len(df_s.index)))
            ax.set_yticklabels(df_s.index, fontsize=13)
            ax.set_xticks(np.arange(len(df_s.columns)))
            ax.set_xticklabels(df_s.columns, rotation=60, fontsize=13)
            ax.set_xlim([-.6, len(df_s.columns) - .4])
            ax.set_ylim([-.6, len(df_s.index) - .4])
            ax.set_title(f'{or_key}\\n{indir}', fontsize=13)

            plt.savefig(f'{outdir}/{basename}.pdf', bbox_inches='tight', pad_inches=0.02, transparent=True)
            plt.close()

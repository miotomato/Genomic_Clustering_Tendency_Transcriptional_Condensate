import os
import re
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

# == Matplotlib and Seaborn Configuration ==
plt.rcParams['font.size'] = 12
plt.rcParams["font.sans-serif"] = ["Arial"]
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams["mathtext.rm"] = "Arial"
sns.set(font_scale=1.1)
sns.set_style("whitegrid", {'axes.grid': False})
sns.set_style("ticks")


# == Scatter Plot Function ==
def scatter_plot(x, y, xlabel, ylabel, figname, title):
    """
    Generates a scatter plot with linear regression and density visualization.

    Args:
        x (Series): X-axis values.
        y (Series): Y-axis values.
        xlabel (str): Label for the X-axis.
        ylabel (str): Label for the Y-axis.
        figname (str): Output file path for the figure.
        title (str): Plot title.

    Returns:
        None
    """
    plt.figure(figsize=(1.8, 2.5))
    plt.scatter(x, y, c='k', s=9)

    # Linear regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    x_sort = np.sort(x)
    plt.plot(x_sort, x_sort * slope + intercept, c='grey', ls='--', lw=0.6)
    plt.text(0.04, 0.88, f'$R^2$ = {r_value**2:.2f}', fontsize=11, transform=plt.gca().transAxes)

    # Plot settings
    plt.axhline(y=0, color='grey', lw=0.5, ls='--')
    plt.axvline(x=0, color='grey', lw=0.5, ls='--')
    plt.ylim([0, 135])
    plt.xlim([-230, 230])
    plt.xticks([-200, 0, 200], ['-200'.rjust(9), '0', '200'.ljust(6)])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)

    if title != 'GM12878':
        plt.yticks([])
        plt.ylabel('')

    plt.savefig(figname, bbox_inches='tight', pad_inches=0.02, transparent=True)
    plt.show()
    plt.close()


# == Output Directory ==
outdir = 'f11_TFBS_CP_cor_peak_nums'
os.makedirs(outdir, exist_ok=True)

# == Read TFBS CP Data ==
tfbp_cp_df = pd.read_csv('../f1_TFMS_TFBS_CP/CP_TFBS_nonBlackList_vs_TFMS/per_data_CP_TFBS_nonBlackList_vs_TFMS.csv', index_col=0)
tfbs_dir = 'f3_TFBS_CP_heatmap_with_motif_SE/_csv/'

# Different types of CP comparisons
tfbs_cp_types = [
    'CP_TFBS_nonBlackList_vs_TFMS',
    'CP_TFBS_nonBlackList_overlap_motif_vs_TFMS',
    'CP_TFBS_nonBlackList_NOT_overlap_motif_vs_TFMS'
]

# Statistical tests to use
test_types = ['T-test-s', 'Wilcoxon-rank-sum-s', 'ks_2samp-s signed']
tfbs_cp_type = 'CP_TFBS_nonBlackList_vs_TFMS'

# == Process and Plot Data ==
for test_type in test_types[2:]:
    basename = f'data_fisher_{tfbs_cp_type}_CP_KStest'
    tfbs_s = pd.read_csv(f'{tfbs_dir}/{basename}_statistics.csv', index_col=0)
    tfbs_p = pd.read_csv(f'{tfbs_dir}/{basename}_P.csv', index_col=0)

    # Iterate through each cell type
    for ct in tfbs_s.columns:
        y = tfbs_s[ct].dropna()
        if len(y) < 2:
            continue

        outdf = pd.DataFrame()
        plt.figure(figsize=(2.6, 2.6))
        color_i = 0

        # Iterate through each factor
        for factor in y.index:
            kept_data = [i for i in tfbp_cp_df.index if re.search(f'{ct}_{factor}_', i)]
            if len(kept_data) < 5 or factor == 'T':
                continue

            # Prepare data for plotting
            px = tfbp_cp_df.loc[kept_data, '#TFBS_nonBlackList']
            px = np.log10(px)
            py = tfbp_cp_df.loc[kept_data, f'log10-dis {test_type}']

            # Perform linear regression
            slope, intercept, r_value, p_value, std_err = stats.linregress(px, py)
            plt.scatter(px, py, color=plt.cm.tab20(color_i), s=9, label=f'{factor}')
            outdf.loc[factor, 'r-value'] = r_value
            outdf.loc[factor, 'p-value'] = p_value
            color_i += 1
            print(factor, py.min(), py.max())

        # Plot settings
        plt.legend(bbox_to_anchor=[1, 1], markerscale=2, fontsize=11, borderaxespad=0.2, labelspacing=0.2,
                   handletextpad=0.2, handlelength=1., loc="upper left", frameon=False)
        if ct == 'MCF-7':
            plt.ylim([-0.15, 0.5])
        else:
            plt.ylim([0, 0.35])

        plt.xlabel('log10 #TFBS')
        plt.ylabel('TFBS CP')
        plt.title(ct)

        if color_i >= 2:
            plt.savefig(f'{outdir}/per_ct_{ct}_{test_type}_no_rp.pdf', bbox_inches='tight', pad_inches=0.02, transparent=True)
            plt.close()

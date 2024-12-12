import os
import glob
import re
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for saving plots
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.interpolate import interpn
from scipy.stats import gaussian_kde
from matplotlib import gridspec

# Set plotting styles
plt.rcParams['font.size'] = 9
sns.set(font_scale=0.9)
sns.set_style("whitegrid", {'axes.grid': False})
sns.set_style("ticks")
plt.rcParams["font.sans-serif"] = ["Arial"]
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams["mathtext.rm"] = "Arial"

def box_plot(pos_len, vals, or_cutoff, mark_p=False):
    """
    Creates a box plot for the given values.

    Parameters:
    - pos_len: Number of positions to plot.
    - vals: Array of values to plot.
    - or_cutoff: Log2 odds ratio cutoff for comparison.
    - mark_p: Boolean to indicate if p-values should be marked on the plot.
    """
    positions = np.arange(pos_len)
    box_step = int(len(vals) / len(positions))
    box_vals = []

    for ii in positions:
        box_val = vals[ii * box_step:(ii + 1) * box_step]
        box_vals.append(box_val)

        # Perform a one-sample t-test against the log2 OR cutoff
        s, p = stats.ttest_1samp(box_val, np.log2(or_cutoff), alternative='greater')
        print(len(box_val), s, p)
        p_val = '*' if p < 0.05 else ''

        # Mark p-values if specified
        if mark_p:
            ax.text(ii, np.percentile(vals, 96), p_val, fontsize=27, color='red', ha='center', va='center')

    # Create the box plot
    plt.boxplot(box_vals, positions=positions, widths=0.6, patch_artist=True,
                boxprops=dict(color='k', facecolor='w', fill=None, lw=1),
                medianprops=dict(color='grey'), showfliers=False)

# ==== Create output directory
outdir = 'f2_TFMS_CP_cor_SE_heatmap'
os.makedirs(outdir, exist_ok=True)

# ==== Define test types
test_types = ['KS_statistics']

for test_type in test_types:
    # Read input CSV data
    df = pd.read_csv('data/TFMS_CP_SE_enrich.csv', index_col=0)

    # Sort the DataFrame based on the specified test type
    if test_type == 'KS_statistics':
        df = df.sort_values(by='log10-dis ks_2samp-s signed', ascending=False)

    # Define a dictionary of values to plot
    value_dic = {
        'Motif_num': df['#TFMS'],
        'Motif_length': df['len-of-TFMS'],
        'KS_statistics': df['log10-dis ks_2samp-s signed'],
        'std': df['log10-dis ks_2samp-s signed std'],
        'log2_OR_at_SE': np.log2(df['enrich-at-SE-fisher-exact-s'])
    }

    # Create the figure and define grid layout
    plt.figure(figsize=(6, 4))
    width_ratios = [9, 0.2]
    height_ratios = [1, 1, 1]
    gs = gridspec.GridSpec(3, 2, width_ratios=width_ratios, height_ratios=height_ratios, wspace=0.1, hspace=0.1)

    # ==== Plot TFMS Cluster Potential (CP) scatter plot
    ax = plt.subplot(gs[0, 0])
    vals = value_dic[test_type]
    ax.scatter(np.arange(len(vals)), vals, s=2, c='k')
    ax.axhline(y=0, color='k', lw=1.2, ls='--')
    ax.set_xlim([-.5, len(vals) - .5])
    ax.xaxis.tick_top()
    ax.set_ylabel('TFMS CP', ha='center')

    # ==== Color bar for the heatmap
    cmap = plt.cm.bwr
    vmax = 1
    vmin = -1 * vmax
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)

    ax = plt.subplot(gs[1, 1])
    cbar = matplotlib.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, orientation='vertical')
    cbar.set_ticks([vmin, 0, vmax])
    cbar.set_label('log2 OR', labelpad=-30, fontsize=9, va='center')
    ax.tick_params(axis='y', direction='out', length=2, width=1, colors='black')

    # ==== Heatmap of SE enrichment
    ax = plt.subplot(gs[1, 0])
    vals = value_dic['log2_OR_at_SE']
    vals = np.transpose(vals.to_frame())
    sns.heatmap(vals, cmap=cmap, cbar=False, vmax=vmax, vmin=vmin, xticklabels=False, yticklabels=False, ax=ax)
    ax.set_ylabel('Enrichment at \n SE (log2 OR)\n', ha='center')

    # ==== Box plot of SE enrichment
    ax = plt.subplot(gs[2, 0])
    vals = value_dic['log2_OR_at_SE']
    or_cutoff = 1.2
    box_plot(20, vals, or_cutoff, mark_p=True)
    ax.axhline(y=np.log2(or_cutoff), color='k', lw=1.2, ls='--')
    ax.set_xticks([])
    ax.set_ylabel('log2 OR', ha='center')

    # ==== Save the figure
    plt.savefig(f'{outdir}/CP_by_{test_type}.pdf', bbox_inches='tight', pad_inches=0.02, transparent=True)
    plt.close()
    print(f'Figure saved to {outdir}/CP_by_{test_type}.pdf')

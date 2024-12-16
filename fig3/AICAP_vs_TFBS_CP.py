import os
import sys
import glob
import re
import bisect
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.interpolate import interpn
from scipy.stats import gaussian_kde
from matplotlib import gridspec
import warnings

# Suppress warnings for cleaner output
warnings.filterwarnings("ignore")

# Set plotting styles
matplotlib.rcParams['font.size'] = 14
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid': False})
sns.set_style("ticks")
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams["mathtext.rm"] = "Arial"

def box_plot(pos_len, df, mark_p=False):
    """
    Create a box plot for the given values.

    Parameters:
    - pos_len: Number of positions to plot.
    - df: DataFrame containing values to plot.
    - mark_p: Boolean to indicate whether to mark p-values.
    """
    positions = np.arange(pos_len)
    box_step = int(len(vals) / len(positions))
    box_vals = []
    for ii in positions:
        box_val = vals[ii * box_step:(ii + 1) * box_step]
        box_vals.append(box_val)
        s, p = stats.ttest_1samp(box_val, 0)
        p_val = '*' if p < 0.05 else ''
        if mark_p:
            ax.text(ii, np.mean(box_val), p_val, fontsize=30, color='r', ha='center', va='center')
    plt.boxplot(box_vals, positions=positions, widths=0.6, patch_artist=True,
                boxprops=dict(color='k', facecolor='w', fill=None, lw=1),
                medianprops=dict(color='grey'), showfliers=False)

def scatter_plot(x, y, xlabel, ylabel, figname, title):
    """
    Create a scatter plot with annotations.

    Parameters:
    - x: X-axis values.
    - y: Y-axis values.
    - xlabel: Label for the x-axis.
    - ylabel: Label for the y-axis.
    - figname: File path to save the figure.
    - title: Title of the plot.
    """
    plt.figure(figsize=(2.2, 2.5))
    plt.scatter(x, y, c='k', s=9)
    
    # Annotate the points with the lowest and highest x values
    label_i = 0
    label_indexes = x[x < 0].sort_values(ascending=True)[:half_len].index
    for label_index in label_indexes:
        plt.scatter(x[label_index], y[label_index], c=colors1[label_i], s=25, label=label_index)
        plt.text(x[label_index], y[label_index], label_index)
        label_i += 1

    label_i = 0
    label_indexes = x[x > 0].sort_values(ascending=False)[:half_len].index
    for label_index in label_indexes:
        plt.scatter(x[label_index], y[label_index], c=colors2[label_i], s=25, label=label_index)
        plt.text(x[label_index], y[label_index], label_index)
        label_i += 1

    plt.axvline(x=0, color='grey', lw=0.5, ls='--')
    plt.xlabel(xlabel)
    plt.title(title, ha='center')
    plt.ylim([ylim_min, ylim_max])
    plt.xlim([0, 10])
    plt.savefig(figname, bbox_inches='tight', pad_inches=0.02, transparent=True)
    plt.close()

def return_colors(pal, half_len, a, b):
    """
    Generate a list of colors based on a colormap.

    Parameters:
    - pal: Colormap to use.
    - half_len: Number of colors to generate.
    - a: Minimum normalization value.
    - b: Maximum normalization value.

    Returns:
    - List of RGBA color values.
    """
    norm = matplotlib.colors.Normalize(vmin=a * half_len, vmax=b * half_len)
    color_map = matplotlib.cm.ScalarMappable(norm=norm, cmap=pal)
    return [color_map.to_rgba(i) for i in np.arange(half_len)]

# ==== Main script ====

# Read AICAP data
idrfile = 'data/public/13059_2021_2456_MOESM2_ESM.xlsx'
idrdf = pd.read_excel(idrfile, sheet_name='condition 2(1,6-HD-2)', index_col=0)

# TFBS info
tfbs_cp_types = ['CP_TFBS_nonBlackList_vs_TFMS',
                 'CP_TFBS_nonBlackList_overlap_motif_vs_TFMS',
                 'CP_TFBS_nonBlackList_NOT_overlap_motif_vs_TFMS']
tfbs_cp_type = tfbs_cp_types[0]

half_len = 5
pal1 = plt.cm.RdPu_r
pal2 = plt.cm.GnBu
colors1 = return_colors(pal1, half_len, -0.2, 1.2)
colors2 = return_colors(pal2, half_len, -0.5, 1.1)
ylim_min = -0.05
ylim_max = 0.35

profileTypes = ['with_motif_SE', 'with_motif_only']
for profileType in profileTypes:
    tfbs_dir = f'f3_TFBS_CP_heatmap_{profileType}/_csv/'
    outdir = f'f8_per_CT_cor_AICAP_cor_CP_{profileType}'
    os.makedirs(f'{outdir}/_csv', exist_ok=True)
    
    tfbs_cp_s = pd.read_csv(f'{tfbs_dir}/data_fisher_{tfbs_cp_type}_CP_KStest_statistics.csv', index_col=0)
    tfbs_se_s = pd.read_csv(f'{tfbs_dir}/data_fisher_{tfbs_cp_type}_SE_statistics.csv', index_col=0)
    
    for ct in tfbs_cp_s.columns:
        y = tfbs_cp_s[ct].dropna()
        x = idrdf.loc[y.index.intersection(idrdf.index), 'AICAP']
        if len(x) < 1:
            continue
        
        df = pd.concat([y, tfbs_se_s[ct].loc[y.index], np.log2(x)], axis=1)
        df.columns = ['TFBS CP', 'log2 OR', 'log2 AICAP']
        df.sort_values(by='TFBS CP', ascending=False).round(2).to_csv(f'{outdir}/_csv/{ct}.csv')
        
        x = -1 * np.log2(x[x < 1].dropna())
        y = y.loc[x.index]
        if len(x) < 3:
            continue
        
        xlabel = '-log2 AICAP'
        ylabel = 'TFBS CP'
        figname = f'{outdir}/{tfbs_cp_type}_{ct}.pdf'
        
        fig = plt.figure(figsize=(3, 3))
        plt.scatter(x, y, c='k', s=9)
        plt.axvline(x=0, color='grey', lw=0.5, ls='--')
        plt.xlabel(xlabel)
        plt.title(ct, ha='center', va='bottom')
        plt.ylim([ylim_min, ylim_max])
        plt.xlim([0, 10])
        plt.savefig(figname, bbox_inches='tight', pad_inches=0.02, transparent=True)
        plt.close()

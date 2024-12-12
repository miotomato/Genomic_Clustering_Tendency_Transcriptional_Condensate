import os
import glob
import re
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

# Set plotting styles
plt.rcParams['font.size'] = 12
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid': False})
sns.set_style("ticks")
plt.rcParams["font.sans-serif"] = ["Arial"]
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams["mathtext.rm"] = "Arial"

def mark_pvalue(compr_pos, positions, box_vals):
    """
    Marks the p-value on a plot for two compared positions.

    Parameters:
    - compr_pos: List specifying the positions to compare and the direction ('t' or other).
    - positions: Positions on the x-axis.
    - box_vals: List of values for each position.
    """
    s, p = stats.ttest_ind(box_vals[compr_pos[0]], box_vals[compr_pos[1]], nan_policy='omit')
    y, h, col = np.percentile(np.append(box_vals[compr_pos[0]], box_vals[compr_pos[1]]), 95) * 0.95, 1.05, 'k'
    y2 = np.percentile(np.append(box_vals[compr_pos[0]], box_vals[compr_pos[1]]), 5) * 0.99
    x1, x2 = positions[compr_pos[0]], positions[compr_pos[1]]
    indicator = 'down' if s > 0 else 'up'
    p_label = f'{indicator} \\n {p:.1e}'

    if p_label[-2] == '0':
        p_label = p_label[:-2] + p_label[-1]
    if p >= 0.05:
        p_label = 'n.s.'

    if compr_pos[2] == 't':
        plt.plot([x1 * 1.03, x1 * 1.03, x2 * 0.97, x2 * 0.97], [y, y * h, y * h, y], lw=1, c=col)
        plt.text((x1 + x2) * 0.5, y * h, p_label, ha='center', va='bottom', color=col, fontsize=12)
    else:
        plt.plot([x1 * 1.03, x1 * 1.03, x2 * 0.97, x2 * 0.97], [y2, y2 * 0.91, y2 * 0.91, y2], lw=1, c=col)
        plt.text((x1 + x2) * 0.5, y2 * 0.95, p_label, ha='center', va='top', color=col, fontsize=12)

def scatter_plot(df, xlabel, ylabel, figname, title):
    """
    Creates a scatter plot with annotations and regression.

    Parameters:
    - df: DataFrame containing data for the plot.
    - xlabel: Column name for x-axis values.
    - ylabel: Column name for y-axis values.
    - figname: File path to save the figure.
    - title: Title of the plot.
    """
    x = df[xlabel]
    y = df[ylabel]

    plt.figure(figsize=(2.2, 2.2))
    plt.scatter(x, y, c='k', s=9)

    # Annotate top TFBS based on TFBS CP rank
    top_y = df['TFBS CP'].sort_values(ascending=False)[:5].index
    for label_index in top_y:
        plt.scatter(x[label_index], y[label_index], c='tab:red', s=30, label=label_index)
        plt.annotate(label_index, xy=(x[label_index], y[label_index]), xytext=(x[label_index] + 0.5, y[label_index]),
                     arrowprops=dict(color='k', lw=0.7, arrowstyle="-"), color='k')

    # Regression line
    slope, intercept, r_value, _, _ = stats.linregress(x, y)
    plt.text(0.6, 0.06, f'$r = {r_value:.2f}$', fontsize=11, transform=plt.gca().transAxes)

    plt.axhline(y=0, color='grey', lw=0.7, ls='--')
    plt.axvline(x=0, color='grey', lw=0.7, ls='--')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title, fontsize=13)

    plt.savefig(figname, bbox_inches='tight', pad_inches=0.02, transparent=True)
    plt.close()

# ==== Main script ====

# Define output directory
outdir = 'f5_per_CT_TFBS_CP_cor_SE'

# TFBS comparison types
tfbs_cp_types = [
    'CP_TFBS_nonBlackList_vs_TFMS',
    'CP_TFBS_nonBlackList_overlap_motif_vs_TFMS',
    'CP_TFBS_nonBlackList_NOT_overlap_motif_vs_TFMS'
]

cp_types = ['TFMS_CP', 'TFBS_CP']
tfbs_dir = 'f3_TFBS_CP_heatmap_with_motif_SE/_csv/'

# Process each TFBS comparison type
for tfbs_cp_type in tfbs_cp_types:
    for cp_type in cp_types:
        os.makedirs(f'{outdir}/{cp_type}', exist_ok=True)

        tfbs_cp_s = pd.read_csv(f'{tfbs_dir}/data_fisher_{tfbs_cp_type}_CP_KStest_statistics.csv', index_col=0)
        tfbs_cp_s_z = pd.read_csv(f'{tfbs_dir}/data_fisher_{tfbs_cp_type}_CP_KStest_statistics_zscored.csv', index_col=0)
        tfbs_se_s = pd.read_csv(f'{tfbs_dir}/data_fisher_{tfbs_cp_type}_SE_statistics.csv', index_col=0)

        # Generate scatter plots for each cell type (CT)
        for ct in tfbs_cp_s.columns:
            y = tfbs_cp_s[ct].dropna()
            if len(y) < 3:
                continue

            x = tfbs_cp_s_z[ct].loc[y.index]
            z = np.log2(tfbs_se_s[ct].loc[y.index])

            df = pd.concat([x, y, z], axis=1)
            df.columns = ['Z-scored CP', 'TFBS CP', 'log2 Odds Ratio']

            figname = f'{outdir}/{cp_type}/{tfbs_cp_type}_{ct}.pdf'
            scatter_plot(df, 'Z-scored CP', 'TFBS CP', figname, ct)

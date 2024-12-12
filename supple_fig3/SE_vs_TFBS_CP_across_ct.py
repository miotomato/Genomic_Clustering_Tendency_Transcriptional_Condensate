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

# Plot settings
matplotlib.rcParams['font.size'] = 12
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid': False})
sns.set_style("ticks")

# Function to mark p-values on boxplots
def mark_pvalue(compr_pos, positions, box_vals):
    s, p = stats.ttest_ind(box_vals[compr_pos[0]], box_vals[compr_pos[1]], nan_policy='omit')
    y, h, col = np.percentile(np.append(box_vals[compr_pos[0]], box_vals[compr_pos[1]]), 95) * 0.95, 1.05, 'k'
    y2 = np.percentile(np.append(box_vals[compr_pos[0]], box_vals[compr_pos[1]]), 5) * 0.99
    x1, x2 = positions[compr_pos[0]], positions[compr_pos[1]]
    indicator = 'down' if s > 0 else 'up'
    p_label = '{} \n {:.1e}'.format(indicator, p)
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


# Function to generate scatter plots
def scatter_plot(df, xlabel, ylabel, figname, title):
    x = df[xlabel]
    y = df[ylabel]

    plt.figure(figsize=(2.2, 2.2))
    plt.scatter(x, y, c='k', s=9)

    # Highlight top-ranked points
    color = 'tab:red'
    label_i = 0
    y_max = abs(max(plt.ylim(), key=abs)) * 0.9
    y_step = y_max / 5
    x_max = abs(max(plt.xlim(), key=abs))

    top_y = df['TFBS CP'].sort_values(ascending=False).index[:5]
    for label_index in top_y:
        plt.scatter(x[label_index], y[label_index], c=color, s=30, label=label_index)
        basex = x[label_index]
        basey = y[label_index]
        dx = x_max / 2.2
        y_adjust = y_max
        y_max = y_adjust - y_step
        plt.annotate(label_index, xy=(basex, basey), xytext=(basex + dx, y_adjust), color='k',
                     arrowprops=dict(color='k', lw=0.7, arrowstyle="-"))
        label_i += 1

    # Linear regression line
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    x_sort = np.sort(x)
    plt.plot(x_sort, x_sort * slope + intercept, c='grey', ls='--', lw=0.6)
    plt.text(0.6, 0.06, r'$r = {:.2f}$'.format(r_value), fontsize=11, transform=plt.gca().transAxes)

    plt.axhline(y=0, color='grey', lw=0.7, ls='--')
    plt.axvline(x=0, color='grey', lw=0.7, ls='--')

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title, fontsize=13)

    plt.savefig(figname, bbox_inches='tight', pad_inches=0.02, transparent=True)
    plt.show()
    plt.close()


# Directory for output
outdir = 'f5_per_CT_TFBS_CP_cor_SE'
os.makedirs(outdir, exist_ok=True)

# TFBS data directory
tfbs_dir = 'f3_TFBS_CP_heatmap_with_motif_SE/_csv/'
tfbs_cp_types = ['CP_TFBS_nonBlackList_vs_TFMS']

cp_types = ['TFMS_CP', 'TFBS_CP']

for tfbs_cp_type in tfbs_cp_types:
    tfbs_cp_s = pd.read_csv(f'{tfbs_dir}/data_fisher_{tfbs_cp_type}_CP_KStest_statistics.csv', index_col=0)
    tfbs_cp_s_z = pd.read_csv(f'{tfbs_dir}/data_fisher_{tfbs_cp_type}_CP_KStest_statistics_zscored.csv', index_col=0)
    tfbs_se_s = pd.read_csv(f'{tfbs_dir}/data_fisher_{tfbs_cp_type}_SE_statistics.csv', index_col=0)

    for cp_type in cp_types:
        os.makedirs(f'{outdir}/{cp_type}', exist_ok=True)

        for ct in tfbs_cp_s.columns:
            y = tfbs_cp_s[ct].dropna()
            if len(y) < 3:
                continue

            x = tfbs_cp_s_z[ct].loc[y.index]
            z = np.log2(tfbs_se_s[ct].loc[y.index])

            df = pd.concat([x, y, z], axis=1)
            df.columns = ['Z-scored CP', 'TFBS CP', 'log2 Odds Ratio']
            df['Z-scored CP rank'] = df['Z-scored CP'].rank(ascending=False)
            df['TFBS CP rank'] = df['TFBS CP'].rank(ascending=False)
            df['avg rank'] = (df['Z-scored CP rank'] + df['TFBS CP rank']) / 2

            df.to_csv(f'{outdir}/{cp_type}/{tfbs_cp_type}_{ct}.csv')

            figname = f'{outdir}/{cp_type}/{tfbs_cp_type}_{ct}.pdf'
            scatter_plot(df, 'log2 Odds Ratio', 'TFBS CP', figname, ct)

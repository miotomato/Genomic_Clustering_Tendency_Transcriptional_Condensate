import os
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

# Plot settings
matplotlib.rcParams['font.size'] = 12
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
sns.set(font_scale=1.1)
sns.set_style("whitegrid", {'axes.grid': False})
sns.set_style("ticks")

def mark_pvalue(compr_pos, positions, box_vals):
    s, p = stats.ttest_ind(box_vals[compr_pos[0]], box_vals[compr_pos[1]], nan_policy='omit')
    y, h, col = np.percentile(np.append(box_vals[compr_pos[0]], box_vals[compr_pos[1]]), 95) * 0.95, 1.05, 'k'
    y2 = np.percentile(np.append(box_vals[compr_pos[0]], box_vals[compr_pos[1]]), 5) * 0.99
    x1, x2 = positions[compr_pos[0]], positions[compr_pos[1]]
    indicator = 'down' if s > 0 else 'up'
    p_label = '{} \n {:.1e}'.format(indicator, p)
    if p >= 0.05:
        p_label = 'n.s.'
    plt.plot([x1 * 1.03, x1 * 1.03, x2 * 0.97, x2 * 0.97], [y, y * h, y * h, y], lw=1, c=col)
    plt.text((x1 + x2) * 0.5, y * h, p_label, ha='center', va='bottom', color=col, fontsize=12)

def scatter_plot(df, xlabel, ylabel, figname, title):
    x = df[xlabel]
    y = df[ylabel]
    
    plt.figure(figsize=(2.2, 2.2))
    plt.scatter(x, y, c='k', s=9)
    
    # Linear regression line
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    plt.text(0.6, 0.06, r'$r = {:.2f}$'.format(r_value), fontsize=11, transform=plt.gca().transAxes)
    plt.axhline(y=0, color='grey', lw=0.7, ls='--')
    plt.axvline(x=0, color='grey', lw=0.7, ls='--')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title, fontsize=13)
    
    plt.savefig(figname, bbox_inches='tight', pad_inches=0.02, transparent=True)
    plt.show()
    plt.close()

# Color palette function
def return_colors(pal, half_len, a, b, c=1):
    norm = matplotlib.colors.Normalize(vmin=a * half_len, vmax=b * half_len)
    color_map = matplotlib.cm.ScalarMappable(norm=norm, cmap=pal)
    return [color_map.to_rgba(i) for i in np.arange(c * half_len)]

# Color settings
half_len = 5
pal1 = plt.cm.RdPu_r
pal2 = plt.cm.GnBu
colors1 = return_colors(pal1, half_len, -0.2, 1.0)
colors2 = return_colors(pal2, half_len, -0.2, 1.0)
colors = colors1 + colors2

# Output directory
outdir = 'f6_per_TF_TFBS_CP_cor_SE'
os.makedirs(outdir, exist_ok=True)

# TFBS data directory
tfbs_dir = 'f3_TFBS_CP_heatmap_with_motif_SE/_csv/'
tfbs_cp_types = ['CP_TFBS_nonBlackList_vs_TFMS']
cp_types = ['TFMS_CP', 'TFBS_CP']

# Main loop to process and generate scatter plots
for tfbs_cp_type in tfbs_cp_types:
    tfbs_cp_s = pd.read_csv(f'{tfbs_dir}/data_fisher_{tfbs_cp_type}_CP_KStest_statistics.csv', index_col=0)
    tfbs_cp_s_z = pd.read_csv(f'{tfbs_dir}/data_fisher_{tfbs_cp_type}_CP_KStest_statistics_zscored.csv', index_col=0)
    tfbs_se_s = pd.read_csv(f'{tfbs_dir}/data_fisher_{tfbs_cp_type}_SE_statistics.csv', index_col=0)

    for cp_type in cp_types:
        os.makedirs(f'{outdir}/{cp_type}', exist_ok=True)
        
        for factor in tfbs_cp_s.index:
            x = tfbs_cp_s_z.loc[factor].dropna()
            if len(x) < 5:
                continue

            y = tfbs_cp_s.loc[factor].loc[x.index]
            z = np.log2(tfbs_se_s.loc[factor].loc[x.index])
            
            df = pd.concat([x, y, z], axis=1)
            df.columns = ['Z-scored CP', 'TFBS CP', 'log2 Odds Ratio']
            
            figname = f'{outdir}/{cp_type}/{tfbs_cp_type}_{factor}.pdf'
            scatter_plot(df, 'log2 Odds Ratio', 'TFBS CP', figname, factor)

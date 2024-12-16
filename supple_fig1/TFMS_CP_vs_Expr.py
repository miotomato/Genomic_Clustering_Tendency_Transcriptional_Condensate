import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from matplotlib import gridspec

# == Matplotlib and Seaborn Configuration ==
matplotlib.rcParams['font.size'] = 18
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams["mathtext.rm"] = "Arial"
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid': False})
sns.set_style("ticks")

# == Function to Add P-Value Annotations ==
def mark_pvalue(compr_pos, positions, box_vals):
    s, p = stats.ttest_ind(box_vals[compr_pos[0]], box_vals[compr_pos[1]], nan_policy='omit')
    y, h, col = np.percentile(np.append(box_vals[compr_pos[0]], box_vals[compr_pos[1]]), 96) * 0.99, 1.05, 'k'
    x1, x2 = positions[compr_pos[0]], positions[compr_pos[1]]
    p_label = '{:.1e}'.format(p)
    if p_label[-2] == '0':
        p_label = p_label[:-2] + p_label[-1]
    if p >= 0.05:
        p_label = 'n.s.'
    plt.plot([x1 * 1.03, x1 * 1.03, x2 * 0.97, x2 * 0.97], [y, y * h, y * h, y], lw=1, c=col)
    plt.text((x1 + x2) * 0.5, y * h, p_label, ha='center', va='bottom', color=col, fontsize=15)

# == Output Directory ==
outdir = 'f12_TFMSCP_cor_expr'
os.makedirs(outdir, exist_ok=True)

# == Load Data ==
df = pd.read_csv('data/TFMS_CP_SE_enrich.csv', index_col=0)
df['TFMSCP'] = df['log10-dis ks_2samp-s signed']
df.index = [i.split('_')[0] for i in df.index]

expr_df = pd.read_csv('data/TFexpr_var_mean.csv', index_col=0)
shared_tf = df.index.intersection(expr_df.index)
outdf = pd.concat([df[['TFMSCP']].loc[shared_tf], expr_df.loc[shared_tf]], axis=1)
outdf = outdf.sort_values(by='TFMSCP', ascending=False)
outdf.to_csv(f'{outdir}/TFMSCP_cor_expr.csv')

# == Plot Scatter Plots ==
colnames = ['FPKM_log_var', 'TPM_log_var']

for colname in colnames:
    plt.figure(figsize=(2.5, 2.5))
    x = outdf['TFMSCP']
    y = outdf[colname]
    plt.scatter(x, y, s=5, color='k')

    # Linear Regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    x_sort = np.sort(x)
    plt.plot(x_sort, x_sort * slope + intercept, c='grey', ls='--', lw=0.6)
    plt.text(0.6, 0.8, f'$R^2$ = {r_value**2:.2f}', fontsize=12, transform=plt.gca().transAxes)

    plt.axvline(x=0, color='k', lw=1.2, ls='--')
    plt.xlabel('TFMS CP')
    plt.ylabel(colname)
    plt.savefig(f'{outdir}/scatter_TFMSCP_cor_Expr_{colname}.pdf', bbox_inches='tight', pad_inches=0.1, transparent=True, dpi=600)
    plt.close()

# == Plot Boxplots Based on TFMSCP Threshold ==
plt.figure(figsize=(7, 3))
gs = gridspec.GridSpec(2, 1, height_ratios=[1, 2], wspace=0.1, hspace=0.1)

# == Top Panel: TFMS CP Values ==
ax = plt.subplot(gs[0, 0])
vals = outdf['TFMSCP']
ax.scatter(np.arange(len(vals)), vals, s=3, c='k')
ax.axhline(y=0, color='k', lw=1.2, ls='--')
ax.xaxis.tick_top()
ax.set_ylabel('TFMS CP', ha='center')

# == Bottom Panel: Boxplot for Expression Values ==
ax = plt.subplot(gs[1, 0])
pos_len = 20
vals = outdf['FPKM_log_var']
cutoff = vals.median()
positions = np.arange(pos_len)
box_step = int(len(vals) / len(positions))
box_vals = [vals[i * box_step:(i + 1) * box_step] for i in positions]

# Plot Boxplots
ax.boxplot(box_vals, positions=positions, widths=0.6, patch_artist=True,
           boxprops=dict(color='k', facecolor='w', fill=None, lw=1),
           medianprops=dict(color='grey'), showfliers=False)

ax.axhline(y=cutoff, color='k', lw=1.2, ls='--')
ax.set_xticks([])
ax.set_ylabel('FPKM Log Variance')
plt.savefig(f'{outdir}/box_TFMS_CP_cor_Expr_FPKM_log_var.pdf', bbox_inches='tight', pad_inches=0.1, transparent=True, dpi=600)
plt.close()

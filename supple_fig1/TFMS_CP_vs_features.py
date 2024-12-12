import os
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import gaussian_kde

# == Plotting Configuration
plt.rcParams['font.size'] = 14
plt.rcParams["font.sans-serif"] = ["Arial"]
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid': False})
sns.set_style("ticks")

# == Function Definitions

def scatter_plot(value_dic, xlabel, ylabel, key_pairs, mark=False):
    """Generate and save a scatter plot with density and optional outlier marking."""
    x = value_dic[key_pairs[0]]
    y = value_dic[key_pairs[1]]
    figname = f"{outdir}/{key_pairs[0]}_vs_{key_pairs[1]}_rmfitline.pdf"

    # Create figure and axes
    fig, ax = plt.subplots(figsize=(2.6, 2.6))
    
    # Calculate density
    xy = np.vstack([x, y])
    z = gaussian_kde(xy)(xy)
    z[np.isnan(z)] = 0.0
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]
    
    # Scatter plot with density
    ax.scatter(x, y, c=z, cmap=plt.cm.GnBu_r, s=3, marker='o')
    
    # Linear regression
    slope, intercept, r_value, _, _ = stats.linregress(x, y)
    x_sort = np.sort(x)
    ax.plot(x_sort, x_sort * slope + intercept, c='grey', ls='--', lw=0.6)
    ax.text(0.65, 0.18, f'$R^2$ = {r_value**2:.2f}', fontsize=11, transform=ax.transAxes)

    # Mark outliers if specified
    if mark:
        y_new = x * slope + intercept
        marked_indexes = (np.abs(y - y_new)).sort_values(ascending=False)[:6].index
        for marked_index in marked_indexes:
            ax.text(x[marked_index], y[marked_index], marked_index.split('_')[0], c='r', fontsize=6)

    # Plot aesthetics
    ax.axhline(y=0, color='k', lw=1.2, ls='--')
    ax.axvline(x=0, color='k', lw=1.2, ls='--')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    # Save the figure
    plt.savefig(figname, bbox_inches='tight', pad_inches=0.02, transparent=True)
    plt.close()

# == Output Directory
outdir = 'f1_TFMS_CP_cor_length_number'
os.makedirs(outdir, exist_ok=True)
os.makedirs('data', exist_ok=True)

# == Load Data
df1 = pd.read_csv('../f1_TFMS_TFBS_CP/CP_TFMS_vs_random_nonBlackList/data_CP_TFMS_vs_random.csv', index_col=0)
df1 = df1[['#TFMS', 'len-of-TFMS']]

df2 = pd.read_csv('../f1_TFMS_TFBS_CP/CP_TFMS_vs_random_nonBlackList/data_TFMS_enrich_at_SE.csv', index_col=0)

# Merge DataFrames
df = pd.concat([df1, df2], axis=1, join='inner')
df.to_csv('data/TFMS_CP_SE_enrich.csv')

# Sort by KS statistics
df = df.sort_values(by='log10-dis ks_2samp-s signed', ascending=False)

# Calculate genome total
df['genome-total'] = df['#TFMS'] * df['len-of-TFMS']

# Create value dictionary for plotting
value_dic = {
    'Motif_num': df['#TFMS'],
    'Motif_length': df['len-of-TFMS'],
    'Genome_total': df['genome-total'],
    'KS_statistics': df['log10-dis ks_2samp-s signed']
}

# == Generate Scatter Plots

# Plot 1: TFMS count vs Motif length
scatter_plot(value_dic, 'TFMS count', 'Motif length', ['Motif_num', 'Motif_length'])

# Plot 2: TFMS CP vs TFMS count
scatter_plot(value_dic, 'TFMS CP', 'TFMS count', ['KS_statistics', 'Motif_num'])

# Plot 3: TFMS CP vs Fraction of genome
scatter_plot(value_dic, 'TFMS CP', 'Fraction of genome', ['KS_statistics', 'Genome_total'])

# Plot 4: TFMS CP vs Motif length
scatter_plot(value_dic, 'TFMS CP', 'Motif length', ['KS_statistics', 'Motif_length'])

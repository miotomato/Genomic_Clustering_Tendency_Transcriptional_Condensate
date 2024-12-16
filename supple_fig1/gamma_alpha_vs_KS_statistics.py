import os
import numpy as np
import pandas as pd
import scipy
import scipy.stats
import matplotlib.pyplot as plt
import seaborn as sns

# == Configure Plotting Parameters
plt.rcParams['font.size'] = 16
plt.rcParams["font.sans-serif"] = ["Arial"]
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid': False})
sns.set_style("ticks")

# == Chromosome List
hg38_chroms = [
    'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
    'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
    'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY'
]

# == Function Definitions

def scatter_plot(x, y, xlabel, ylabel, figname):
    """Generate and save a scatter plot with linear regression line."""
    plt.figure(figsize=(2.6, 2.6))
    plt.scatter(x, y, c='k', s=3)

    # Perform linear regression
    slope, intercept, r_value, _, _ = scipy.stats.linregress(x, y)
    x_sort = np.sort(x)
    plt.plot(x_sort, x_sort * slope + intercept, c='grey', ls='--', lw=0.9)
    plt.text(0.55, 0.85, f'$R$ = {r_value:.2f}', fontsize=12, transform=plt.gca().transAxes)

    # Plot aesthetics
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.axhline(y=0, color='k', lw=0.5, ls='--')
    plt.axvline(x=0, color='k', lw=0.5, ls='--')

    # Save the figure
    plt.savefig(figname, bbox_inches='tight', pad_inches=0.02, transparent=True)
    plt.show()
    plt.close()

# == Output Directory
outdir = 'f2_gamma_alpha_vs_KS'
os.makedirs(outdir, exist_ok=True)

# == Load Data
df1 = pd.read_csv('f1_TFMS_gamma_alpha/TFMS_gamma_alpha.csv', index_col=0)
df1 = df1[['alpha', 'scale']]
df2 = pd.read_csv('data/TFMS_CP_SE_enrich.csv', index_col=0)
df = pd.concat([df1, df2], axis=1, join='inner')
df = df[['#TFMS', 'alpha', 'scale', 'log10-dis ks_2samp-s signed']]

# == Calculate Gamma Distribution Metrics and KS Scores
for index in df.index:
    print(index)

    # Retrieve parameters
    alpha = df.loc[index, 'alpha']
    scale = df.loc[index, 'scale']
    size = int(df.loc[index, '#TFMS'])

    # Calculate f(alpha) using the gamma distribution
    g = scipy.special.gamma(alpha)
    g = scale * g ** (1 / (alpha - 1))
    f_alpha = scipy.stats.gamma.cdf(g, a=alpha, scale=scale) - (1 - np.exp(-g / scale))
    df.loc[index, 'f(alpha)'] = f_alpha

    # Generate random samples for KS test
    data1 = scipy.stats.gamma.rvs(a=alpha, scale=scale, size=size)
    data2 = scipy.stats.expon.rvs(scale=scale, size=size)

    # KS test on raw data
    s, _ = scipy.stats.ks_2samp(data1, data2)
    df.loc[index, 'ks from distribution'] = s

    # KS test on log-transformed data
    s, _ = scipy.stats.ks_2samp(np.log10(data1), np.log10(data2))
    df.loc[index, 'ks from log10 distribution'] = s

# == Save the Processed Data
df.to_csv(f'{outdir}/data_gamma_fit.csv')

# == Generate Scatter Plots

# Scatter of alpha vs KS score
x = df['alpha']
y = df['log10-dis ks_2samp-s signed']
xlabel, ylabel = r'Gamma $\alpha$', 'KS score'
figname = f'{outdir}/gamma_alpha_vs_KS.pdf'
scatter_plot(x, y, xlabel, ylabel, figname)

# Scatter of log2(alpha) vs KS score
x = np.log2(df['alpha'])
y = df['log10-dis ks_2samp-s signed']
xlabel, ylabel = r'log2 Gamma $\alpha$', 'KS score'
figname = f'{outdir}/gamma_log_alpha_vs_KS.pdf'
scatter_plot(x, y, xlabel, ylabel, figname)

# Scatter of f(alpha) vs KS score
x = df['f(alpha)']
y = df['log10-dis ks_2samp-s signed']
xlabel, ylabel = r'f($\alpha$)', 'KS score'
figname = f'{outdir}/gamma_Fa_vs_KS.pdf'
scatter_plot(x, y, xlabel, ylabel, figname)

# Scatter of alpha vs f(alpha)
x = df['alpha']
y = df['f(alpha)']
xlabel, ylabel = r'Gamma $\alpha$', r'f($\alpha$)'
figname = f'{outdir}/gamma_alpha_vs_Fa.pdf'
scatter_plot(x, y, xlabel, ylabel, figname)

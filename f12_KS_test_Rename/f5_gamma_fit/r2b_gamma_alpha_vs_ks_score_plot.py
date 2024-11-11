import os,sys,argparse,glob,re,bisect
import numpy as np
import pandas as pd
from collections import Counter
import operator
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
#matplotlib.rcParams['agg.path.chunksize'] = 10000
matplotlib.rcParams['font.size']=16
import seaborn as sns
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid' : False})
import scipy
import scipy.optimize
sns.set_style("ticks")
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
from scipy import stats
import time


hg38_chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',
             'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',
             'chr18','chr19','chr20','chr21','chr22','chrX','chrY']



def scatter_plot(x,y,xlabel,ylabel,figname):
    
    plt.figure(figsize=(2.5,2.5))
    # xy = np.vstack([x,y])
    # z = gaussian_kde(xy)(xy)
    # z[np.where(np.isnan(z))] = 0.0
    # idx = z.argsort()
    # x, y, z = x[idx], y[idx], z[idx]
    # g = plt.scatter(x,y,c=z,cmap = plt.cm.GnBu_r,s=3,marker='o')
    plt.scatter(x,y,c='k',s=3)
    # ==== linear regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)       
    x_sort = np.sort(x)
    plt.plot(x_sort,x_sort*slope+intercept,c = 'grey',ls='--',lw=.9)
    plt.text(.55,.85,'$R$ = {:.2f}'.format(r_value),fontsize=12,transform=plt.axes().transAxes)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.axhline(y=0,color='k',lw=.5,ls='--')
    plt.axvline(x=0,color='k',lw=.5,ls='--')
    plt.savefig('{}'.format(figname),bbox_inches='tight',pad_inches=0.02,transparent=True)
    plt.show()
    plt.close()




outdir = 'f2_gamma_alpha_vs_KS'

df = pd.read_csv('{}/data_gamma_fit.csv'.format(outdir),index_col=0)
# df = df.dropna()


# scatter of alphs vs. KS score
x = df['alpha']
y = df['log10-dis ks_2samp-s signed']
xlabel,ylabel = r'Gamma $k$', 'TFMS CP'
figname = '{}/gamma_alpha_vs_KS.pdf'.format(outdir)
scatter_plot(x,y,xlabel,ylabel,figname)


# scatter of log alphs vs. KS score
x = np.log2(df['alpha'])
y = df['log10-dis ks_2samp-s signed']
xlabel,ylabel = r'log2 Gamma $k$', 'TFMS CP'
figname = '{}/gamma_log_alpha_vs_KS.pdf'.format(outdir)
scatter_plot(x,y,xlabel,ylabel,figname)


# scatter of function alpha vs KS score
x = df['f(alpha)']
y = df['log10-dis ks_2samp-s signed']
xlabel,ylabel = r'f($\alpha$)', 'KS score'
figname = '{}/gamma_Fa_vs_KS.pdf'.format(outdir)
scatter_plot(x,y,xlabel,ylabel,figname)


# scatter of function alpha vs alpha
x = df['alpha']
y = df['f(alpha)']
xlabel,ylabel = r'Gamma $\alpha$', r'f($\alpha$)'
figname = '{}/gamma_alpha_vs_Fa.pdf'.format(outdir)
scatter_plot(x,y,xlabel,ylabel,figname)




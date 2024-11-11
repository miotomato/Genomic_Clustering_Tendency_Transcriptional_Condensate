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



def scatter_plot(df,columns,outdir):
    
    xlabel,ylabel = columns[0], columns[1]
    plt.figure(figsize=(2.6,2.6))
    # xy = np.vstack([x,y])
    # z = gaussian_kde(xy)(xy)
    # z[np.where(np.isnan(z))] = 0.0
    # idx = z.argsort()
    # x, y, z = x[idx], y[idx], z[idx]
    # g = plt.scatter(x,y,c=z,cmap = plt.cm.GnBu_r,s=3,marker='o')
    x = df[xlabel]
    y = df[ylabel]
    plt.scatter(x,y,c='k',s=3)
    # ==== linear regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)       
    x_sort = np.sort(x)
    plt.plot(x_sort,x_sort*slope+intercept,c = 'grey',ls='--',lw=.9)
    plt.text(.55,1.15,'$r$ = {:.2f}'.format(r_value),fontsize=12,transform=plt.axes().transAxes)
    plt.text(.55,1.05,'$p$ = {:.1e}'.format(p_value),fontsize=12,transform=plt.axes().transAxes)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    # plt.axhline(y=0,color='k',lw=.5,ls='--')
    # plt.axvline(x=0,color='k',lw=.5,ls='--')
    figname = '{}/fig_{}_vs_{}.pdf'.format(outdir,xlabel,ylabel)
    plt.savefig('{}'.format(figname),bbox_inches='tight',pad_inches=0.02,transparent=True)
    plt.show()
    plt.close()




outdir = 'f2_gamma_alpha_vs_entropy'
os.makedirs(outdir,exist_ok=True)

df1 = pd.read_csv('f1_TFMS_entropy/TFMS_entropy.csv',index_col=0)
df2 = pd.read_csv('../f5_gamma_fit/f2_gamma_alpha_vs_KS/data_gamma_fit.csv',index_col=0)
df = pd.concat([df1,df2],axis=1,join='inner')

df = df[['entropy_sum', 'entropy_avg', 'motif_len', '#TFMS', 'alpha', 
        'log10-dis ks_2samp-s signed']]
df.columns = ['entropy_sum', 'entropy_avg', 'motif_len', 'TFMS_num', 'alpha', 'TFMS_CP']


# scatter for correlation    
columns = ['motif_len', 'entropy_sum']
scatter_plot(df,columns,outdir)

columns = ['motif_len', 'entropy_avg']
scatter_plot(df,columns,outdir)

columns = ['motif_len', 'alpha']
scatter_plot(df,columns,outdir)

columns = ['motif_len', 'TFMS_CP']
scatter_plot(df,columns,outdir)

columns = ['entropy_sum', 'alpha']
scatter_plot(df,columns,outdir)

columns = ['entropy_avg', 'alpha']
scatter_plot(df,columns,outdir)

columns = ['entropy_sum', 'TFMS_CP']
scatter_plot(df,columns,outdir)

columns = ['entropy_avg', 'TFMS_CP']
scatter_plot(df,columns,outdir)
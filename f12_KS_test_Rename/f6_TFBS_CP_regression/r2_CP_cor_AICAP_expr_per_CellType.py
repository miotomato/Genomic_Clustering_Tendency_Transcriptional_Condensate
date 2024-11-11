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
from scipy.stats import gamma


hg38_chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',
             'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',
             'chr18','chr19','chr20','chr21','chr22','chrX','chrY']



def scatter_plot(df,cellType,xlabel,ylabel,outdir):
    
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
    plt.text(.85,1.05,'$R$ = {:.2f}'.format(r_value),fontsize=12,transform=plt.axes().transAxes)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(cellType)
    # plt.axhline(y=0,color='k',lw=.5,ls='--')
    # plt.axvline(x=0,color='k',lw=.5,ls='--')
    figname = '{}/{}_vs_{}_{}.pdf'.format(outdir,xlabel,ylabel,cellType)
    plt.savefig('{}'.format(figname),bbox_inches='tight',pad_inches=0.02,transparent=True)
    plt.show()
    plt.close()




indir = 'f1_CP_AICAP_expr_correlation'
outdir = 'f2_CP_cor_AICAP_expr_per_CellType'
os.makedirs(outdir,exist_ok=True)


df_all = pd.read_csv('{}/data_all.csv'.format(indir),index_col=0)
cellTypes = list(set([ii.split('_')[0] for ii in df_all.index]))
columns = ['TFMS_CP', 'TFMS_alpha', 'TFBS_CP', 'TFBS_ZCP', 'TFBS_CP_avgRank',
       'TFBS_SE', 'Expr_TPM', 'Expr_log2_TPM_Plus1', 'negLog2_IDR']


for cellType in cellTypes[:]:
    df = pd.read_csv('{}/data_{}.csv'.format(indir,cellType),index_col=0)

    xlabel = 'TFBS_CP'
    ylabel = 'TFMS_CP'
    scatter_plot(df,cellType,xlabel,ylabel,outdir)

    xlabel = 'TFBS_CP'
    ylabel = 'TFBS_SE'
    scatter_plot(df,cellType,xlabel,ylabel,outdir)
    
    xlabel = 'TFBS_CP'
    ylabel = 'Expr_log2_TPM_Plus1'
    scatter_plot(df,cellType,xlabel,ylabel,outdir)

    xlabel = 'TFBS_CP'
    ylabel = 'negLog2_IDR'
    scatter_plot(df,cellType,xlabel,ylabel,outdir)


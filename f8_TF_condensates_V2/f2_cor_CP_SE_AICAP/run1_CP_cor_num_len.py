import os,sys,argparse
# import fileinput,time
import glob
import re,bisect
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=14
import seaborn as sns
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid' : False})
import scipy
import scipy.optimize
sns.set_style("ticks")
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams["mathtext.rm"] = "Arial"
from scipy.interpolate import interpn
from scipy.stats import gaussian_kde



def scatter_plot(x,y,xlabel,ylabel,figname):

    # ==== scatter plot with density
    plt.figure(figsize=(2.6,2.6))
    xy = np.vstack([x,y])
    z = gaussian_kde(xy)(xy)
    z[np.where(np.isnan(z))] = 0.0
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]
    g = plt.scatter(x,y,c=z,cmap = plt.cm.GnBu_r,s=3,marker='o')
    # plt.scatter(x,y,c='k',s=1)
    # ==== linear regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)       
    x_sort = np.sort(x)
    plt.plot(x_sort,x_sort*slope+intercept,c = 'grey',ls='--',lw=.6)
    plt.text(.55,1.02,'$R^2$ = {:.2f}'.format(r_value**2),fontsize=13,transform=plt.axes().transAxes)

    plt.axhline(y=0,color='k',lw=1.2,ls='--')
    plt.axvline(x=0,color='k',lw=1.2,ls='--')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig('{}'.format(figname),bbox_inches='tight',pad_inches=0.02,transparent=True)
    plt.show()
    plt.close()
    


# == check the TF motifs
outdir = 'f1_TFMS_CP_cor_length_number'
os.makedirs(outdir,exist_ok=True)


len_df = pd.read_csv('data/data_motif_Length.csv',index_col=0)
num_df = pd.read_csv('../f1_TF_cluster_potential/f2_TFMS_TFBS_CP/CP_TFMS_vs_random/data_CP_TFMS_vs_random.csv',index_col=0)
df = pd.concat([len_df,num_df],axis=1)
df.to_csv('data/data_TFMS_CP_motif_Num_Len.csv')
df = df[df['log10 distance t-stats']<500]

# == plot 
x = df['motif_Length']
y = df['#TFMS']
xlabel = 'Length of motif'
ylabel = 'Number of TFMS'
figname = outdir+os.sep+'motif_Len_vs_Num.pdf'
scatter_plot(x,y,xlabel,ylabel,figname)


# == plot 
x = df['log10 distance t-stats']
y = df['motif_Length']
xlabel = 'Cluster potential'
ylabel = 'Length of motif'
figname = outdir+os.sep+'motif_CP_vs_Len.pdf'
scatter_plot(x,y,xlabel,ylabel,figname)


# == plot 
x = df['log10 distance t-stats']
# x = np.log2(x.clip(1))
y = df['#TFMS']
# y = np.log10(y)
xlabel = 'Cluster potential'
ylabel = 'Number of TFMS'
figname = outdir+os.sep+'motif_CP_vs_Num.pdf'
scatter_plot(x,y,xlabel,ylabel,figname)



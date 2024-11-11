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


# == check the TF motifs
outdir = 'f1_IDR_vs_cluster_potential'
os.makedirs(outdir,exist_ok=True)

# read data of cluster potential
infile = '../f1_TF_motif_cluster_cor_SE/f4_cluster_enrichment_at_SE/f3_SE_OR_vs_cluster_potential_figs/data.csv'
df = pd.read_csv(infile,index_col=0)
# df.index = [i.split('_')[0] for i in df.index]
cutoff_OR2 = df['cutoff of OR>2']
OR_100bp = df['OR of 100bp']
OR_200bp = df['OR of 200bp']
log10_stats = df['log10 distance t-stats']
# raw_stats = df['distance t-stats']
or_at_SE = df.loc[df.index,'fisher_exact_s']

# read TF AICAP
idrfile = 'data/13059_2021_2456_MOESM2_ESM.xlsx'
idrdf = pd.read_excel(idrfile,sheet_name='condition 2(1,6-HD-2)',index_col=0)

# get shared TF/index
shared_tfs = df.index.intersection(idrdf.index)

# == plot 
plt.figure(figsize=(2.6,2.6))
x = log10_stats.loc[shared_tfs]
y = idrdf.loc[shared_tfs].AICAP
y = np.log2(y)
plt.scatter(x,y,c='k',s=11)
plt.axhline(y=0,color='k',lw=1.2,ls='--')
plt.axvline(x=0,color='k',lw=1.2,ls='--')
plt.xlabel('Cluster potential')
plt.ylabel('log2 AICAP')
# plt.ylim([0,2])
# plt.xlim([-330,680])
plt.savefig('{}/IDR_vs_cluster_potential.pdf'.format(outdir),bbox_inches='tight',pad_inches=0.02,transparent=True)
plt.show()
plt.close()


# == plot 
plt.figure(figsize=(2.6,2.6))
x = or_at_SE.loc[shared_tfs]
y = idrdf.loc[shared_tfs].AICAP
y = np.log2(y)
plt.scatter(x,y,c='k',s=11)
plt.axhline(y=0,color='k',lw=1.2,ls='--')
plt.axvline(x=1,color='k',lw=1.2,ls='--')
plt.xlabel('OR at SE')
plt.ylabel('log2 AICAP')
# plt.ylim([0,2])
# plt.xlim([-330,680])
plt.savefig('{}/IDR_vs_OR_st_SE.pdf'.format(outdir),bbox_inches='tight',pad_inches=0.02,transparent=True)
plt.show()
plt.close()


# == plot 
from matplotlib import gridspec
plt.figure(figsize=(3.2,2.6))
width_ratios = [4,.5]
gs = gridspec.GridSpec(2,2,width_ratios=width_ratios,wspace=.2,hspace=.2)
vmin=-2; vmax=2
# ==== scatter plot
ax = plt.subplot(gs[:,0])   
x = log10_stats.loc[shared_tfs]
y = or_at_SE.loc[shared_tfs]
z = idrdf.loc[shared_tfs].AICAP
z = np.log2(z).clip(vmin,vmax)
# plt.scatter(x,y,c='k',s=11)
ax.scatter(x,y,c=z,cmap=plt.cm.PiYG,s=11)
ax.axhline(y=1,color='k',lw=1.2,ls='--')
ax.axvline(x=0,color='k',lw=1.2,ls='--')
ax.set_xlabel('Cluster potential')
ax.set_ylabel('Odds ratio at SE')
# plt.ylim([0,2])
# plt.xlim([-330,680])
# ==== color bar
ax = plt.subplot(gs[0,1]) 
norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
cbar = matplotlib.colorbar.ColorbarBase(ax,cmap=plt.cm.PiYG,norm=norm,orientation='vertical')
# cbar.set_ticks([-1,0,2])
cbar.set_label('log2 AICAP',labelpad=-55,fontsize=13)

plt.savefig('{}/OR_at_SE_vs_cluster_potential.pdf'.format(outdir),bbox_inches='tight',pad_inches=0.02,transparent=True)
plt.show()
plt.close()


# save the data
outdf = pd.concat([df,idrdf[['AICAP']]],axis=1,join='inner')
# outdf.to_csv('{}/data_add_AICAP.csv'.format(outdir))


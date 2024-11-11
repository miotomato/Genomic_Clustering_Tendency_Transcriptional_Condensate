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
matplotlib.rcParams['font.size']=11
import seaborn as sns
sns.set(font_scale=1.0)
sns.set_style("whitegrid", {'axes.grid' : False})
import scipy
import scipy.optimize
sns.set_style("ticks")
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
from matplotlib import gridspec


# == check the TF motifs
outdir = 'f2_fisher_stats_figs'
os.makedirs(outdir,exist_ok=True)
indir = 'f1_cluster_enrichment_at_SE'
infiles = glob.glob('f1_cluster_enrichment_at_SE/*csv')
infiles.sort()

# read the data
df_stats, df_pvalue = pd.DataFrame(),pd.DataFrame()
for infile in infiles[:]:
    basename = os.path.basename(infile).split('.csv')[0]
    if basename=='all_hg38_SE':
        continue
    df_tmp = pd.read_csv(infile,index_col=0)
    s = df_tmp[['fisher_exact_s']].rename(columns = {'fisher_exact_s':basename})
    p = df_tmp[['fisher_exact_p']].rename(columns = {'fisher_exact_p':basename})
    df_stats = pd.concat([df_stats,s],axis=1)
    df_pvalue = pd.concat([df_pvalue,p],axis=1)
df_stats.to_csv('{}/data_fisher_s.csv'.format(outdir))
df_pvalue.to_csv('{}/data_fisher_p.csv'.format(outdir))

# emrichment at all merged SE
df_all = pd.read_csv('f1_cluster_enrichment_at_SE/all_hg38_SE.csv',index_col=0)
df_all = df_all[['fisher_exact_s']]


# == plot the heatmap
plt.figure(figsize=(2.6,5))
width_ratios = [4,.4]
gs = gridspec.GridSpec(1,2,width_ratios=width_ratios,wspace=0.1,)

# SE from each of the 86 samples
ax = plt.subplot(gs[0,0])    
ranked_row = df_stats.sum(axis=1).sort_values(ascending=False).index
df = df_stats.loc[ranked_row]
# df = np.log2(df)
g = sns.heatmap(df,cmap=plt.cm.bwr,cbar=False\
                ,vmax=2,vmin=0\
                ,ax=ax\
                ,xticklabels=False,yticklabels=False\
                ,cbar_kws={"orientation": "vertical","use_gridspec":False},rasterized=True)
ax.set_xlabel('SE of {} samples '.format(df.shape[1]),ha='center',fontsize=12)
# ax.set_ylabel('{} TFs'.format(df.shape[0]))
ax.set_ylabel(' '.format(df.shape[0]))
# ax.axhline(y=10,color='k',lw=1.2,ls='--')
# == across all SE
ax = plt.subplot(gs[0,1])       
g = sns.heatmap(df_all.loc[ranked_row],cmap=plt.cm.bwr,cbar=True\
                ,vmax=2 ,vmin=0\
                ,ax=ax\
                ,xticklabels=False,yticklabels=False\
                ,cbar_kws={"orientation": "vertical","use_gridspec":False})
ax.set_xlabel('  Union SE',ha='center',fontsize=12)
ax.set_ylabel(' '.format(df.shape[0]),fontsize=12)
cbar = g.collections[0].colorbar
cbar.ax.set_position([.97,0.55,.8,.3]) 
cbar.set_ticks([0,1,2])
cbar.set_label('Odds Ratio',labelpad=-33,fontsize=11)
# ax.yaxis.set_label_position('left')
plt.savefig('{}/cluster_TF_OR_at_SE.pdf'.format(outdir),bbox_inches='tight',pad_inches=0.02,transparent=True)
plt.show()
plt.close()
df_all.loc[ranked_row].to_csv('{}/data_at_union_SE.csv'.format(outdir))








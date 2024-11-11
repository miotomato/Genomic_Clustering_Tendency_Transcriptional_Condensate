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
sns.set(font_scale=1.1)
sns.set_style("whitegrid", {'axes.grid' : False})
import scipy
import scipy.optimize
sns.set_style("ticks")
matplotlib.rcParams["font.sans-serif"] = ["Arial"]

from scipy.interpolate import interpn
from scipy.stats import gaussian_kde


# == check the TF motifs
outdir = 'f4_stats_figs'
os.makedirs(outdir,exist_ok=True)

infile = 'f3_closest_distribution_stats//distance_to_closest_site_vs_control_stats.csv'
df = pd.read_csv(infile,index_col=0)
cutoff_OR2 = df['cutoff of OR>2']
OR_100bp = df['OR of 100bp']
OR_200bp = df['OR of 200bp']
log10_stats = df['log10 distance t-stats']
raw_stats = df['distance t-stats']


# == plot the distribution of true/control data
plt.figure(figsize=(2.6,2.6))
sns.distplot(cutoff_OR2[cutoff_OR2<=500],kde=False)
plt.ylabel('Number of TFs')
plt.xlabel('Cutoff of Odds Ratio>2')
plt.savefig('{}/hist_OddsRatio.pdf'.format(outdir),bbox_inches='tight',pad_inches=0.02,transparent=True)
plt.show()
plt.close()


# == plot the distribution of true/control data
plt.figure(figsize=(2.6,2.6))
# x,y = log10_stats,cutoff_OR2
# xy = np.vstack([x,y])
# z = gaussian_kde(xy)(xy)
# z[np.where(np.isnan(z))] = 0.0
# idx = z.argsort()
# x, y, z = x[idx], y[idx], z[idx]
# g = plt.scatter(x,y,c=z,cmap = plt.cm.GnBu_r,s=3,marker='o')
plt.scatter(cutoff_OR2,log10_stats,color='k',s=5)
# write_out_index = np.append(y[y>5].index, x[x>500].index)
# for ii in write_out_index:
#     xi = x[ii]
#     yi = y[ii]
#     plt.text(xi, yi, ii.split('_')[0],fontsize=10)               
plt.axhline(y=1,color='k',lw=1.2,ls='--')
plt.axvline(x=0,color='k',lw=1.2,ls='--')
plt.ylabel('Cluster potential')
plt.xlabel('Cutoff of OR>2')
plt.savefig('{}/stats_vs_OddsRatio.pdf'.format(outdir),bbox_inches='tight',pad_inches=0.02,transparent=True)
plt.show()
plt.close()


# plot the rank
plt.figure(figsize=(2.6,2.6))
vals = log10_stats.sort_values(ascending=False)
g = plt.scatter(np.arange(len(vals)),vals,s=3,c='k')
plt.axhline(y=0,color='k',lw=1.2,ls='--')
plt.xlabel('Rank')
plt.ylabel('Cluster potential')
plt.savefig('{}/stats_rank.pdf'.format(outdir),bbox_inches='tight',pad_inches=0.02,transparent=True)
plt.show()
plt.close()

# save with new index
df.index = [i.split('_')[0] for i in df.index]
df = df.sort_values(by='log10 distance t-stats',ascending=False)
df.to_csv(outdir+os.sep+'_new_index_stats_data.csv')


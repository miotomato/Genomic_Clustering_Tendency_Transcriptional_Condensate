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
outdir = 'f3_SE_OR_vs_cluster_potential_figs'
os.makedirs(outdir,exist_ok=True)

# read data of cluster potential
infile = '../f2_TF_fimo_jarspar_cluster/f3_closest_distribution_stats/distance_to_closest_site_vs_control_stats.csv'
df = pd.read_csv(infile,index_col=0)
df.index = [i.split('_')[0] for i in df.index]
cutoff_OR2 = df['cutoff of OR>2']
OR_100bp = df['OR of 100bp']
OR_200bp = df['OR of 200bp']
log10_stats = df['log10 distance t-stats']
raw_stats = df['distance t-stats']

# emrichment/OR at all merged SE
df_OR = pd.read_csv('f1_cluster_enrichment_at_SE/all_hg38_SE.csv',index_col=0)
or_at_SE = df_OR.loc[df.index,'fisher_exact_s']


# == plot the distribution of true/control data
plt.figure(figsize=(2.6,2.6))
# scatter with density
x,y = log10_stats, or_at_SE
xy = np.vstack([x,y])
z = gaussian_kde(xy)(xy)
z[np.where(np.isnan(z))] = 0.0
idx = z.argsort()
x, y, z = x[idx], y[idx], z[idx]
g = plt.scatter(x,y,c=z,cmap = plt.cm.GnBu_r,s=3,marker='o')
plt.axhline(y=1,color='k',lw=1.2,ls='--')
plt.axvline(x=0,color='k',lw=1.2,ls='--')
plt.ylabel('Odds ratio at SE')
plt.xlabel('Cluster potential')
plt.ylim([0.25,1.85])
plt.xlim([-330,680])
plt.savefig('{}/SE_OR_vs_cluster_potential.pdf'.format(outdir),bbox_inches='tight',pad_inches=0.02,transparent=True)
plt.show()
plt.close()


# == plot the distribution of true/control data
plt.figure(figsize=(2.6,2.6))
# scatter with density
x,y = log10_stats, or_at_SE
sub_index = y[y>1].index
x,y = x.loc[sub_index], y[sub_index]
# xy = np.vstack([x,y])
# z = gaussian_kde(xy)(xy)
# z[np.where(np.isnan(z))] = 0.0
# idx = z.argsort()
# x, y, z = x[idx], y[idx], z[idx]
# g = plt.scatter(x,y,c=z,cmap = plt.cm.GnBu_r,s=3,marker='o')
g = plt.scatter(x,y,c='k',s=3,marker='o')
# regression
slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)       
x_sort = np.sort(x)
plt.plot(x_sort,x_sort*slope+intercept,c = 'grey',ls='--',lw=.6)
plt.text(.55,1.03,'$R^2$ = {:.2f}'.format(r_value**2),fontsize=13,transform=plt.axes().transAxes)
plt.axhline(y=1,color='k',lw=1.2,ls='--')
plt.axvline(x=0,color='k',lw=1.2,ls='--')
plt.ylim([0.25,1.85])
plt.xlim([-330,680])
plt.ylabel('Odds ratio at SE')
plt.xlabel('Cluster potential')
plt.savefig('{}/SE_OR_GT1_vs_cluster_potential.pdf'.format(outdir),bbox_inches='tight',pad_inches=0.03,transparent=True)
plt.show()
plt.close()


# save the data
pd.concat([df,df_OR],axis=1).to_csv(outdir+os.sep+'data.csv')





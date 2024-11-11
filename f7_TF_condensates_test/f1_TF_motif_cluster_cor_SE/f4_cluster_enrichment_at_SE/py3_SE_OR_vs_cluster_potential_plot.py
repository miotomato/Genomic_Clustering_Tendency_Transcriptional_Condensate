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
matplotlib.rcParams['font.size']=12
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



def mark_pvalue(compr_pos,positions,box_vals):
    s,p = stats.ttest_ind(box_vals[compr_pos[0]],box_vals[compr_pos[1]],nan_policy='omit')
    y, h, col = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),96)*.99 ,1.05, 'k'
    y2 = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),3)*0.99
    x1,x2 = positions[compr_pos[0]],positions[compr_pos[1]]
    # p_label='*'
    p_label='{:.1e}'.format(p)
    if p_label[-2]=='0':
        p_label = p_label[:-2]+p_label[-1]
    if p>=0.05:
        p_label = 'NS'
    if compr_pos[2] == 't':
        plt.plot([x1*1.03, x1*1.03, x2*0.97, x2*0.97], [y, y*h, y*h, y], lw=1, c=col)
        plt.text((x1+x2)*.5, y*h, p_label, ha='center', va='bottom', color=col,fontsize=12)
    else:
        plt.plot([x1*1.03, x1*1.03, x2*0.97, x2*0.97], [y2, y2*.91, y2*.91, y2], lw=1, c=col)
        plt.text((x1+x2)*.5, y2*.95, p_label, ha='center', va='top', color=col,fontsize=12)


# == check the TF motifs
outdir = 'f3_SE_OR_vs_cluster_potential_figs'
os.makedirs(outdir,exist_ok=True)
# ==== read data of cluster potential
infile = '../f2_TF_fimo_jarspar_cluster/f3_closest_distribution_stats/distance_to_closest_site_vs_control_stats.csv'
df = pd.read_csv(infile,index_col=0)
df.index = [i.split('_')[0] for i in df.index]
cutoff_OR2 = df['cutoff of OR>2']
OR_100bp = df['OR of 100bp']
OR_200bp = df['OR of 200bp']
cluster_potential = df['log10 distance t-stats']
raw_stats = df['distance t-stats']
# ==== emrichment/OR at all merged SE
df_OR = pd.read_csv('f1_cluster_enrichment_at_SE/all_hg38_SE.csv',index_col=0)
or_at_SE = df_OR.loc[df.index,'fisher_exact_s']
# ==== save the data
df = pd.concat([df,df_OR],axis=1)
df.to_csv(outdir+os.sep+'data.csv')


# == plot the distribution of true/control data
plt.figure(figsize=(2.6,2.6))
# scatter with density
x,y = cluster_potential, or_at_SE
xy = np.vstack([x,y])
z = gaussian_kde(xy)(xy)
z[np.where(np.isnan(z))] = 0.0
idx = z.argsort()
x, y, z = x[idx], y[idx], z[idx]
g = plt.scatter(x,y,c=z,cmap = plt.cm.GnBu_r,s=3,marker='o')
plt.axhline(y=1,color='k',lw=1.2,ls='--')
plt.axvline(x=0,color='k',lw=1.2,ls='--')
plt.ylabel('Enrichment of \n clustered TFBM at SE')
plt.xlabel('Cluster potential')
# plt.ylim([0.,2])
# plt.xlim([-330,680])
plt.savefig('{}/SE_OR_vs_cluster_potential.pdf'.format(outdir),bbox_inches='tight',pad_inches=0.02,transparent=True)
plt.show()
plt.close()


# == plot the distribution of true/control data
plt.figure(figsize=(2.6,2.6))
# scatter with density
x,y = cluster_potential, or_at_SE
sub_index = y[y>1].index
x,y = x.loc[sub_index], y[sub_index]
g = plt.scatter(x,y,c='k',s=3,marker='o')
# regression
slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)       
x_sort = np.sort(x)
plt.plot(x_sort,x_sort*slope+intercept,c = 'grey',ls='--',lw=.6)
plt.text(.55,1.03,'$R^2$ = {:.2f}'.format(r_value**2),fontsize=13,transform=plt.axes().transAxes)
plt.axhline(y=1,color='k',lw=1.2,ls='--')
plt.axvline(x=0,color='k',lw=1.2,ls='--')
# plt.ylim([0.,2])
# plt.xlim([-330,680])
plt.ylabel('Enrichment of \n clustered TFBM at SE')
plt.xlabel('Cluster potential')
plt.savefig('{}/SE_OR_GT1_vs_cluster_potential.pdf'.format(outdir),bbox_inches='tight',pad_inches=0.03,transparent=True)
plt.show()
plt.close()



# == compr enrichment at SE among groups w/ diff. cluster potential
plt.figure(figsize=(2.6,2.6))
x,y = cluster_potential, or_at_SE
groups = [[-400,0],[0,50],[50,100],[100,150],[150,200],[200,250],[250,300],[300,800]]
xticklabels = ['$<0$','$[0,50)$','[50,100)','[100,150)','[150,200)','[200,250)','[250,300)','$≥ 300$']
# groups = [[-400,0],[0,100],[100,200],[200,300],[300,700]]
# xticklabels = ['$<0$','[0,100)','[100,200)','[200,300)','$≥ 300$']
box_vals = []
positions = np.arange(len(groups))
for group in groups:
    group_min = group[0]
    group_max = group[1]
    sub_index = x[(x>=group_min)&(x<group_max)].index;print(len(sub_index))
    box_vals.append(y.loc[sub_index].values)
# box plot
g = plt.boxplot(box_vals,positions=positions,widths = .6,patch_artist=True,\
            boxprops=dict(color='k',facecolor='w',fill=None,lw=1),\
            medianprops=dict(color='grey'),showfliers=False)    
plt.axhline(y=1,color='k',lw=1.2,ls='--')
plt.ylabel('Enrichment of \n clustered TFBM at SE')
plt.xlabel('Cluster potential')
plt.axes().set_xticks(positions)
plt.axes().set_xticklabels(xticklabels,rotation=60, ha='right',fontsize=13,color='k')
plt.savefig('{}/SE_OR_compr_diff_cluster_potential.pdf'.format(outdir),bbox_inches='tight',pad_inches=0.03,transparent=True)
plt.show()
plt.close()



# == compr cluster potential w/ diff. enrichment at SE 
plt.figure(figsize=(1.8,2.6))
x,y = cluster_potential, or_at_SE
groups = [[-4,1],[1,1.5],[1.5,3],]
xticklabels = ['$<1$','$[1,1.5)$','$≥1.5$']
box_vals = []
positions = np.arange(len(groups))
for group in groups:
    group_min = group[0]
    group_max = group[1]
    sub_index = y[(y>=group_min)&(y<group_max)].index;print(len(sub_index))
    box_vals.append(x.loc[sub_index].values)
# box plot
g = plt.boxplot(box_vals,positions=positions,widths = .6,patch_artist=True,\
            boxprops=dict(color='k',facecolor='w',fill=None,lw=1),\
            medianprops=dict(color='grey'),showfliers=False)    
for compr_pos in [[0,1,'t'],[0,2,'t'],]:
    mark_pvalue(compr_pos,positions,box_vals)
plt.axhline(y=1,color='k',lw=1.2,ls='--')
plt.xlabel('Enrichment of \n clustered TFBS at SE')
plt.ylabel('Cluster potential')
plt.axes().set_xticks(positions)
plt.axes().set_xticklabels(xticklabels,rotation=30, ha='right',fontsize=13,color='k')
plt.savefig('{}/cluster_potential_compr_diff_enrich_at_SSE.pdf'.format(outdir),bbox_inches='tight',pad_inches=0.03,transparent=True)
plt.show()
plt.close()







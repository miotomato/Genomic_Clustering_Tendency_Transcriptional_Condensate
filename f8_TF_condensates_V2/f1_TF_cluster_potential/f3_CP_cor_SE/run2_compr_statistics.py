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
sns.set(font_scale=1.)
sns.set_style("whitegrid", {'axes.grid' : False})
import scipy
import scipy.optimize
sns.set_style("ticks")
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
from matplotlib import gridspec


def mark_pvalue(compr_pos,positions,box_vals):
    s,p = stats.ttest_ind(box_vals[compr_pos[0]],box_vals[compr_pos[1]],nan_policy='omit')
    y, h, col = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),95)*.95 ,1.05, 'k'
    y2 = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),5)*0.99
    x1,x2 = positions[compr_pos[0]],positions[compr_pos[1]]
    indicator = 'down' if s>0 else 'up'
    p_label='{} \n {:.1e}'.format(indicator,p)
    if p_label[-2]=='0':
        p_label = p_label[:-2]+p_label[-1]
    if p>=0.05:
        p_label = 'n.s.'
    if compr_pos[2] == 't':
        plt.plot([x1*1.03, x1*1.03, x2*0.97, x2*0.97], [y, y*h, y*h, y], lw=1, c=col)
        plt.text((x1+x2)*.5, y*h, p_label, ha='center', va='bottom', color=col,fontsize=12)
    else:
        plt.plot([x1*1.03, x1*1.03, x2*0.97, x2*0.97], [y2, y2*.91, y2*.91, y2], lw=1, c=col)
        plt.text((x1+x2)*.5, y2*.95, p_label, ha='center', va='top', color=col,fontsize=12)





outdir = 'f2_compr_statistics_box'
os.makedirs(outdir,exist_ok=True)


indirs = ['CP_TFBS_vs_TFMS','CP_TFBS_with_motif_vs_TFMS','CP_TFBS_without_motif_vs_TFMS']
group_info = {'CP':'log10-distance CP t-stats', 'SE':'fisher_exact_s'}
or_thre = {'CP':9,'SE':6}  
label_name = {'CP':'log2 T-statistics','SE':'log2 Odds Ratio'}

xticklabels = ['All peaks','Peaks w/ motif','Peaks w/o motif']
positions = [0,1,2]
colors = ['grey','grey','grey']

for or_key in ['CP','SE'][:]:
    col = group_info[or_key]
    box_vals = []
    for indir in indirs[:]:
        infile = '../f2_TFMS_TFBS_CP/{}/per_TF_per_Celltype_{}.csv'.format(indir,indir)
        df = pd.read_csv(infile,index_col = 0)
        vals = np.log2(df[col].replace(0,np.nan)).dropna().values
        box_vals.append(vals)
        
        
    ## box with scatter
    plt.figure(figsize=(3,3))
    g = plt.boxplot(box_vals,positions=positions,widths = .6,patch_artist=True,\
                boxprops=dict(color='k',facecolor='w',fill=None,lw=1),\
                medianprops=dict(color='grey'),showfliers=False)    
                
    for patch, color in zip(g['boxes'], colors):
        patch.set_facecolor(color)

    # scatter_X = []
    # for position_id in np.arange(len(positions)):
    #     scatter_x = np.random.normal(positions[position_id],0.07,len(box_vals[position_id]))
    #     plt.scatter(scatter_x,box_vals[position_id],color=colors[position_id],s=20,zorder=0,alpha=0.99,rasterized=True)

    # positions=[0,1,2]
    for compr_pos in [[0,1,'b'], [1,2,'b'],[0,2,'t']]:
        mark_pvalue(compr_pos,positions,box_vals)

    plt.axes().set_xticklabels(xticklabels,rotation=30, ha='right',fontsize=16,color='k')
    plt.ylabel(label_name[or_key],fontsize=16)
    plt.title(or_key)
    plt.savefig('{}/{}.pdf'.format(outdir,or_key),bbox_inches='tight',pad_inches=0.02,transparent=True)
    plt.show()
    plt.close()
    
        






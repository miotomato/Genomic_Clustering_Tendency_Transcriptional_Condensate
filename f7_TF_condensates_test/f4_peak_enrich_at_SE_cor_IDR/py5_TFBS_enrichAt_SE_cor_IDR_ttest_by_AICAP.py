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
matplotlib.rcParams['font.size']=16
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
    y, h, col = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),95)*1.01 ,1.05, 'k'
    y2 = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),2)*0.99
    x1,x2 = positions[compr_pos[0]],positions[compr_pos[1]]
    if p<0.05:
        p_label='{:.1e}'.format(p)
        if p_label[-2]=='0':
            p_label = p_label[:-2]+p_label[-1]
    else:
        p_label='n.s.'
    
    if compr_pos[2] == 't':
        plt.plot([x1*1.03, x1*1.03, x2*0.97, x2*0.97], [y, y*h, y*h, y], lw=1, c=col)
        plt.text((x1+x2)*.5, y*h, p_label, ha='center', va='bottom', color=col,fontsize=16)
    else:
        plt.plot([x1*1.03, x1*1.03, x2*0.97, x2*0.97], [y2, y2*.91, y2*.91, y2], lw=1, c=col)
        plt.text((x1+x2)*.5, y2*.95, p_label, ha='center', va='top', color=col,fontsize=16)





# == check the TF motifs
outdir = 'f5_enrich_at_SE_cor_IDR_ttest_sep_by_AICAP'
os.makedirs(outdir,exist_ok=True)

# read TF AICAP
idrfile = '../f2_TF_motif_IDR_SE_mutation/data/13059_2021_2456_MOESM2_ESM.xlsx'
idrdf = pd.read_excel(idrfile,sheet_name='condition 2(1,6-HD-2)',index_col=0)
# data with SE enrichment
enrich_data_df = pd.read_csv('f1_slurm_log/cistrome_id_log.csv',index_col=0)   
shared_tfs = enrich_data_df.index.intersection(idrdf.index)
enrich_data_df = enrich_data_df.loc[shared_tfs]
enrich_data_df = enrich_data_df.loc[:,enrich_data_df.notnull().sum(axis=0)>2]
enrich_data_df.to_csv(outdir+os.sep+'_cistrome_id_log_AICAP_filtered.csv')

# get shared TF/index
permutation_types = ['sampleMergePeak']
for celltype in enrich_data_df.columns[:]:
    stat_df = pd.read_csv('f2_TFBS_enrichAt_SE_stats/{}_stats.csv'.format(celltype),index_col=0)   
    factor_df = enrich_data_df[celltype].dropna()
    factors = factor_df.index
    for permutation_type in permutation_types:
        rank = stat_df.loc[factors,'{}_relative_rank'.format(permutation_type)]
        aicap = np.log2(idrdf.loc[factors,'AICAP'])
        
        # box plot of AICAP sep by zscore
        plt.figure(figsize=(1.6,2.2))
        low_group = aicap<0
        high_group = aicap>0
        box_vals = [rank[low_group],rank[high_group]]
        # sns.distplot(rank[low_group],kde=False)
        # sns.distplot(rank[high_group],kde=False)
        positions = [1,2]
        colors = ['k','k']
        g = plt.boxplot(box_vals,positions=positions,widths = .6,patch_artist=True,\
                    boxprops=dict(color='k',facecolor='w',fill=None,lw=1),\
                    medianprops=dict(color='grey'),showfliers=False)    
                    
        for patch, color in zip(g['boxes'], colors):
            patch.set_facecolor(color)
        
        scatter_X = []
        for position_id in np.arange(len(positions)):
            scatter_x = np.random.normal(positions[position_id],0.07,len(box_vals[position_id]))
            plt.scatter(scatter_x,box_vals[position_id],color=colors[position_id],s=20,zorder=0,alpha=0.99,rasterized=False)
        for compr_pos in [[0,1,'t']]:
            mark_pvalue(compr_pos,positions,box_vals)

        plt.axes().set_xticklabels(['log2 AICAP<0','log2 AICAP>0',],rotation=45,fontsize=14)
        plt.ylabel('Relative rank')
        plt.title(celltype)
        plt.savefig('{}/{}_{}.pdf'.format(outdir,permutation_type,celltype),bbox_inches='tight',pad_inches=0.02,transparent=True)
        plt.show()
        plt.close()
        
        

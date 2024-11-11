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



def mark_pvalue(compr_pos,positions,box_vals):
    s,p = stats.ttest_ind(box_vals[compr_pos[0]],box_vals[compr_pos[1]],nan_policy='omit')
    y, h, col = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),96)*.99 ,1.05, 'k'
    y2 = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),3)*0.99
    x1,x2 = positions[compr_pos[0]],positions[compr_pos[1]]
    p_label='{:.1e}'.format(p)
    if p_label[-2]=='0':
        p_label = p_label[:-2]+p_label[-1]
    if p>=0.1 or np.isnan(p):
        p_label = 'n.s.'
    if compr_pos[2] == 't':
        plt.plot([x1*1.03, x1*1.03, x2*0.97, x2*0.97], [y, y*h, y*h, y], lw=1, c=col)
        plt.text((x1+x2)*.5, y*h, p_label, ha='center', va='bottom', color=col,fontsize=12)
    else:
        plt.plot([x1*1.03, x1*1.03, x2*0.97, x2*0.97], [y2, y2*.91, y2*.91, y2], lw=1, c=col)
        plt.text((x1+x2)*.5, y2*.95, p_label, ha='center', va='top', color=col,fontsize=12)
    return p



# == check the TF motifs
indir = 'f3_enrhch_at_SE_cor_IDR_mergePeak'
outdir = 'f4_enrich_at_SE_cor_IDR_ttest_sep_by_zscore'
os.makedirs(outdir,exist_ok=True)

idrfile = '../f2_TF_motif_IDR_SE_mutation/data/13059_2021_2456_MOESM2_ESM.xlsx'
sheet_data = {'condition1_1-6HD':['condition 1(1,6-HD-1)','AICAP'],
              'condition2_1-6HD':['condition 2(1,6-HD-2)','AICAP'],
              'condition2_2-5HD':['condition 2(2,5-HD)','AICAP(25HD)']}
permutation_types = ['sampleMergePeak']

for aicap_key in sheet_data.keys():
# for aicap_key in ['condition2_1-6HD']:
    sheet_name = sheet_data[aicap_key][0]
    sheet_col = sheet_data[aicap_key][1]
    idrdf = pd.read_excel(idrfile,sheet_name=sheet_name,index_col=0)
    df_summary = pd.read_csv('{}/_{}_cistrome_id_log_AICAP_filtered.csv'.format(indir,aicap_key),index_col=0)      
    
    # get shared TF/index
    for celltype in df_summary.columns:
        for permutation_type in permutation_types:  
            df = pd.read_csv('{}/_{}_{}_{}.csv'.format(indir,aicap_key,permutation_type,celltype),index_col=0)
            zscore = df['zscore']
            log2_aicap = df['log2 AICAP']
            
            # box plot of AICAP sep by zscore
            plt.figure(figsize=(1.3,2.2))
            # ==== select top 5 vs. others
            group1 = zscore.sort_values(ascending=False)[:5].index
            group2 = zscore.index.difference(group1)
            xticklabels = ['Others','Top5 z-score']
            
            # ==== select top and bottom x%
            # selected_number = int(zscore.shape[0]*.30)
            # group1 = zscore.sort_values(ascending=False)[:selected_number].index
            # group2 = zscore.sort_values(ascending=False)[-1*selected_number:].index
            # xticklabels = ['Low 30% z-score','High 30% z-score']
            
            box_vals = [log2_aicap[group2],log2_aicap[group1]]
            positions = [1,2]
            colors = ['k','k']
            g = plt.boxplot(box_vals,positions=positions,widths = .5,patch_artist=True,\
                        boxprops=dict(color='k',facecolor='w',fill=None,lw=1),\
                        medianprops=dict(color='grey'),showfliers=False)    
            for patch, color in zip(g['boxes'], colors):
                patch.set_facecolor(color)
        
            scatter_X = []
            for position_id in np.arange(len(positions)):
                scatter_x = np.random.normal(positions[position_id],0.07,len(box_vals[position_id]))
                plt.scatter(scatter_x,box_vals[position_id],color=colors[position_id],s=20,zorder=0,alpha=0.99,rasterized=False)
            for compr_pos in [[0,1,'t']]:
                pvalue = mark_pvalue(compr_pos,positions,box_vals)
                # print(pvalue)
    
            plt.axhline(y=0,color='k',lw=1.2,ls='--')
            # plt.axvline(x=0,color='k',lw=1.2,ls='--')
            plt.axes().set_xticklabels(xticklabels,rotation=30,ha='right',fontsize=14)
            # plt.ylabel('log2 AICAP')
            plt.title(celltype)
            plt.savefig('{}/{}_{}_{}.pdf'.format(outdir,aicap_key,permutation_type,celltype),bbox_inches='tight',pad_inches=0.02,transparent=True)
            plt.show()
            plt.close()
            

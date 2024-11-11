import os,sys,argparse,glob
import numpy as np
import pandas as pd
# import find_overlap_keep_info_NOT_sep_strand_asimport
import re
from scipy import stats
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=18
import seaborn as sns
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid' : False,'grid.color': 'grey'})
sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k'})
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams['mathtext.fontset']='custom'
matplotlib.rcParams['mathtext.rm']="Arial"
# import CTCF_TALL_modules_new


def mark_pvalue(compr_pos,positions,box_vals):
    s,p = stats.ttest_ind(box_vals[compr_pos[0]],box_vals[compr_pos[1]],nan_policy='omit')
    y, h, col = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),99.95)*0.95 ,1.1, 'k'
    y2 = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),0)*0.9
    x1,x2 = positions[compr_pos[0]],positions[compr_pos[1]]
    p_label='{:.1e}'.format(p)
    if p_label[-2]=='0':
        p_label = p_label[:-2]+p_label[-1]
    if p<1:
        if p>0.05:
            p_label = 'n.s.'
        if compr_pos[2] == 't':
            plt.plot([x1*1.03, x1*1.03, x2*0.97, x2*0.97], [y, y*h, y*h, y], lw=1, c=col)
            plt.text((x1+x2)*.5, y*1.2, p_label, ha='center', va='center', color=col,fontsize=14)
        else:
            plt.plot([x1*1.03, x1*1.03, x2*0.97, x2*0.97], [y2, y2*1.1, y2*1.1, y2], lw=1, c=col)
            plt.text((x1+x2)*.5, y2*1.25, p_label, ha='center', va='top', color=col,fontsize=14)
    return s,p



# == main 
indir = 'f3_diff_ATAC_overlap_SE'
outdir = 'f4_diff_atac_figs'
os.makedirs(outdir,exist_ok=True)

# read matched names between TCGA ATAC and SE data
name_match = pd.read_excel('TCGA_ATAC_SE_cancerType_match.xlsx',index_col=0)   
name_match = name_match.dropna()
pyfile="find_overlap_keep_info_NOT_sep_strand_revised.py"

for cancertype in name_match.index[:]:
    cancertype_SE = name_match.loc[cancertype].SE
    cancertype_SE_rename = name_match.loc[cancertype].SE_rename
    overlap_file = '{}/{}_ATAC_overlap_{}_SE.bed'.format(indir,cancertype,cancertype_SE)
    nonoverlap_file = '{}/{}_ATAC_NOT_overlap_{}_SE.bed'.format(indir,cancertype,cancertype_SE)
    overlap_df = pd.read_csv(overlap_file,sep='\t')
    nonoverlap_df = pd.read_csv(nonoverlap_file,sep='\t')
    
    # read the diff atac-seq sig on and off SE
    box_vals = [overlap_df.stats,nonoverlap_df.stats]
    positions = [1,2]
    colors = [ 'red','darkgrey']
    
    plt.figure(figsize=(2.,2.7))
    g = plt.boxplot(box_vals,positions=positions,widths = .55,patch_artist=True,\
                boxprops=dict(color='k',facecolor='w',fill=None,lw=1),\
                medianprops=dict(color='k'),showfliers=False)
    
    scatter_X = []
    for position_id in np.arange(len(positions)):
        scatter_x = np.random.normal(positions[position_id],0.06,len(box_vals[position_id]))
        plt.scatter(scatter_x,box_vals[position_id],color=colors[position_id],
                    s=22,zorder=0,alpha=1,marker='o',rasterized=True)
    
    for compr_pos in [[0,1,'t']]:
        s,p = mark_pvalue(compr_pos,positions,box_vals)
        # print(cancertype,flag,compr_pos,s,p)
    
    plt.title(cancertype)
    plt.ylabel('$\Delta$ ATAC-seq score \n {} over others'.format(cancertype))
    plt.axes().set_xticklabels(['{:,} ATAC-seq peaks \n overlap w/ {} SE'.format(overlap_df.shape[0],cancertype_SE_rename),
                                '{:,} ATAC-seq peaks \n NOT overlap w/ {} SE'.format(nonoverlap_df.shape[0],cancertype_SE_rename)],
                               rotation=45,ha='center',fontsize=14)
    plt.axhline(y=0,color='grey',linestyle='--',linewidth=.7)
    plt.savefig('{}/{}_diff_ATAC_compr.pdf'.format(outdir,cancertype),
                bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
    plt.show()
    plt.close()



 
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
sns.set(font_scale=1.1)
sns.set_style("whitegrid", {'axes.grid' : False})
import scipy
import scipy.optimize
sns.set_style("ticks")
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams["mathtext.rm"] = "Arial"
from scipy.interpolate import interpn
from scipy.stats import gaussian_kde



def scatter_plot(x,y,xlabel,ylabel,figname,title):

    # ==== scatter plot with density
    plt.figure(figsize=(1.8,2.5))
    plt.scatter(x,y,c='k',s=9)
    # label the TFs
    label_i=0
    label_indexes = y.sort_values(ascending=False)[:half_len].index
    label_indexes = label_indexes.append(x.sort_values(ascending=False)[:half_len].index)
    label_indexes = y.loc[set(label_indexes)].sort_values(ascending=False).index
    for label_index in label_indexes:
        plt.scatter(x[label_index],y[label_index],c=colors[label_i],
                    s=25,label=label_index)
        label_i+=1
    plt.legend(bbox_to_anchor=[.0,-.3],
           markerscale=1.2,fontsize=11,borderaxespad=0.2,labelspacing=.2,
           handletextpad=0.2,handlelength=1.,loc="upper left",frameon=False)
    # ==== linear regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)       
    x_sort = np.sort(x)
    plt.plot(x_sort,x_sort*slope+intercept,c = 'grey',ls='--',lw=.6)
    plt.text(.04,.88,'$R^2$ = {:.2f}'.format(r_value**2),fontsize=11,transform=plt.axes().transAxes)
    plt.axhline(y=0,color='grey',lw=.5,ls='--')
    plt.axvline(x=0,color='grey',lw=.5,ls='--')
    plt.ylim([0,135])
    plt.xlim([-230,230])
    plt.axes().set_xticks([-200,0,200])
    plt.axes().set_xticklabels(['-200'.rjust(9),0,'200'.ljust(6)])

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    if not title=='GM12878':
        plt.yticks([])
        plt.ylabel('')

    plt.savefig('{}'.format(figname),bbox_inches='tight',pad_inches=0.02,transparent=True)
    plt.show()
    plt.close()
    




# == check the TF motifs
outdir = 'f12_sameCT_TF_top10k_CP_cor_allPeak_CP'
os.makedirs(outdir,exist_ok=True)

# ==== TFMS CP
# tfms = pd.read_csv('data/TFMS_CP_SE_enrich.csv',index_col=0)
# tfms.index = [i.split('_')[0] for i in tfms.index]

# ==== TFBS info
tfbp_cp_df_top10k = pd.read_csv('../f1_TFMS_TFBS_CP/CP_TFBS_top10k_vs_TFMS/per_data_CP_TFBS_top10k_vs_TFMS.csv',index_col=0)
tfbp_cp_df_all = pd.read_csv('../../../f9_TF_condensates_V3/f1_TF_cluster_potential/f1_TFMS_TFBS_CP/CP_TFBS_vs_TFMS/per_data_CP_TFBS_vs_TFMS.csv',index_col=0)
# tfbs_dir = 'f3_TFBS_CP_heatmap/_csv/'
tfbs_cp_types = ['CP_TFBS_vs_TFMS',
                  'CP_TFBS_overlap_motif_vs_TFMS',
                  'CP_TFBS_NOT_overlap_motif_vs_TFMS']


for tfbs_cp_type in tfbs_cp_types[:1]:
    # basename = '{}_{}_{}'.format(or_key,fisher_P_key,indir)
    # basename = 'data_fisher_{}_CP_RankSum'.format(tfbs_cp_type)
    # tfbs_s = pd.read_csv('{}/{}_statistics.csv'.format(tfbs_dir,basename),index_col=0)
    # tfbs_p = pd.read_csv('{}/{}_P.csv'.format(tfbs_dir,basename),index_col=0)
    df = pd.read_csv('../../data/cistrome/cistrome2019_data_Count_with_SE_motif.csv',sep=',',index_col=0)
    df = df.replace(0,np.nan)
    # ==== TFBS CP for each CT
    for ct in df.columns[:]:
        y = df[ct].dropna()
        if len(y)<=3:
            continue
        outdf = pd.DataFrame()
        plt.figure(figsize=(2.6,2.6))
        color_i=0
        for factor in y.index[:]:
            kept_data = [i for i in tfbp_cp_df_top10k.index if re.search('{}_{}_'.format(ct,factor),i)]
            if len(kept_data)<=3:
                continue
            # print(ct,factor,len(kept_data))
            # px = tfbp_cp_df.loc[kept_data,'#TFBS_top10k']
            # px = np.log10(px)+color_i
            px = tfbp_cp_df_top10k.loc[kept_data,'log10-dis Wilcoxon-rank-sum-s']
            all_peak_df_index = ['_'.join(i.split('_')[::-1]) for i in kept_data]
            py = tfbp_cp_df_all.loc[all_peak_df_index,'log10-dis Wilcoxon-rank-sum-s']
            
            plt.scatter(px,py,color=plt.cm.tab20(color_i),s=9,label=factor)
            slope, intercept, r_value, p_value, std_err = stats.linregress(px, py)      
            outdf.loc[factor,'r-value'] = r_value
            color_i+=1
        
        # print(ct,color_i)
        plt.legend(bbox_to_anchor=[1,1],
           markerscale=2,fontsize=11,borderaxespad=0.2,labelspacing=.2,
           handletextpad=0.2,handlelength=1.,loc="upper left",frameon=False)
        plt.xlabel('Top10k TFBS CP')
        plt.ylabel('All peaks TFBS CP')
        plt.title(ct)
        if color_i>=2:
            plt.savefig('{}/{}.pdf'.format(outdir,ct),bbox_inches='tight',pad_inches=0.02,transparent=True)
            outdf.to_csv('{}/_{}.csv'.format(outdir,ct))
        plt.show()
        plt.close()
            
            
    
    
    
    
    
    
    
    
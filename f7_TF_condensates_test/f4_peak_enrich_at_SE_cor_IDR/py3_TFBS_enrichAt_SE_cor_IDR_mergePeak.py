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
sns.set(font_scale=1.)
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
outdir = 'f3_enrhch_at_SE_cor_IDR_mergePeak'
os.makedirs(outdir,exist_ok=True)

# read TF AICAP
idrfile = '../f2_TF_motif_IDR_SE_mutation/data/13059_2021_2456_MOESM2_ESM.xlsx'
sheet_data = {'condition1_1-6HD':['condition 1(1,6-HD-1)','AICAP'],
              'condition2_1-6HD':['condition 2(1,6-HD-2)','AICAP'],
              'condition2_2-5HD':['condition 2(2,5-HD)','AICAP(25HD)']}
# permutation_types = ['shuffleGenome','sampleUDHS','sampleMergePeak']
permutation_types = ['sampleMergePeak']

for aicap_key in sheet_data.keys():
    sheet_name = sheet_data[aicap_key][0]
    sheet_col = sheet_data[aicap_key][1]
    idrdf = pd.read_excel(idrfile,sheet_name=sheet_name,index_col=0)
    # idrdf = idrdf[idrdf['AICAP t-test Pvalue']<0.1]
    # data with SE enrichment
    enrich_data_df = pd.read_csv('f1_slurm_log/cistrome_id_log.csv',index_col=0)   
    shared_tfs = enrich_data_df.index.intersection(idrdf.index)
    enrich_data_df = enrich_data_df.loc[shared_tfs]
    enrich_data_df = enrich_data_df.loc[:,enrich_data_df.notnull().sum(axis=0)>5]
    enrich_data_df.to_csv(outdir+os.sep+'_{}_cistrome_id_log_AICAP_filtered.csv'.format(aicap_key))
    # lable_total = []
    
    # get shared TF/index
    colors = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple',
              'tab:brown','tab:pink','tab:grey','tab:olive','tab:cyan']
    
    for celltype in enrich_data_df.columns[:]:
        stat_df = pd.read_csv('f2_TFBS_enrichAt_SE_stats/{}_stats.csv'.format(celltype),index_col=0)   
        factor_df = enrich_data_df[celltype].dropna()
        factors = factor_df.index
        for permutation_type in permutation_types:
            zscore = stat_df.loc[factors,'{}_zscore'.format(permutation_type)]
            aicap = idrdf.loc[factors,sheet_col]
            log2_aicap = np.log2(aicap)
            zscore_rank = zscore.rank(ascending=False)
            aicap_rank = log2_aicap.rank(ascending=True)
            avg_rank = pd.concat([zscore_rank,aicap_rank],axis=1).mean(axis=1)
            # ==== write out the data
            df_out = pd.DataFrame()
            df_out['zscore']= zscore
            df_out['log2 AICAP']= log2_aicap
            df_out['avg rank']= avg_rank
            df_out = df_out.sort_values(by='avg rank',ascending=True)
            df_out.to_csv('{}/_{}_{}_{}.csv'.format(outdir,aicap_key,permutation_type,celltype))
            
            # == plot 
            plt.figure(figsize=(2.,2.))
            x = zscore
            y = log2_aicap
            plt.scatter(x,y,c='k',s=11)
            plt.axhline(y=0,color='k',lw=1.2,ls='--')
            plt.axvline(x=0,color='k',lw=1.2,ls='--')
            plt.xlabel('Enrichment (z-score) at SE',fontsize=13)
            plt.ylabel('log2 AICAP',fontsize=13)
            plt.title(celltype)
            # highlight top few factors
            # label_indexes1 = x.sort_values(ascending=False)[:5].index
            # label_indexes2 = y.sort_values()[:5].index
            # label_indexes = label_indexes1.append(label_indexes2).unique()
            # label_indexes = y[label_indexes].sort_values(ascending=False).index
            label_i=0
            label_indexes = avg_rank.sort_values()[:10].index
            for label_index in label_indexes:
                # plt.text(x[label_index],y[label_index],label_index,fontsize=7)
                plt.scatter(x[label_index],y[label_index],c=colors[label_i],
                            s=25,label=label_index)
                label_i+=1
            plt.legend(bbox_to_anchor=[.98,1.05],
                   markerscale=1.2,fontsize=10,borderaxespad=0.2,labelspacing=.2,
                   handletextpad=0.2,handlelength=1.,loc="upper left",frameon=False)
            plt.savefig('{}/{}_{}_{}.pdf'.format(outdir,aicap_key,permutation_type,celltype),bbox_inches='tight',pad_inches=0.02,transparent=True)
            plt.show()
            plt.close()
            
            # lable_total.append(label_indexes.values)
            # len(set([i for j in lable_total for i in j])) = 37
            
            
            

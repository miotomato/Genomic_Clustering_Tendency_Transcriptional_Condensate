import sys,argparse
import os,glob
import numpy as np
import pandas as pd
from scipy import stats
import re,bisect

import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=16
import seaborn as sns
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k'})
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams["mathtext.rm"] = "Arial"
from scipy import stats
import warnings
warnings.filterwarnings("ignore")



chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',\
             'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',\
             'chr18','chr19','chr20','chr21','chr22','chrX','chrY']

    

indir = 'f6_avg_RP_compr_mutation'
outdir = 'f6b_avg_RP_compr_mutation_figs'
os.makedirs(outdir,exist_ok=True)
 

# genome_control = {'Rate':13474.4*6285/3298912062,
#                   'Count':76088647/3298912062}
                   
name_match = pd.read_excel('../../../f9_TF_condensates_V3/data/TCGA/TCGA-ATAC_SE_cancerType_match.xlsx',index_col=0)   
name_match = name_match.dropna()
cts = ['MCF-7','HCT-116','HeLa','LNCaP','U87','HepG2']
# datadir = '../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged/'
# tfbs_dir = '../f1_TF_cluster_potential/f2_cor_CP_SE_AICAP/f3_TFBS_CP_heatmap//_csv/'
# tfbs_cp_s = pd.read_csv('{}/data_fisher_CP_TFBS_vs_TFMS_CP_RankSum_statistics.csv'.format(tfbs_dir),index_col=0)


rank_dir = '../../f1_TF_cluster_potential/f2_cor_CP_SE_AICAP/f9_per_CT_TFBS_CP_cor_zscore_CP_with_motif_SE/TFBS_CP/'

halflifes = [10000,5000,20000]

flags = ['percentile_T.merge',
          'percentile_C.merge',
          'percentile_T_ExtendMerge.merge',
          'percentile_T.merge.SE_overlapped',
          'percentile_C.merge.SE_overlapped',
          'percentile_T_ExtendMerge.merge.SE_overlapped',]

flags = ['percentile_T.merge',
         'percentile_T_ExtendMerge.merge',]


thre = None
for halflife in halflifes[:1]:
    for flag in flags[:]:
        for cancertype in ['BRCA', 'CESC', 'COAD','LIHC', 'PRAD'][:]:
            # os.makedirs(outdir+os.sep+cancertype,exist_ok=True)
            ct = name_match.loc[cancertype,'SE']
            rank_df = pd.read_csv('{}/_CP_TFBS_all_vs_TFMS_{}.csv'.format(rank_dir,ct),index_col=0)
            factors = rank_df.index
            patient_df = pd.read_csv('f5_avg_RP_across_regions/{}/{}_avg_RP_halflife_{}_on_{}_{}_{}.csv'.format(cancertype,cancertype,halflife,ct,factors[0],flag),index_col=0,header=None)    
                
            df = pd.read_csv('{}/{}_mutation_by_RP_halflife_{}_{}.csv'.format(indir,cancertype,halflife,flag),index_col=0)

            data_col = ['% high RP w/ mutations','% low RP w/ mutations',]            
            data_col = ['# high RP w/ mutations','# low RP w/ mutations',]            
        
            positions = np.arange(df.shape[0])
            plt.figure(figsize=(df.shape[0]/3.5,2.6))
            a = plt.bar(positions-.2,df[data_col[0]],width=.3,color='salmon')
            b = plt.bar(positions+.2,df[data_col[1]],width=.3,color='royalblue')
            # plt.bar(positions,all_DNAme_none,bottom=all_DNAme_de,width=.9,linewidth=0,color='lightgrey')
            plt.xlim([-1,df.shape[0]])
            # plt.ylim([0.5,1])
            # plt.axhline(y=genome_control[mutationType],color='k',lw=1.2,ls='--')
            plt.ylabel('# patients w/ muration \n (Total = {})'.format(patient_df.shape[0]))
            plt.legend([a,b],['Top 50% patients w/ high RP','Bottom 50% patients w/ low RP'],fontsize=10,bbox_to_anchor=[1,1.2],
                       borderaxespad=0.1,labelspacing=.1,loc="upper right",frameon=False)
            plt.axes().set_xticks(positions)
            plt.axes().set_xticklabels(df.index,rotation=60, ha='right',fontsize=12,color='k')
            if cancertype in [ 'CESC', 'LIHC', 'PRAD']:
                plt.axes().set_yticks([0,1])
            plt.title(cancertype,loc='left')
            figname = '{}/{}_mutation_by_RP_halflife_{}_{}.pdf'.format(outdir,cancertype,halflife,flag )
            plt.savefig(figname,bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
            plt.show()
            plt.close()
            
            s,p = stats.ttest_ind(df[data_col[0]],df[data_col[1]])
            print(s,p)
    
    
    




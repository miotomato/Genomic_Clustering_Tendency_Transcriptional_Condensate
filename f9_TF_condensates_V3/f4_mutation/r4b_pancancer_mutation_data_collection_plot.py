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

    

indir = 'f4_pancancer_mutation_data_collection'
outdir = 'f4b_pancancer_mutation_data_collection_fig'
os.makedirs(outdir,exist_ok=True)
 

genome_control = {'Rate':100*13474.4*6285/3298912062,
                  'Count':100*76088647/3298912062}
                   
cts = ['MCF-7','HCT-116','HeLa','LNCaP','U87','HepG2']
datadir = '../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged/'
tfbs_dir = '../f1_TF_cluster_potential/f2_cor_CP_SE_AICAP/f3_TFBS_CP_heatmap//_csv/'
# tfbs_cp_s = pd.read_csv('{}/data_fisher_CP_TFBS_vs_TFMS_CP_RankSum_statistics.csv'.format(tfbs_dir),index_col=0)

treat_flags = ['percentile_T','percentile_T_ExtendMerge']

for mutationType in ['Rate','Count']:
    for treat_flag in treat_flags[:]:
        for ct in cts[:]:
            # os.makedirs(outdir+os.sep+ct,exist_ok=True)
            # factors = tfbs_cp_s[ct].dropna().index
            df = pd.read_csv(indir+os.sep+'{}_murationRate.csv'.format(ct),index_col=0)
            df = df.drop(df.index.intersection(['T']))
            data_col = ['{} %mutation{} per bp'.format(treat_flag,mutationType),
                        'percentile_C %mutation{} per bp'.format(mutationType)]            
        
            positions = np.arange(df.shape[0])
            plt.figure(figsize=(df.shape[0]/3.5,2.6))
            a = plt.bar(positions-.2,100*df[data_col[0]],width=.3,color='salmon')
            b = plt.bar(positions+.2,100*df[data_col[1]],width=.3,color='royalblue')
            # plt.bar(positions,all_DNAme_none,bottom=all_DNAme_de,width=.9,linewidth=0,color='lightgrey')
            plt.xlim([-1,df.shape[0]])
            # plt.ylim([0.5,1])
            plt.axhline(y=genome_control[mutationType],color='k',lw=1.2,ls='--')
            plt.ylabel('Muration Rate')
            plt.legend([a,b],['C-TFBS','NC-TFBS'],fontsize=10,bbox_to_anchor=[1,1.2],
                       borderaxespad=0.1,labelspacing=.1,loc="upper right",frameon=False)
            plt.axes().set_xticks(positions)
            plt.axes().set_xticklabels(df.index,rotation=60, ha='right',fontsize=12,color='k')
            plt.title(ct,loc='left')
            figname = '{}/{}_{}_mutation{}.pdf'.format(outdir,ct,treat_flag,mutationType)
            plt.savefig(figname,bbox_inches='tight',pad_inches=0.1,dpi=600)
            plt.show()
            plt.close()
            
            s,p = stats.ttest_ind(100*df[data_col[0]],100*df[data_col[1]])
            print(s,p)
    
    
    




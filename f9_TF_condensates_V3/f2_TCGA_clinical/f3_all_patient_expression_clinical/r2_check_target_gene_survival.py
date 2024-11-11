import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import json
#from GenomeData import *
import re,bisect
from lifelines.statistics import logrank_test
from lifelines import KaplanMeierFitter

import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=12
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=1)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style('ticks')

from scipy import stats


   
# ==== main  
outdir = 'f2_target_gene_clinical'
os.makedirs(outdir,exist_ok=True)
# os.makedirs(outdir+os.sep+'fig',exist_ok=True)
indir = '../f2_ATAC_RP_clinical/f1_gene_fpkm-uq_clinical_by_TFRP_topCP/_csv/'

# factors_dic = {'BRCA':['ERG','E2F1'],
#                 'COAD':['SRF','JUND']}

factors_dic = {'BRCA':['NR2F2','ERG','MAX','CEBPB','ELF1'],
                'COAD':['JUND','CEBPB','MAX','FOSL1','SRF']}


# factors_dic = {'BRCA':['CoTarget_2'],
#                 'COAD':['CoTarget_2']}


for cancertype in ['BRCA','COAD'][:]:   
    clinical_df = pd.read_csv('f1_clinical_per_gene/{}_logrank_info.csv'.format(cancertype),index_col=0)
    factors = factors_dic[cancertype]
    target_genes = []
    
    out_df = pd.DataFrame()
    for factor in factors[:2]:       
        infile = '{}/{}_{}_targets.csv'.format(indir,cancertype,factor)       
        infile_df = pd.read_csv(infile,index_col=0)
        infile_df.columns = ['{} {}'.format(factor,i) for i in infile_df.columns]
        target_genes = np.append(target_genes,infile_df.index)
        out_df = pd.concat([out_df,infile_df],axis=1)
        
    # df = clinical_df.loc[out_df.index].drop_duplicates()
    out_df = pd.concat([out_df,clinical_df.loc[out_df.index]],axis=1)
    out_df = out_df.sort_values(by='log rank p')
    out_df.to_csv('{}/{}_clinical.csv'.format(outdir,cancertype))    
    
    a,b = clinical_df.shape[0],clinical_df[clinical_df['log rank p']<0.05].shape[0]
    print(cancertype,'genome', np.round(b/a,3))
    c,d = out_df.shape[0],out_df[out_df['log rank p']<0.05].shape[0]
    print(cancertype,'target genes', np.round(d/c,3))
    
    s,p = stats.fisher_exact([[d,c-d],[b,a-b]])
    print('odds ratio=',s.round(3),'p=',p.round(3),'\n')
    
    
    
    # expr_df = pd.read_csv('../../data/TCGA//GDC-PANCAN.htseq_fpkm-uq.{}.tsv'.format(cancertype),sep='\t',index_col=0)



# x=expr_df.iloc[:,11]
# y=expr_df.iloc[:,111]
# plt.scatter(x,y)







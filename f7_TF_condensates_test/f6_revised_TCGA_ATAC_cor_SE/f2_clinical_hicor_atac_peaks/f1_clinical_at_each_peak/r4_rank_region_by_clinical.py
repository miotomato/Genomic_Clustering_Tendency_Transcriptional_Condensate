import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import json
#from GenomeData import *
import re,bisect
from lifelines.statistics import logrank_test
from lifelines import KaplanMeierFitter

# import matplotlib
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt
# matplotlib.rcParams['font.size']=16
# matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
# matplotlib.rcParams["font.family"] = "sans-serif"
# import seaborn as sns
# sns.set(font_scale=1.2)
# sns.set_style("whitegrid", {'axes.grid' : False})
# sns.set_style('ticks')




    
# ==== main  
indir = 'f2_caseID_each_peak_vs_clinical'
outdir = 'f4_rank_region_by_clinical_pvalue'
os.makedirs(outdir,exist_ok=True)
name_match = pd.read_excel('../../data/TCGA/TCGA-ATAC_SE_cancerType_match.xlsx',index_col=0)   
name_match = name_match.dropna()
filtered_df = pd.read_excel('../../data/TCGA/TCGA-ATAC_clustered_samples.xlsx',index_col=0)   
sig_file = '../../data/TCGA/mynorm_TCGA-ATAC_PanCan_Log2_QuantileNorm_Counts_plus5.caseID.avg.txt'


for cancertype in ['BRCA','COAD'][:]:
    cancertype_SE = name_match.loc[cancertype].SE
    cancertype_SE_rename = name_match.loc[cancertype].SE_rename
    filtered_id = filtered_df[filtered_df.cohort==cancertype].case_id
    # ==== read TCGA peaks overlap with SE data
    with open(sig_file) as sig_inf:
        sig_df = pd.read_csv(sig_inf,sep='\t',index_col=3)
    sig_df = sig_df.iloc[:,:3]
    
    # == read log rank info
    clinical_df = pd.read_csv('{}/{}_logrank_info.csv'.format(indir,cancertype),index_col=0)
    clinical_df['neg_logP'] = -1*np.log10(clinical_df['log rank p'])
    
    sig_df = sig_df.loc[clinical_df.index]
    sig_df.insert(3,'name',sig_df.index)
    sig_df.insert(4,'score',clinical_df['neg_logP'])
    sig_df = sig_df.sort_values(by='score',ascending=False)
    sig_df.to_csv('{}/{}_clinical_ranked_regions.bed'.format(outdir,cancertype),
                                          sep='\t',header=None,index=False)
        

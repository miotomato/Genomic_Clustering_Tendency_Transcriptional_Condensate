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
indir1 = '../../f1_clinical_at_each_SE/f1_ATAC_overlap_SE_caseID'
indir2 = '../../f1_clinical_at_each_SE/f2_caseID_each_SE_vs_clinical'
outdir = './'
os.makedirs(outdir,exist_ok=True)
name_match = pd.read_excel('../../../data/TCGA/TCGA-ATAC_SE_cancerType_match.xlsx',index_col=0)   
name_match = name_match.dropna()
filtered_df = pd.read_excel('../../../data/TCGA/TCGA-ATAC_clustered_samples.xlsx',index_col=0)   


for cancertype in ['BRCA','COAD'][:]:
    cancertype_SE = name_match.loc[cancertype].SE
    cancertype_SE_rename = name_match.loc[cancertype].SE_rename
    filtered_id = filtered_df[filtered_df.cohort==cancertype].case_id
    # ==== read TCGA peaks overlap with SE data
    overlap_file = '{}/{}_ATAC_overlap_{}_SE_caseID.bed'.format(indir1,cancertype,cancertype_SE_rename)
    with open(overlap_file) as sig_inf:
        sig_df = pd.read_csv(sig_inf,sep='\t',index_col=3)
    sig_df = sig_df.iloc[:,:6]
    
    # == read log rank info
    clinical_df = pd.read_csv('{}/{}_logrank_info.csv'.format(indir2,cancertype),index_col=0)
    clinical_df = clinical_df[clinical_df['log rank p']<0.05]
    clinical_df = clinical_df[clinical_df['treat time']<clinical_df['ctrl time']]
    
    sig_df.insert(3,'name',sig_df.index)
    sig_df.loc[clinical_df.index].to_csv('{}/{}_clinical_hicor_regions.bed'.format(outdir,cancertype),
                                         sep='\t',header=None,index=False)
        

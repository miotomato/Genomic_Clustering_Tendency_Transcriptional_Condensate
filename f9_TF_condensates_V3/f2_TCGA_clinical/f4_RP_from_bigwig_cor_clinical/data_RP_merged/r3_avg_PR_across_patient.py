import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *
from scipy import stats
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#matplotlib.rcParams['font.size']=16
#matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
#matplotlib.rcParams["font.family"] = "sans-serif"
#import seaborn as sns
#sns.set(font_scale=2)
#sns.set_style("whitegrid", {'axes.grid' : False})
import re,bisect
#plus = re.compile('\+')
#minus = re.compile('\-')


# ==== main
indir='f2_avg_RP_per_sample'
outdir = 'f3_avg_RP_per_sample_across_patients'
os.makedirs(outdir,exist_ok=True)

# atac_id_info = pd.read_csv('../../data/TCGA/TCGA_identifier_mapping.txt',sep='\t',index_col=3)
filtered_df = pd.read_excel('../../../data/TCGA/TCGA-ATAC_clustered_samples.xlsx',index_col=0)   
# data_types = ['tcga_atac','HCT-116.merge','MCF-7.merge']

name_match = pd.read_excel('../../../data/TCGA/TCGA-ATAC_SE_cancerType_match.xlsx',index_col=0)   
name_match = name_match.dropna()

tfbs_dir = '../../..//f1_TF_cluster_potential/f2_cor_CP_SE_AICAP/f3_TFBS_CP_heatmap//_csv/'
tfbs_cp_s = pd.read_csv('{}/data_fisher_CP_TFBS_vs_TFMS_CP_RankSum_statistics.csv'.format(tfbs_dir),index_col=0)

halflifes = [10000,5000,20000]

flags = ['percentile_T.merge',
         'percentile_C.merge',
         'percentile_T_ExtendMerge.merge',
         'percentile_T.merge.SE_overlapped',
         'percentile_C.merge.SE_overlapped',
         'percentile_T_ExtendMerge.merge.SE_overlapped',]

# for cancertype in name_match.index[:]:
for halflife in halflifes[:]:
    for flag in flags[:]:
        for cancertype in ['BRCA', 'CESC', 'COAD','LIHC', 'PRAD'][:]:
            os.makedirs(outdir+os.sep+cancertype,exist_ok=True)
            ct = name_match.loc[cancertype,'SE']
            # filtered_ids = filtered_df[filtered_df.cohort==cancertype].index
            factors = tfbs_cp_s[ct].dropna().index

            for factor in factors[:]:
                print(halflife, flag, cancertype, factor)
                df = pd.read_csv('{}/{}/{}_avg_RP_halflife_{}_on_{}_{}_{}.csv'.format(indir,cancertype,cancertype,halflife,ct,factor,flag),index_col=0)    
                df.mean(axis=1).to_csv('{}/{}/{}_avg_RP_halflife_{}_on_{}_{}_{}.csv'.format(outdir,cancertype,cancertype,halflife,ct,factor,flag))    








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
indir='f1_get_RP_from_bedGraph'
outdir='f2_avg_RP_per_sample'
os.makedirs(outdir,exist_ok=True)

# atac_id_info = pd.read_csv('../../data/TCGA/TCGA_identifier_mapping.txt',sep='\t',index_col=3)
filtered_df = pd.read_excel('../../../f9_TF_condensates_V3/data/TCGA/TCGA-ATAC_clustered_samples.xlsx',index_col=0)   
# data_types = ['tcga_atac','HCT-116.merge','MCF-7.merge']

name_match = pd.read_excel('../../../f9_TF_condensates_V3/data/TCGA/TCGA-ATAC_SE_cancerType_match.xlsx',index_col=0)   
name_match = name_match.dropna()

tfbs_dir = '../..//f1_TF_cluster_potential/f2_cor_CP_SE_AICAP/f3_TFBS_CP_heatmap//_csv/'
tfbs_cp_s = pd.read_csv('{}/data_fisher_CP_TFBS_all_vs_TFMS_CP_KStest_statistics.csv'.format(tfbs_dir),index_col=0)

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
#         for cancertype in ['BRCA', 'CESC', 'COAD','GBMx', 'LIHC', 'PRAD'][:]:
        for cancertype in ['BRCA', 'CESC', 'COAD', 'LIHC', 'PRAD'][:]:
            # os.makedirs(outdir+os.sep+cancertype,exist_ok=True)
            os.makedirs(outdir+os.sep+cancertype,exist_ok=True)
            ct = name_match.loc[cancertype,'SE']
            filtered_ids = filtered_df[filtered_df.cohort==cancertype].index
            factors = tfbs_cp_s[ct].dropna().index
        
            for factor in factors[:]:
                # ==== averaged RP from all replicates per cancer per data type
                df_cancertype = pd.DataFrame()
                for filtered_id in filtered_ids[:]:
                    # ==== average RP from all replicates per patient
                    prename = '{}_{}_*_halflife_{}_on_{}_{}_{}.csv'.format(cancertype,filtered_id,halflife,ct,factor,flag)
                    infiles = glob.glob('{}/{}/{}'.format(indir,cancertype,prename))
                    if len(infiles)==0:
                        continue
                    print(cancertype,factor,flag,halflife,filtered_id,len(infiles)) 
                    print('\n'.join(infiles),'\n')
    
                    df_id = pd.DataFrame()
                    for infile in infiles:
                        df = pd.read_csv(infile,sep='\t',index_col = 0)
                        df_id = pd.concat([df_id,df],axis=1)
                    df_id = df_id.mean(axis=1).rename(filtered_id)
                    # ==== average RP per patient 
                    df_cancertype = pd.concat([df_cancertype,df_id],axis=1)
                
                df_cancertype.to_csv('{}/{}/{}_avg_RP_halflife_{}_on_{}_{}_{}.csv'.format(outdir,cancertype,cancertype,halflife,ct,factor,flag))    
        
        

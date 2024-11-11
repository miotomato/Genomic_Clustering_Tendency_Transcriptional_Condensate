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
indir = 'f5_avg_RP_across_regions'
outdir = 'f6_avg_RP_compr_mutation'
os.makedirs(outdir,exist_ok=True)

# atac_id_info = pd.read_csv('../../data/TCGA/TCGA_identifier_mapping.txt',sep='\t',index_col=3)
filtered_df = pd.read_excel('../../../f9_TF_condensates_V3/data/TCGA/TCGA-ATAC_clustered_samples.xlsx',index_col=0)   
# data_types = ['tcga_atac','HCT-116.merge','MCF-7.merge']

name_match = pd.read_excel('../../../f9_TF_condensates_V3/data/TCGA/TCGA-ATAC_SE_cancerType_match.xlsx',index_col=0)   
name_match = name_match.dropna()

tfbs_dir = '../..//f1_TF_cluster_potential/f2_cor_CP_SE_AICAP/f3_TFBS_CP_heatmap_with_motif_SE//_csv/'
tfbs_cp_s = pd.read_csv('{}/data_fisher_CP_TFBS_all_vs_TFMS_CP_KStest_statistics.csv'.format(tfbs_dir),index_col=0)

rank_dir = '../../f1_TF_cluster_potential/f2_cor_CP_SE_AICAP/f9_per_CT_TFBS_CP_cor_zscore_CP_with_motif_SE/TFBS_CP/'
mutation_dir = '../../../f9_TF_condensates_V3/data/TCGA/mutation/icgc_mutation_bed_hg38'


halflifes = [10000,5000,20000]

flags = ['percentile_T.merge',
          'percentile_C.merge',
          'percentile_T_ExtendMerge.merge',
          'percentile_T.merge.SE_overlapped',
          'percentile_C.merge.SE_overlapped',
          'percentile_T_ExtendMerge.merge.SE_overlapped',]

# flags = ['percentile_T.merge',
#          'percentile_T.merge.SE_overlapped',]


thre = None
# for cancertype in name_match.index[:]:
for halflife in halflifes[:thre]:
    for flag in flags[:thre]:
        print(halflife,flag)
        for cancertype in ['BRCA', 'CESC', 'COAD','LIHC', 'PRAD'][:thre]:
            # os.makedirs(outdir+os.sep+cancertype,exist_ok=True)
            ct = name_match.loc[cancertype,'SE']
            # filtered_ids = filtered_df[filtered_df.cohort==cancertype].index
            rank_df = pd.read_csv('{}/_CP_TFBS_all_vs_TFMS_{}.csv'.format(rank_dir,ct),index_col=0)
            factors = rank_df.index
            # factors = tfbs_cp_s[ct].dropna().index
            
            mutation_file = '{}/{}_mutation.bed'.format(mutation_dir,cancertype)
            mutation_df = pd.read_csv(mutation_file,sep='\t',header=None)
            mutation_df.columns = ['chr','start','end','ID','gene','strand']
            
            outdf = pd.DataFrame()
            for factor in factors[:]:
                # print(halflife, flag, cancertype, factor)
                df = pd.read_csv('{}/{}/{}_avg_RP_halflife_{}_on_{}_{}_{}.csv'.format(indir,cancertype,cancertype,halflife,ct,factor,flag),index_col=0,header=None)    
                # df.mean(axis=0).to_csv('{}/{}/{}_avg_RP_halflife_{}_on_{}_{}_{}.csv'.format(outdir,cancertype,cancertype,halflife,ct,factor,flag),header=None)   
                low_RP_patients = df.sort_values(by=1)[:int(df.shape[0]*.5)].index
                high_RP_patients = df.sort_values(by=1)[-int(df.shape[0]*.5):].index
                 
                low_RP_df = mutation_df[(mutation_df['ID'].isin(low_RP_patients))&
                            (mutation_df['gene']==factor)]
                high_RP_df = mutation_df[(mutation_df['ID'].isin(high_RP_patients))&
                            (mutation_df['gene']==factor)]
                
                outdf.loc[factor,'# low RP total'] = len(low_RP_patients)
                outdf.loc[factor,'# low RP w/ mutations'] = low_RP_df.shape[0]
                outdf.loc[factor,'% low RP w/ mutations'] = low_RP_df.shape[0]/len(low_RP_patients)
                
                outdf.loc[factor,'# high RP total'] = len(high_RP_patients)
                outdf.loc[factor,'# high RP w/ mutations'] = high_RP_df.shape[0]
                outdf.loc[factor,'% high RP w/ mutations'] = high_RP_df.shape[0]/len(high_RP_patients)
                
                s,p = stats.fisher_exact([[high_RP_df.shape[0],len(high_RP_patients)-high_RP_df.shape[0]],
                                          [low_RP_df.shape[0],len(low_RP_patients)-low_RP_df.shape[0]]])
                
                outdf.loc[factor,'fisher-exact-s'] = s
                outdf.loc[factor,'fisher-exact-p'] = p
            
            outdf.to_csv('{}/{}_mutation_by_RP_halflife_{}_{}.csv'.format(outdir,cancertype,halflife,flag))
                # print()




import sys,argparse
import os,glob
import numpy as np
import pandas as pd
from scipy import stats
import re,bisect

chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',\
             'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',\
             'chr18','chr19','chr20','chr21','chr22','chrX','chrY']


indir = 'f5_pancancer_mutation_ByRate_on_Motif_Site'  
outdir = 'f6_pancancer_mutation_data_collection'
os.makedirs(outdir,exist_ok=True)
            
            
cts = ['MCF-7','HCT-116','HeLa','LNCaP','U87','HepG2']
datadir = '../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged/'
tfbs_dir = '../f1_TF_cluster_potential/f2_cor_CP_SE_AICAP/f3_TFBS_CP_heatmap//_csv/'
tfbs_cp_s = pd.read_csv('{}/data_fisher_CP_TFBS_vs_TFMS_CP_RankSum_statistics.csv'.format(tfbs_dir),index_col=0)

# mutation_file='data/pancancer_mutation.bed'
# df = pd.DataFrame()
for ct in cts[:]:
    # os.makedirs(outdir+os.sep+ct,exist_ok=True)
    factors = tfbs_cp_s[ct].dropna().index
    outdf = pd.DataFrame()
    for factor in factors[:]:
        for flag in ['percentile_T','percentile_T_ExtendMerge','percentile_C']:
            mergefile = '{}/{}/{}_{}_{}.merge.bed'.format(datadir,ct,ct,factor,flag)
            if os.path.isfile(mergefile):
                outfile = '{}/{}/{}_{}_{}_mutatioinCount.bed'.format(indir,ct,ct,factor,flag)

                df = pd.read_csv(outfile,sep='\t')
                # df.columns = ['chr','start','end','id','count']
                total_len = (df['end']-df['start']).sum()
                total_mutation_count = df['mutationCount'].sum()
                total_mutation_rate = df['mutationRate'].sum()*6285
                # df['mutation rate per bp'] = df['count']/(df['end']-df['start'])
                outdf.loc[factor,'{} total_len'.format(flag)] = total_len
                outdf.loc[factor,'{} total mutationCount'.format(flag)] = total_mutation_count
                outdf.loc[factor,'{} %mutationCount per bp'.format(flag)] = total_mutation_count/total_len
                outdf.loc[factor,'{} total mutationRate'.format(flag)] = total_mutation_rate
                outdf.loc[factor,'{} %mutationRate per bp'.format(flag)] = total_mutation_rate/total_len
                
    outdf.to_csv(outdir+os.sep+'{}_murationRate.csv'.format(ct))










import os,sys,argparse
# import fileinput,time
import glob
import re,bisect
import pandas as pd
import numpy as np
import matplotlib
from scipy import stats




indir = '../f3_clustered_TFBS/f2_bedfiles_merged'
merge_indir = '../f3_clustered_TFBS/f3_coBinding_merge/'
outdir = 'f3_coBinding_CTCF_specificity'
os.makedirs(outdir,exist_ok=True)



selected_factors = {}
tfbs_cp_dir = '../f2_cor_CP_SE_AICAP/f9_per_CT_TFBS_CP_cor_zscore_CP/TFBS_CP/'
for ct in ['MCF-7','HCT-116'][:]:
    df = pd.read_csv('{}/_CP_TFBS_all_vs_TFMS_{}.csv'.format(tfbs_cp_dir,ct),index_col=0)
    selected_factors['{} top_TFBSCP'.format(ct)] = df['TFBS CP rank'].sort_values().iloc[:3].index
    selected_factors['{} top_zscored_TFBSCP'.format(ct)] = df['avg rank'].sort_values().iloc[:3].index



for treat_flag in ['percentile_T','percentile_T_ExtendMerge'][:1]:
    for ct in ['MCF-7','HCT-116'][:]:
        for factorType in ['top_TFBSCP','top_zscored_TFBSCP'][1:]:
            factors = selected_factors['{} {}'.format(ct,factorType)]
            # print(factors)
            # bedfiles = ['{}/{}/{}_{}_{}.merge.bed'.format(indir,ct,ct,i,treat_flag) for i in factors]
            subdir = '{}_{}'.format(ct,factorType)
            os.makedirs(outdir+os.sep+subdir,exist_ok=True)
            
            mergeFile = '{}/{}/{}.merge.bed'.format(merge_indir,subdir,treat_flag)
            
            
            for factor in factors:
                afile = mergeFile
                bfile = '{}/{}/{}_{}_{}.merge.bed'.format(indir,ct,ct,factor,treat_flag)
                outfile = '{}/{}/{}.merge.{}.overlapped.bed'.format(outdir,subdir,treat_flag,factor)
                commandLine = 'bedtools intersect -a {} -b {} -wa -c > {}'.format(afile,bfile,outfile)
                print(commandLine)
                
            for prename in ['CBS', 'CBS_occupancy_GT3', 
                            'CBS_occupancy_GT3_thre_0.1',
                            'CBS_occupancy_GT3_thre_0.2',
                            'CBS_occupancy_GT3_thre_0.3',
                            'CBS_thre_0.1',
                            'CBS_thre_0.2',
                            'CBS_thre_0.3',
                            'CBS_thre_0.8']:
                afile = mergeFile
                bfile = 'data/union_{}.bed'.format(prename)
                outfile = '{}/{}/{}.merge.{}.overlapped.bed'.format(outdir,subdir,treat_flag,prename)
                commandLine = 'bedtools intersect -a {} -b {} -wa -c > {}'.format(afile,bfile,outfile)
                print(commandLine)
            
        
        






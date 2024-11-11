import os,sys,argparse
# import fileinput,time
import glob
import re,bisect
import pandas as pd
import numpy as np
import matplotlib
from scipy import stats

def merge_files(bedfiles,outdir,ct,factor,flag,excluded=False):

    print(ct,factor)
    se_file = '../../data/SE_hg38/{}.bed'.format(ct)
    catfile = '{}/{}/{}_{}_{}.cat.bed'.format(outdir,ct,ct,factor,flag)
    sortfile = '{}/{}/{}_{}_{}.sort.bed'.format(outdir,ct,ct,factor,flag)
    mergefile = '{}/{}/{}_{}_{}.merge.bed'.format(outdir,ct,ct,factor,flag)
    mergefile_SE = '{}/{}/{}_{}_{}.merge.SE_overlapped.bed'.format(outdir,ct,ct,factor,flag)
    commandLine = 'cat {} > {}'.format(' \\\n'.join(bedfiles),catfile)
    print(commandLine);os.system(commandLine)
    commandLine = 'bedtools sort -i {} > {}'.format(catfile,sortfile)
    print(commandLine);os.system(commandLine)
    commandLine = 'bedtools merge -i {} > {}'.format(sortfile,mergefile)
    print(commandLine);os.system(commandLine)
    
    if excluded:
        mergefile_tmp = '{}/{}/{}_{}_{}.merge.tmp.bed'.format(outdir,ct,ct,factor,flag)
        commandLine = 'bedtools intersect -a {} -b {} -wa -v > {}'.format(mergefile,excluded,mergefile_tmp)
        print(commandLine);os.system(commandLine)
        commandLine = 'mv {} {}'.format(mergefile_tmp,mergefile)
        print(commandLine);os.system(commandLine)

    commandLine = 'bedtools intersect -a {} -b {} -wa -u > {}\n'.format(mergefile,se_file,mergefile_SE)
    print(commandLine);os.system(commandLine)
    commandLine = 'rm {} {}'.format(catfile,sortfile)
    print(commandLine);os.system(commandLine)
    
    return mergefile






indir = 'f1_clustered_TFBS_percentile-5'
outdir = 'f2_bedfiles_merged'
os.makedirs(outdir,exist_ok=True)

# name_match = pd.read_excel('../../data/TCGA/TCGA-ATAC_SE_cancerType_match.xlsx',index_col=0)   
# name_match = name_match.dropna()
   

cts = ['MCF-7','HCT-116','HeLa','LNCaP','U87','HepG2']
tfbs_dir = '..//f2_cor_CP_SE_AICAP/f3_TFBS_CP_heatmap//_csv/'
tfbs_cp_s = pd.read_csv('{}/data_fisher_CP_TFBS_all_vs_TFMS_CP_KStest_statistics.csv'.format(tfbs_dir),index_col=0)
    

for ct in cts[:]:
    factors = tfbs_cp_s[ct].dropna().index
    os.makedirs(outdir+os.sep+ct,exist_ok=True)
    
    for factor in factors[:]:
        bedfiles_T = glob.glob('{}/{}/{}_{}_*_T.bed'.format(indir,ct,ct,factor))
        bedfiles_T_merge = glob.glob('{}/{}/{}_{}_*_T_ExtendMerge.bed'.format(indir,ct,ct,factor))
        bedfiles_C = glob.glob('{}/{}/{}_{}_*_C.bed'.format(indir,ct,ct,factor))
        
        if len(bedfiles_T)>0:
            mergefile_T = merge_files(bedfiles_T,outdir,ct,factor,'percentile_T')
            mergefile_E = merge_files(bedfiles_T_merge,outdir,ct,factor,'percentile_T_ExtendMerge')
            mergefile_C = merge_files(bedfiles_C,outdir,ct,factor,'percentile_C',excluded=mergefile_T)
        







import os,sys,argparse
# import fileinput,time
import glob
import re,bisect,os
import pandas as pd
import numpy as np
import matplotlib
from scipy import stats



outdir = 'f0_TFBS_merge'
os.makedirs(outdir,exist_ok=True)


# data_dir = '/Volumes/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor'
# data_dir = '/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor'
data_dir = '../../data/cistrome/peaks_all/'
tfbs_dir = '../f2_cor_CP_SE_AICAP/f3_TFBS_CP_heatmap/_csv'   
# tfbs_cp_types = ['TFBS','TFBS_overlap_motif','TFBS_NOT_overlap_motif']
tfbs_cp_types = ['CP_TFBS_all_vs_TFMS',
                 'CP_TFBS_overlap_motif_vs_TFMS',
                 'CP_TFBS_NOT_overlap_motif_vs_TFMS']



for tfbs_cp_type in tfbs_cp_types[:1]:
    tfbs_cp_s = pd.read_csv('{}/data_fisher_{}_CP_RankSum_statistics.csv'.format(tfbs_dir,tfbs_cp_type),index_col=0)
    df = pd.read_csv('../f1_TFMS_TFBS_CP/{}/per_data_{}.csv'.format(tfbs_cp_type,tfbs_cp_type),index_col=0)
    # df = df[(df['#TFBS']>2000)]
    
    for ct in ['MCF-7','HCT-116','HeLa','LNCaP','U87','HepG2'][:]:
        outdir_ct = '{}/{}'.format(outdir,ct)
        os.makedirs(outdir_ct,exist_ok=True)
        factors = tfbs_cp_s[ct].dropna().index
        
        for factor in factors[:]:            
            prename = '{}_{}_'.format(ct,factor)            
            ids = [i for i in df.index if re.search(prename,i)]
            # ids = [i.split('_')[2] for i in ids]
            peak_files = ['{}/{}.bed'.format(data_dir,i) for i in ids]
            print(ct,factor,len(ids))
            
            cat_file = '{}/{}_{}.cat.bed'.format(outdir_ct,ct,factor)
            sort_file = '{}/{}_{}.sort.bed'.format(outdir_ct,ct,factor)
            merge_file = '{}/{}_{}.merge.bed'.format(outdir_ct,ct,factor)
            se_file = '../../data/SE_hg38/{}.bed'.format(ct)
            se_overlapped_file = '{}/{}_{}.merge.SE_overlapped.bed'.format(outdir_ct,ct,factor)
            
            commandLine = 'cat {} > {}\n'.format(' \\\n'.join(peak_files),cat_file)
            print(commandLine);os.system(commandLine)
            commandLine = 'bedtools sort -i {} > {}'.format(cat_file,sort_file)
            print(commandLine);os.system(commandLine)
            commandLine = 'bedtools merge -i {} > {}'.format(sort_file,merge_file)
            print(commandLine);os.system(commandLine)
            commandLine = 'bedtools intersect -a {} -b {} -wa -u > {}\n'.format(merge_file,se_file,se_overlapped_file)
            print(commandLine);os.system(commandLine)
    









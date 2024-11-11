import os,sys,argparse
# import fileinput,time
import glob
import re,bisect
import pandas as pd
import numpy as np
import matplotlib
from scipy import stats



def get_lines(infile):

    with open(infile,'rb') as f:
        lines = 0
        buf_size = 1024*1024
        buf = f.raw.read(buf_size)
        while buf:
            lines += buf.count(b'\n')
            buf = f.raw.read(buf_size)
    return lines




outdir = 'f2_bedfiles_merged'
os.makedirs(outdir,exist_ok=True)
 
# tfbs_prenames = ['TFBS','TFBS_overlap_motif','TFBS_NOT_overlap_motif']
# percentile_thres = ['percentile-5','percentile-1']
cts = ['MCF-7','HCT-116','HeLa','LNCaP','U87','HepG2']
tfbs_dir = '../../f1_TF_cluster_potential/f2_cor_CP_SE_AICAP/f3_TFBS_CP_heatmap//_csv/'
tfbs_cp_s = pd.read_csv('{}/data_fisher_CP_TFBS_all_vs_TFMS_CP_KStest_statistics.csv'.format(tfbs_dir),index_col=0)


df = pd.DataFrame()
for ct in cts[:]:
    factors = tfbs_cp_s[ct].dropna().index
    for factor in factors[:]:
        for flag in ['percentile_T','percentile_T_ExtendMerge','percentile_C']:
            mergefile = '{}/{}/{}_{}_{}.merge.bed'.format(outdir,ct,ct,factor,flag)
            mergefile_SE = '{}/{}/{}_{}_{}.merge.SE_overlapped.bed'.format(outdir,ct,ct,factor,flag)

            if os.path.isfile(mergefile):
                index = '{} {}'.format(ct, factor)
                a = get_lines(mergefile)
                b = get_lines(mergefile_SE)
                df.loc[index,'# {}'.format(flag)] = a
                df.loc[index,'# {} on SE'.format(flag)] = b
                df.loc[index,'% {} on SE'.format(flag)] = b/a
            
df.to_csv('{}/data_merged_SE_overlapped.csv'.format(outdir))
    








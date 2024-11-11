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




outdir = 'f2_CI_intersect_TFBS'
os.makedirs(outdir,exist_ok=True)
 
tfbs_dir = '../../f1_TF_cluster_potential/f2_cor_CP_SE_AICAP/f3_TFBS_CP_heatmap_with_motif_SE//_csv/'
tfbs_cp_s = pd.read_csv('{}/data_fisher_CP_TFBS_all_vs_TFMS_CP_KStest_statistics.csv'.format(tfbs_dir),index_col=0)
merge_file_dir = '../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged'

cts = ['MCF-7','HCT-116','HeLa','LNCaP','U87','HepG2']
cellTypes = ['LNCaP','H1','HepG2','MCF7','Panc1','GM12878','HeLa','IMR90','K562','HCT116','Jurkat']
genomic_KB_distances = [20,50,100,200,500]

ci_dir = '../../../f11_TF_condensates_KS_test/f3_public_data/f2_TFBS_CI/f1_CI/'

for ct in cts[:]:
    ct_rename = ''.join(ct.split('-'))
    
    
    factors = tfbs_cp_s[ct].dropna().index
    for factor in factors[:]:
        for flag in ['percentile_T','percentile_C','percentile_T_ExtendMerge'][:]:
            mergefile = '{}/{}/{}_{}_{}.merge.bed'.format(merge_file_dir,ct,ct,factor,flag)
            if os.path.isfile(mergefile):
                
                for genomic_KB_distance in genomic_KB_distances:
                    ci_file = '{}/{}_CI_{}KB.bed'.format(ci_dir,ct_rename,genomic_KB_distance)
                    if os.path.isfile(ci_file):
                        
                        os.makedirs(outdir+os.sep+ct_rename,exist_ok=True)
                        outfile = '{}/{}/{}_{}_{}_overlap_CI_{}KB.bed'.format(outdir,ct_rename,ct_rename,factor,flag,genomic_KB_distance)
                        commandLine = 'bedtools intersect \\\n-a {} \\\n-b {} \\\n-wa -u > {}\n'.format(ci_file,mergefile,outfile)
                        print(commandLine);os.system(commandLine)







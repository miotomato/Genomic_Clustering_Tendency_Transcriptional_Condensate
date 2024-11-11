import os,sys,argparse
# import fileinput,time
import glob
import re,bisect,os
import pandas as pd
import numpy as np
import matplotlib
from scipy import stats



outdir = 'f8_diff_ATAC_overlap_TFBS_clustered'
os.makedirs(outdir,exist_ok=True)

tfbs_dir = '../../f1_TF_cluster_potential/f2_cor_CP_SE_AICAP/f3_TFBS_CP_heatmap/_csv'   
# tfbs_cp_types = ['CP_TFBS_vs_TFMS',
#                  'CP_TFBS_overlap_motif_vs_TFMS',
#                  'CP_TFBS_NOT_overlap_motif_vs_TFMS']

atac_file = '../../../f9_TF_condensates_V3/data/TCGA/tcga_atac.bed'
mergefile_dir = '../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged/'
# percentile_mergefile_dir = '../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merge/'



# for tfbs_cp_type in tfbs_cp_types[:1]:
if 1:
    tfbs_cp_s = pd.read_csv('{}/data_fisher_CP_TFBS_all_vs_TFMS_CP_KStest_statistics.csv'.format(tfbs_dir),index_col=0)
    
    for ct in ['MCF-7','HCT-116','HeLa','LNCaP','U87','HepG2'][:]:
        factors = tfbs_cp_s[ct].dropna().sort_values(ascending=False).index
        
        for factor in factors[:]:
            for flag in ['percentile_T','percentile_T_ExtendMerge','percentile_C']:
                mergefile = '{}/{}/{}_{}_{}.merge.bed'.format(mergefile_dir,ct,ct,factor,flag)
                                  
                if os.path.isfile(mergefile):                    
                    outfile = '{}/atac_overlap_{}_{}_{}.bed'.format(outdir,ct,factor,flag)
                    commandLine = 'bedtools intersect -a {} -b {} -wa -c > {}'.format(atac_file,mergefile,outfile)
                    print(commandLine)
                    os.system(commandLine)
    

        se_file = '../../data/SE_hg38/{}.bed'.format(ct)
        outfile = '{}/atac_overlap_{}_SE.bed'.format(outdir,ct)
        commandLine = 'bedtools intersect -a {} -b {} -wa -c > {}\n'.format(atac_file,se_file,outfile)
        print(commandLine)
        os.system(commandLine)








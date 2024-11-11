import os,sys,argparse
# import fileinput,time
import glob
import re,bisect,os
import pandas as pd
import numpy as np
import matplotlib
from scipy import stats



outdir = 'f5_diff_ATAC_overlap_TFBS'
os.makedirs(outdir,exist_ok=True)

tfbs_dir = '../../f1_TF_cluster_potential/f2_cor_CP_SE_AICAP/f3_TFBS_CP_heatmap/_csv'   
tfbs_cp_types = ['CP_TFBS_vs_TFMS',
                 'CP_TFBS_overlap_motif_vs_TFMS',
                 'CP_TFBS_NOT_overlap_motif_vs_TFMS']

atac_file = '../../data/TCGA/tcga_atac.bed'
mergefile_dir = '../../f1_TF_cluster_potential/f3_clustered_TFBS/f0_TFBS_merge/'
percentile_mergefile_dir = '../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merge/'



for tfbs_cp_type in tfbs_cp_types[:1]:
    tfbs_cp_s = pd.read_csv('{}/data_fisher_{}_CP_RankSum_statistics.csv'.format(tfbs_dir,tfbs_cp_type),index_col=0)
    
    for ct in ['MCF-7','HCT-116','HeLa','LNCaP','U87','HepG2'][:2]:
        factors = tfbs_cp_s[ct].dropna().sort_values(ascending=False).index
        
        for factor in factors[:]:                        
            merge_file = '{}/{}/{}_{}.merge.bed'.format(mergefile_dir,ct,ct,factor)
            se_overlapped_file = '{}/{}/{}_{}.merge.SE_overlapped.bed'.format(mergefile_dir,ct,ct,factor)
            percentile_merge_file = '{}/{}_{}.merge.bed'.format(percentile_mergefile_dir,ct,factor)
            
            bfile = merge_file
            outfile = '{}/atac_overlap_{}_{}.bed'.format(outdir,ct,factor)
            commandLine = 'bedtools intersect -a {} -b {} -wa -c > {}\n'.format(atac_file,bfile,outfile)
            # print(commandLine)#;os.system(commandLine)
    
            bfile = se_overlapped_file
            outfile = '{}/atac_overlap_{}_{}_overlap_SE.bed'.format(outdir,ct,factor)
            commandLine = 'bedtools intersect -a {} -b {} -wa -c > {}\n'.format(atac_file,bfile,outfile)
            # print(commandLine)#;os.system(commandLine)
    
            bfile = percentile_merge_file
            if os.path.isfile(bfile):
                outfile = '{}/atac_overlap_{}_{}_percentile5.bed'.format(outdir,ct,factor)
                commandLine = 'bedtools intersect -a {} -b {} -wa -c > {}\n'.format(atac_file,bfile,outfile)
                print(commandLine);os.system(commandLine)
    

        se_file = '../../data/SE_hg38/{}.bed'.format(ct)
        outfile = '{}/atac_overlap_{}_SE.bed'.format(outdir,ct)
        commandLine = 'bedtools intersect -a {} -b {} -wa -c > {}\n'.format(atac_file,se_file,outfile)
        # print(commandLine)#;os.system(commandLine)








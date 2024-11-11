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


    
tfbs_prenames = ['TFBS_all','TFBS_overlap_motif','TFBS_NOT_overlap_motif']
percentile_thres = ['percentile-5','percentile-1']
cts = ['MCF-7','HCT-116','HeLa','LNCaP','U87','HepG2']


for percentile_thre in percentile_thres[:1]:
    outdir = 'f1_clustered_TFBS_{}'.format(percentile_thre)
    
    for tfbs_prename in tfbs_prenames[:1]:
        cp_prename = 'CP_{}_vs_TFMS'.format(tfbs_prename)
        df = pd.read_csv('{}/{}_new.csv'.format(outdir,cp_prename),index_col=0)
        
        for index in df.index[:]:
            data_id = index.split('_')[2]
            tf = index.split('_')[1]
            ct = index.split('_')[0]
            if ct not in cts:
                continue
            
            outname = '{}_{}_{}'.format(ct,tf,data_id)
            tfbs_prename_old = 'TFBS' if tfbs_prename == 'TFBS_all' else tfbs_prename
            data_file = '../../../f8_TF_condensates_V2/f1_TF_cluster_potential/f1_bedtools_closest//data_{}/{}_sort_peaks.narrowPeak.tsv'.format(tfbs_prename_old,data_id)
            se_file = '../../data/SE_hg38/{}.bed'.format(ct)
            bed_T_file = '{}/{}/{}_T.bed'.format(outdir,ct,outname)
            bed_T_SE = '{}/{}/{}_T_on_SE.bed'.format(outdir,ct,outname)
            bed_C_file = '{}/{}/{}_C.bed'.format(outdir,ct,outname)
            bed_C_SE = '{}/{}/{}_C_on_SE.bed'.format(outdir,ct,outname)
            merge_file = '{}/{}/{}_T_ExtendMerge.bed'.format(outdir,ct,outname)
            merge_SE = '{}/{}/{}_T_ExtendMerge_on_SE.bed'.format(outdir,ct,outname)
            
            a = get_lines(bed_T_file)
            b = get_lines(bed_T_SE)
            df.loc[index,'#TFBS <= {}'.format(percentile_thre)] = a
            df.loc[index,'#TFBS <= {} on SE'.format(percentile_thre)] = b
            df.loc[index,'%TFBS <= {} on SE'.format(percentile_thre)] = b/a
            
            a = get_lines(merge_file)
            b = get_lines(merge_SE)
            df.loc[index,'#TFBS <= {} merged'.format(percentile_thre)] = a
            df.loc[index,'#TFBS <= {} merged on SE'.format(percentile_thre)] = b
            df.loc[index,'%TFBS <= {} merged on SE'.format(percentile_thre)] = b/a
            
            a = get_lines(bed_C_file)
            b = get_lines(bed_C_SE)
            df.loc[index,'#TFBS > {}'.format(percentile_thre)] = a
            df.loc[index,'#TFBS > {} on SE'.format(percentile_thre)] = b
            df.loc[index,'%TFBS > {} on SE'.format(percentile_thre)] = b/a
    
        df.to_csv('{}/{}_new_add_data.csv'.format(outdir,cp_prename))
    








import os,sys,argparse
# import fileinput,time
import glob
import re,bisect
import pandas as pd
import numpy as np
import subprocess
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


def run_bedtools_intersect(afile,bfile):
    cl = 'bedtools intersect -a {} -b {} -u -wa|wc -l'.format(afile,bfile)
    # print(cl)
    wa_count = subprocess.check_output(cl,shell=True).decode(sys.stdout.encoding).strip()
    return int(wa_count)
    


# ==== main
    
outdir_par = 'f5_cellType_specific_peaks_cor_SE_over_AllPeaks'
os.makedirs(outdir_par,exist_ok=True)

list_df = pd.read_csv('f2_selected_data/pairwise_sharedTF_count_id.csv',index_col=0)
for celltype_pairs in list_df.index[:]:
    # data for each cell type pairs
    print(celltype_pairs)
    outdir = outdir_par+os.sep+celltype_pairs
    os.makedirs(outdir,exist_ok=True)
    data_df = pd.read_csv('f2_selected_data/selected_data_IDs_{}.csv'.format(celltype_pairs),index_col=0)   

    df = pd.DataFrame()
    df_s = pd.DataFrame()
    df_p = pd.DataFrame()
    
    for factor in data_df.index[:]:
        factor_df = data_df.loc[factor].dropna()#;print(factor)
        for celltype in factor_df.index[:]:
            celltype_specific_peak_file='f4_cellType_specific_peaks/{}/{}_{}_specific.bed'.format(celltype_pairs,factor,celltype)
            celltype_all_peak_file = 'f3_union_peaks_across_cellTypes/{}/{}_merge_intersect_{}.bed'.format(celltype_pairs,factor,celltype)
            # celltype_all_peak_file = 'f3_union_peaks_across_cellTypes/{}/_{}_{}.bed'.format(celltype_pairs,factor,celltype)
            celltype_SE_file = '../f1_TF_motif_cluster_cor_SE/f4_cluster_enrichment_at_SE/data/SE_hg38/{}.bed'.format(celltype)
            if os.path.isfile(celltype_SE_file) and os.path.isfile(celltype_specific_peak_file):            
                # == whether the motif overlap with SE
                all_total = get_lines(celltype_all_peak_file)
                all_overlap = run_bedtools_intersect(celltype_all_peak_file,celltype_SE_file)
                # == whether the cluster motif overlap with SE
                sub_total = get_lines(celltype_specific_peak_file)
                sub_overlap = run_bedtools_intersect(celltype_specific_peak_file,celltype_SE_file)
                # == fisher exact test
                a = sub_overlap
                b = sub_total-sub_overlap
                c = all_overlap - sub_overlap
                d = all_total - sub_total - all_overlap + sub_overlap
                s,p = stats.fisher_exact([[a, b], [c, d]])
                                          
                df.loc['{}_{}'.format(factor,celltype),'all_peaks_total'] = all_total
                df.loc['{}_{}'.format(factor,celltype),'all_peaks_overlap_SE'] = all_overlap
                df.loc['{}_{}'.format(factor,celltype),'celltype_specific_peaks_total'] = sub_total
                df.loc['{}_{}'.format(factor,celltype),'celltype_specific_peaks_overlap_SE'] = sub_overlap
                df.loc['{}_{}'.format(factor,celltype),'fisher_exact_s'] = s
                df.loc['{}_{}'.format(factor,celltype),'fisher_exact_p'] = p
                df_s.loc[factor,celltype] = s
                df_p.loc[factor,celltype] = p
            
    # df.index.name='TF'
    df.to_csv('{}/fisher_details.csv'.format(outdir))
    df_s.to_csv('{}/fisher_s.csv'.format(outdir))
    df_p.to_csv('{}/fisher_p.csv'.format(outdir))
        
                              
                          


import os,sys,argparse,glob
import numpy as np
import pandas as pd
# import find_overlap_keep_info_NOT_sep_strand_asimport




# == main 
indir = 'f2_diff_ATAC_each_CancerType'
outdir = 'f3_diff_ATAC_overlap_SE'
os.makedirs(outdir,exist_ok=True)

# read matched names between TCGA ATAC and SE data
name_match = pd.read_excel('../data/TCGA/TCGA-ATAC_SE_cancerType_match.xlsx',index_col=0)   
name_match = name_match.dropna()
pyfile="find_overlap_keep_info_NOT_sep_strand_revised.py"

for cancertype in name_match.index:
    cancertype_SE = name_match.loc[cancertype].SE
    cancertype_SE_rename = name_match.loc[cancertype].SE_rename
    tcga_file = indir+os.sep+'{}_diff_atac.csv'.format(cancertype)
    se_file = '../data/SE_hg38/{}.bed'.format(cancertype_SE)
    # == get TCGA peaks overlap/NOT-overlap with SE data
    overlap_file = '{}/{}_ATAC_overlap_{}_SE.bed'.format(outdir,cancertype,cancertype_SE_rename)
    nonoverlap_file = '{}/{}_ATAC_NOT_overlap_{}_SE.bed'.format(outdir,cancertype,cancertype_SE_rename)
    print('python {} -a {} -b {} -s hg38 -p {} -q {} -e2 0'.format(pyfile,tcga_file,se_file,overlap_file,nonoverlap_file))
    os.system('python {} -a {} -b {} -s hg38 -p {} -q {} -e2 0'.format(pyfile,tcga_file,se_file,overlap_file,nonoverlap_file))


  


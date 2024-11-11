import os,sys,argparse,glob
import numpy as np
import pandas as pd
# import find_overlap_keep_info_NOT_sep_strand_asimport




# == main 
outdir = 'f1_ATAC_overlap_SE_caseID'
os.makedirs(outdir,exist_ok=True)

# read matched names between TCGA ATAC and SE data
name_match = pd.read_excel('../../data/TCGA/TCGA-ATAC_SE_cancerType_match.xlsx',index_col=0)   
name_match = name_match.dropna()
pyfile="find_overlap_keep_info_NOT_sep_strand_revised.py"
# tcga_file = '../../data/TCGA/mynorm_TCGA-ATAC_PanCan_Log2_QuantileNorm_Counts_plus5.caseID.avg.head5k.txt'
tcga_file = '../../data/TCGA/mynorm_TCGA-ATAC_PanCan_Log2_QuantileNorm_Counts_plus5.caseID.avg.txt'

for cancertype in name_match.index[:]:
    cancertype_SE = name_match.loc[cancertype].SE
    cancertype_SE_rename = name_match.loc[cancertype].SE_rename
    se_file = '../../data/SE_hg38/{}.bed'.format(cancertype_SE)
    # == get TCGA peaks overlap/NOT-overlap with SE data
    overlap_file = '{}/{}_ATAC_overlap_{}_SE_caseID.bed'.format(outdir,cancertype,cancertype_SE_rename)
    print('python {} -a {} -b {} -s hg38 -p {} -e2 0'.format(pyfile,tcga_file,se_file,overlap_file))
    os.system('python {} -a {} -b {} -s hg38 -p {} -e2 0'.format(pyfile,tcga_file,se_file,overlap_file))


  


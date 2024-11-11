import os,sys,argparse,glob
import numpy as np
import pandas as pd
# import find_overlap_keep_info_NOT_sep_strand_asimport
import re
from scipy import stats
# import matplotlib
# # matplotlib.use('Agg')
# import matplotlib.pyplot as plt
# matplotlib.rcParams['font.size']=16
# import seaborn as sns
# sns.set(font_scale=1.2)
# sns.set_style("whitegrid", {'axes.grid' : False,'grid.color': 'grey'})
# sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k'})
# matplotlib.rcParams["font.sans-serif"] = ["Arial"]
# matplotlib.rcParams['mathtext.fontset']='custom'
# matplotlib.rcParams['mathtext.rm']="Arial"
# import CTCF_TALL_modules_new



# main 
indir = 'f2_diff_ATAC_each_CancerType'
outdir = 'f5_diff_ATAC_bedfile'
os.makedirs(outdir,exist_ok=True)



# read cancer types
name_match = pd.read_excel('../data/TCGA/TCGA-ATAC_SE_cancerType_match.xlsx',index_col=0)   
name_match = name_match.dropna()
filtered_df = pd.read_excel('../data/TCGA/TCGA-ATAC_clustered_samples.xlsx',index_col=0)   

# == differential atac for each cancer type
for cancertype in name_match.index[:]:
    cancertype_SE = name_match.loc[cancertype].SE
    cancertype_SE_rename = name_match.loc[cancertype].SE_rename
    filtered_id = filtered_df[filtered_df.cohort==cancertype].case_id
    filtered_id = ['{}_{}'.format(cancertype,i) for i in filtered_id]
    
    # == cancer data
    cancer_df = pd.read_csv(indir+os.sep+'{}_diff_atac.csv'.format(cancertype),sep='\t')
    cancer_df = cancer_df[['#seqnames', 'start', 'end', 'name', 'stats', 'pvalue']]
    cancer_df.to_csv(outdir+os.sep+'{}_diff_atac.bed'.format(cancertype),sep='\t',header=None,index=False)
    
    
 


 



 
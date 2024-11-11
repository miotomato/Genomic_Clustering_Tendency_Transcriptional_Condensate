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
outdir = 'f2_diff_ATAC_each_CancerType'
os.makedirs(outdir,exist_ok=True)


# == read atac data                 
# infile = '../data/TCGA/mynorm_TCGA-ATAC_PanCan_Log2_QuantileNorm_Counts_plus5.caseID.avg.head5k.txt'
infile = '../data/TCGA/mynorm_TCGA-ATAC_PanCan_Log2_QuantileNorm_Counts_plus5.caseID.avg.txt'
df = pd.read_csv(infile,sep='\t')
df_info = df.iloc[:,:7]
df_data = df.iloc[:,7:]

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
    cancer_col = df_data.columns.intersection(filtered_id);print(cancertype,len(cancer_col))
    cancer_data = df_data[cancer_col]
    # == control data
    control_col = df_data.columns.difference(cancer_col)
    control_data = df_data[control_col]
    s,p = stats.ttest_ind(cancer_data,control_data,axis=1)
    cancer_df = pd.concat([df_info,pd.DataFrame(s,columns = ['stats'])],axis=1)
    cancer_df = pd.concat([cancer_df,pd.DataFrame(p,columns = ['pvalue'])],axis=1)
    cancer_df.to_csv(outdir+os.sep+'{}_diff_atac.csv'.format(cancertype),index=False,sep='\t')
    
    
 


 



 
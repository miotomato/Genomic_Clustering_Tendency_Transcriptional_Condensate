import os,sys,argparse,glob
import numpy as np
import pandas as pd
# import find_overlap_keep_info_NOT_sep_strand_asimport
import re
from scipy import stats
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=16
import seaborn as sns
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid' : False,'grid.color': 'grey'})
sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k'})
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams['mathtext.fontset']='custom'
matplotlib.rcParams['mathtext.rm']="Arial"
# import CTCF_TALL_modules_new



# main 
infile = 'mynorm_TCGA-ATAC_PanCan_all_SE_overlapped.head5k.bed'
infile = 'mynorm_TCGA-ATAC_PanCan_all_SE_overlapped.bed'
outdir = 'f2_diff_ATAC_each_CancerType'
os.makedirs(outdir,exist_ok=True)

# == read atac data                 
df = pd.read_csv(infile,sep='\t')
df_info = df.iloc[:,:7]
df_data = df.iloc[:,7:]

# == differential atac for each cancer type
cancertypes =  set([i.split('_')[0] for i in df_data.columns])
for cancertype in cancertypes:
    cancer_col = [i for i in df_data.columns if i.startswith(cancertype)]
    cancer_data = df_data[cancer_col]
    # == control data
    control_col = df_data.columns.difference(cancer_col)
    control_data = df_data[control_col]
    s,p = stats.ttest_ind(cancer_data,control_data,axis=1)
    cancer_df = pd.concat([df_info,pd.DataFrame(s,columns = ['stats'])],axis=1)
    cancer_df = pd.concat([cancer_df,pd.DataFrame(p,columns = ['pvalue'])],axis=1)
    cancer_df.to_csv(outdir+os.sep+'{}_diff_atac.csv'.format(cancertype),index=False,sep='\t')
    
    
 


 



 
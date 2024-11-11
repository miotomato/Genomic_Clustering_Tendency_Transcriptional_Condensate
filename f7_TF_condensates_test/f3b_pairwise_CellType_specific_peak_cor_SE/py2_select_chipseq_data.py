import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect

# import matplotlib
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt
# matplotlib.rcParams['font.size']=16
# import seaborn as sns
# sns.set(font_scale=1.2)
# sns.set_style("whitegrid", {'axes.grid' : False})
# sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k'})
# matplotlib.rcParams["font.sans-serif"] = ["Arial"]
# matplotlib.rcParams['mathtext.fontset'] = 'custom'
# matplotlib.rcParams["mathtext.rm"] = "Arial"


outdir = 'f2_selected_data'
os.makedirs(outdir,exist_ok=True)

# read the celltype-factor count matrix 
df = pd.read_csv('f1_cellType_Factor_data_Count/cistrome2019_data_details.csv',index_col=0)
df_count = pd.read_csv('f1_cellType_Factor_data_Count/cistrome2019_data_Count.csv',index_col=0)  

# read SE info
se_dir = '../../data/SE_hg38/'
# se_dir = '/Volumes/zanglab/rs8cp/superEnhancer/human_data/hg38_bed/'
se_files = glob.glob('{}/*.bed'.format(se_dir))
se_celltypes = [os.path.basename(i).split('.bed')[0] for i in se_files]

# get cell types with SE data
shared_celltypes = df_count.columns.intersection(se_celltypes)
df_count = df_count[shared_celltypes]
allf=set()
# shared TFs for pair-wise cell types
cols = len(df_count.columns)
out_df = pd.DataFrame()
out_df_id = pd.DataFrame()
for ii in np.arange(cols)[:]:
    for jj in np.arange(ii,cols)[:]:
        # get the cell type name
        col_i = df_count.columns[ii]
        col_j = df_count.columns[jj]
        basename = '{}_vs_{}'.format(col_i,col_j)
        # get the TFs w/ data in both cell types
        pair_df = df_count[[col_i,col_j]]
        pair_df_sum = pair_df.astype(bool).sum(axis=1)
        pair_df_sum = pair_df_sum[pair_df_sum==2]
        out_df.loc[col_i,col_j] = pair_df_sum.shape[0]
        
        if not ii==jj:
            out_df_id.loc[basename,'#datasets'] = pair_df_sum.shape[0]
            # save chip-seq id
            selected_data_id = pd.DataFrame()
            selected_data_GSM = pd.DataFrame()
            for factor in pair_df_sum.index[:]:
                for celltype in [col_i,col_j]:
                    df_sub = df[(df.Factor==factor)&(df.Cell_line==celltype)]
                    if df_sub.shape[0]>0:
                        kept_index = df_sub.sort_values(by='PeaksFoldChangeAbove10',ascending=False).index[0]
                        selected_data_id.loc[factor,celltype]=int(kept_index)
                        gsmID = df_sub.loc[kept_index].GSMID
                        selected_data_GSM.loc[factor,celltype]=gsmID
                allf.add(factor)
            # == save the data
            selected_data_id.to_csv('{}/selected_data_IDs_{}.csv'.format(outdir,basename))           
            selected_data_GSM.to_csv('{}/selected_data_GSMID_{}.csv'.format(outdir,basename))           

       
out_df.to_csv(outdir+os.sep+'pairwise_sharedTF_count.csv')
out_df_id.to_csv(outdir+os.sep+'pairwise_sharedTF_count_id.csv')
    





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


# cistrome_dir = '/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new'
# cistrome_dir = '/Volumes/zanglab/zw5j/data/cistrome/cistrome_db_2019_new'
indir = 'f3_union_peaks_across_cellTypes'
outdir = 'f4_cellType_specific_peaks'
os.makedirs(outdir,exist_ok=True)

data_df = pd.read_csv('f2_selected_data/selected_data_IDs.csv',index_col=0)   
factors = data_df.index

# for each factor, merge the intersection info and get cell type specific regions
for factor in factors[:]:
    factor_df = data_df.loc[factor].dropna()#;print(factor)
    # read the merged chr-start-end region
    merge_df = pd.read_csv('{}/{}.merge.bed'.format(indir,factor),sep='\t',header=None)
    merge_df.columns = ['chr','start','end']
    for celltype in factor_df.index[:]:
        intersect_df = pd.read_csv('{}/{}_merge_intersect_{}.bed'.format(indir,factor,celltype),sep='\t',header=None)       
        intersect_df.columns = ['chr','start','end',celltype]
        merge_df = pd.concat([merge_df,intersect_df[[celltype]]],axis=1)
    # save the merged data
    merge_df.to_csv('{}/{}_union.csv'.format(outdir,factor),index=False)
    
    # get the cell type specific regions
    for celltype in factor_df.index[:]:
        other_celltypes = factor_df.index[:].difference([celltype])
        celltype_df = merge_df[(merge_df[other_celltypes].sum(axis=1)==0)&(merge_df[celltype]!=0)]
        # celltype_df.to_csv('{}/{}_{}_specific.csv'.format(outdir,factor,celltype),index=False)
        celltype_df[['chr','start','end']].to_csv('{}/{}_{}_specific.bed'.format(outdir,factor,celltype),index=False,header=None,sep='\t')
        
        

            





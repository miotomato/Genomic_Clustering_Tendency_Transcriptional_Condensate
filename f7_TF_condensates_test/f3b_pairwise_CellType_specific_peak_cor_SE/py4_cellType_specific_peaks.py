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


def get_lines(infile):
    with open(infile,'rb') as f:
        lines = 0
        buf_size = 1024*1024
        buf = f.raw.read(buf_size)
        while buf:
            lines += buf.count(b'\n')
            buf = f.raw.read(buf_size)
    return lines



indir_par = 'f3_union_peaks_across_cellTypes'
outdir_par = 'f4_cellType_specific_peaks'
os.makedirs(outdir_par,exist_ok=True)

list_df = pd.read_csv('f2_selected_data/pairwise_sharedTF_count_id.csv',index_col=0)
for celltype_pairs in list_df.index[:]:
    # data for each cell type pairs
    print(celltype_pairs)
    indir = indir_par+os.sep+celltype_pairs
    outdir = outdir_par+os.sep+celltype_pairs
    os.makedirs(outdir,exist_ok=True)
    data_df = pd.read_csv('f2_selected_data/selected_data_IDs_{}.csv'.format(celltype_pairs),index_col=0)   
    
    # for each factor, merge the intersection info and get cell type specific regions
    for factor in data_df.index[:]:
        factor_df = data_df.loc[factor].dropna()#;print(factor)
        # read the merged chr-start-end region
        merge_file = '{}/{}.merge.bed'.format(indir,factor)
        if not get_lines(merge_file) >0:
            continue
        print(factor)
        merge_df = pd.read_csv(merge_file,sep='\t',header=None)
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
            
            

            





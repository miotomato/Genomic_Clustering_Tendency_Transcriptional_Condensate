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


cistrome_dir = '/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new'
# cistrome_dir = '/Volumes/zanglab/zw5j/data/cistrome/cistrome_db_2019_new'
outdir = 'f3_union_peaks_across_cellTypes'
os.makedirs(outdir,exist_ok=True)

data_df = pd.read_csv('f2_selected_data/selected_data_IDs.csv',index_col=0)   
factors = data_df.index

for factor in factors:
    factor_df = data_df.loc[factor].dropna()
    factor_ids = factor_df.values.astype(int)
    if factor.startswith('H'):
        factor_dir='{}/human_hm'.format(cistrome_dir)
    else:
        factor_dir='{}/human_factor'.format(cistrome_dir)
    peak_files = ['{}/{}_sort_peaks.narrowPeak.bed'.format(factor_dir,factor_id) for factor_id in factor_ids]
#     for peak_file in peak_files:
#         if not os.path.isfile(peak_file):
#             print(peak_file)
#     continue
    # merge all data from the same factor
    commandline = 'cat {} > {}/{}.cat.bed'.format(' \\\n'.join(peak_files),outdir,factor)
    print(commandline)
    commandline = 'bedtools sort -i {}/{}.cat.bed > {}/{}.sort.bed'.format(outdir,factor,outdir,factor)
    print(commandline)
    commandline = 'bedtools merge -i {}/{}.sort.bed > {}/{}.merge.bed'.format(outdir,factor,outdir,factor)
    print(commandline)
    # bedtools intersect
    for celltype in factor_df.index:
        peak_file = '{}/{}_sort_peaks.narrowPeak.bed'.format(factor_dir,int(factor_df.loc[celltype]))
        commandline = 'bedtools intersect -a {}/{}.merge.bed -b {} -c > {}/{}_merge_intersect_{}.bed'.format(outdir,factor,peak_file,outdir,factor,celltype)
        print(commandline)
    for celltype in factor_df.index:
        peak_file = '{}/{}_sort_peaks.narrowPeak.bed'.format(factor_dir,int(factor_df.loc[celltype]))
        commandline = 'cp {} {}/_{}_{}.bed'.format(peak_file,outdir,factor,celltype)
        print(commandline)
    # add new line to the end per factor    
    print()
        
    
    
            





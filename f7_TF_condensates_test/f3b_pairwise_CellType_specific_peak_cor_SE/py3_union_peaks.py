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


def return_cistrome_peak_file(cistrome_id):
    # cistrome_dir = '/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019'
    for batch in ['human_hm','human_factor']:
        df = pd.read_csv('{}/{}.txt'.format(cistrome_dir,batch),sep='\t',index_col=0)
        if cistrome_id in df.index:
            peak_file = '{}/{}/{}'.format(cistrome_dir,batch,df.loc[cistrome_id,'File'])
            return peak_file


cistrome_dir = '/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019'
# cistrome_dir = '/Volumes/zanglab/zw5j/data/cistrome/cistrome_db_2019'
pardir = 'f3_union_peaks_across_cellTypes'
os.makedirs(pardir,exist_ok=True)


list_df = pd.read_csv('f2_selected_data/pairwise_sharedTF_count_id.csv',index_col=0)
for celltype_pairs in list_df.index[:8]:
    outdir = pardir+os.sep+celltype_pairs
    os.makedirs(outdir,exist_ok=True)
    
    outfile = '{}/run_{}.slurm'.format(pardir,celltype_pairs)
    outf = open(outfile,'w')
    outf.write('''#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=8000
#SBATCH -t 4:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
''')            
    outf.write('#SBATCH -o {}/slurm_{}.out\n\n'.format(pardir,celltype_pairs))
    
    data_df = pd.read_csv('f2_selected_data/selected_data_IDs_{}.csv'.format(celltype_pairs),index_col=0)   
    for factor in data_df.index[:]:
        outf.write('echo {} \n'.format(factor))
        factor_df = data_df.loc[factor].dropna()
        factor_ids = factor_df.values.astype(int)
        peak_files = [return_cistrome_peak_file(factor_id) for factor_id in factor_ids]
        # merge all data from the same factor
        commandline = 'cat {} > {}/{}.cat.bed\n'.format(' \\\n'.join(peak_files),outdir,factor)
        outf.write(commandline)
        commandline = 'bedtools sort -i {}/{}.cat.bed > {}/{}.sort.bed\n'.format(outdir,factor,outdir,factor)
        outf.write(commandline)
        commandline = 'bedtools merge -i {}/{}.sort.bed > {}/{}.merge.bed\n'.format(outdir,factor,outdir,factor)
        outf.write(commandline)
        # bedtools intersect
        for celltype in factor_df.index:
            peak_file = return_cistrome_peak_file(factor_df.loc[celltype])
            commandline = 'bedtools intersect -a {}/{}.merge.bed -b {} -c > {}/{}_merge_intersect_{}.bed\n'.format(outdir,factor,peak_file,outdir,factor,celltype)
            outf.write(commandline)
        for celltype in factor_df.index:
            peak_file = return_cistrome_peak_file(factor_df.loc[celltype])
            commandline = 'cp {} {}/_{}_{}.bed\n'.format(peak_file,outdir,factor,celltype)
            outf.write(commandline)
        # add new line to the end per factor    
        outf.write('\n')
    outf.close()
            
        
        
            





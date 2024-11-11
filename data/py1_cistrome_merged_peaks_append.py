import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect
import random
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


data_dir = '/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019'
# data_dir = '/Volumes/zanglab/zw5j/data/cistrome/cistrome_db_2019'
outdir = 'cistrome_data'
os.makedirs(outdir,exist_ok=True)
excluede = ['HBG1, HBG2','SMAD2/3','T','H3K9ac, H3K14ac']

out_df = pd.DataFrame()
for batch in ['human_factor','human_hm']:  
    batch_dir = '{}/{}'.format(data_dir,batch)
    df = pd.read_csv('{}/{}.txt'.format(data_dir,batch),sep='\t',index_col=0)
    factors = sorted(df.Factor.unique())
    factors = [i for i in factors if not i in excluede]
    print('\n====\n',batch,len(factors),'\n====\n',)

    for factor in ['H3K27ac','H3K27me3','H3K4me3']:
        factor_df = df[df.Factor==factor]
        factor_ids = factor_df.index.astype(int)
        if not len(factor_ids)>=5:
            continue
        out_df.loc[factor,'#datasets']=len(factor_ids)
        peak_files = ['{}/{}'.format(batch_dir,factor_df.loc[factor_id,'File']) for factor_id in factor_ids]
        peak_files = random.sample(peak_files,1000)
        print(factor,len(peak_files))
        # peak_files = ['{}/{}_sort_peaks.narrowPeak.bed'.format(batch_dir,factor_id) for factor_id in factor_ids]
        # for peak_file in peak_files:
        #     if not os.path.isfile(peak_file):
        #         print(peak_file)
        
        # merge all data from the same factor
        commandline = 'cat {} > {}/{}.cat.bed'.format(' \\\n'.join(peak_files),outdir,factor)
        os.system(commandline)
        print(commandline)
        commandline = 'bedtools sort -i {}/{}.cat.bed > {}/{}.sort.bed'.format(outdir,factor,outdir,factor)
        os.system(commandline)
        print(commandline)
        commandline = 'bedtools merge -i {}/{}.sort.bed > {}/{}.merge.bed'.format(outdir,factor,outdir,factor)
        os.system(commandline)
        print(commandline)
        print()

out_df.to_csv('cistrome_data_summary.csv')
    


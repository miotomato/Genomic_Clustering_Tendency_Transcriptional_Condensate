import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect
import glob
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


def write_slurm(slurm_dir,cistrome_id,factor,celltype):
    basename = '{}_{}_{}'.format(cistrome_id,factor,celltype)
    slurmfile = slurm_dir+os.sep+'run_{}.slurm'.format(basename)
    with open(slurmfile,'w') as slurmout:
        slurmout.write('''#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
''')

        #######################
        # REMEMBER TO CHANGE THE MEMORY !!!
        #########################
            
        slurmout.write('#SBATCH -o {}/slurm_{}.out\n\n'.format(slurm_dir,basename))
        slurmout.write('time python py1_TFBS_enrichAt_SE.py -i {} -f {} -c {}\n\n'.format(cistrome_id,factor,celltype))



slurm_dir = 'f1_slurm'
outdir = 'f1_slurm_log'
os.makedirs(slurm_dir,exist_ok=True)
os.makedirs(outdir,exist_ok=True)

# read SE cell types
SE_dir = '../../data/SE_hg38/'
SE_files = glob.glob('{}/*.bed'.format(SE_dir))
SE_celltypes = [os.path.basename(i).split('.bed')[0] for i in SE_files]

# read cistrome data
data_dir = '/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new'
# data_dir = '/Volumes/zanglab/zw5j/data/cistrome/cistrome_db_2019_new'
excluded_factors = ['HBG1, HBG2','SMAD2/3','T','H3K9ac, H3K14ac','H2AZK4ac, H2AZK7ac, H2AZK11ac']

out_df = pd.DataFrame()
out_df_GSM = pd.DataFrame()
for batch in ['human_factor','human_hm'][:]:  
    df = pd.read_csv('{}/{}_full_QC.txt'.format(data_dir,batch),sep='\t',index_col=0)
    df = df.dropna()
    # df = df[df.Cell_line.str.lower().isin([i.lower() for i in SE_celltypes])]
    df = df[df.Cell_line.isin(SE_celltypes)]
    df = df[~df.Factor.isin(excluded_factors)]
    df = df[df.PeaksFoldChangeAbove10>200]
    df.to_csv('{}/{}.csv'.format(outdir,batch))
    factors = sorted(df.Factor.unique())
    celltypes = sorted(df.Cell_line.unique())
    print('\n====\n',batch,len(factors),'\n====\n',)
    
    for factor in factors:
        for celltype in celltypes:
            # check if the existance of merge peak file
            if not os.path.isfile('../../data/cistrome_data/{}.merge.bed'.format(factor)):
                continue
            # for each factor/cell-type, only use the ChIP-seq data with most number of peaks
            df_sub = df[(df.Factor==factor)&(df.Cell_line==celltype)]
            # if celltype=='HeLa':
            #     df_sub = df[(df.Factor==factor)&
            #                 ((df.Cell_line=='HeLa')|(df.Cell_line=='HeLa-S3'))]
            if df_sub.shape[0]>0:
                cistrome_id = df_sub.sort_values(by='PeaksFoldChangeAbove10',ascending=False).index[0]
                write_slurm(slurm_dir,cistrome_id,factor,celltype)
                gsmID = df_sub.loc[cistrome_id].GSMID
                out_df.loc[factor,celltype] = cistrome_id
                out_df_GSM.loc[factor,celltype] = gsmID
    
# save the data
out_df.to_csv('{}/cistrome_id_log.csv'.format(outdir))
out_df_GSM.to_csv('{}/cistrome_GSM_log.csv'.format(outdir))
        
        
        
    


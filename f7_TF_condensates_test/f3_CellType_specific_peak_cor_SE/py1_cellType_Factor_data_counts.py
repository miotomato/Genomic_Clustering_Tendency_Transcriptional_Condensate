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



# data_dir = '/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new'
data_dir = '/Volumes/zanglab/zw5j/data/cistrome/cistrome_db_2019_new'
outdir = 'f1_cellType_Factor_data_Count'
os.makedirs(outdir,exist_ok=True)
    
df1 = pd.read_csv('{}/human_factor_full_QC.txt'.format(data_dir),sep='\t',index_col=0)
df2 = pd.read_csv('{}/human_hm_full_QC.txt'.format(data_dir),sep='\t',index_col=0)
df = pd.concat([df1,df2])
df = df.dropna()
df = df[df.PeaksFoldChangeAbove10>500]
df.to_csv('{}/cistrome2019_data_details.csv'.format(outdir))           

celltypes = sorted(df.Cell_line.unique())
factors = sorted(df.Factor.unique())
df_out = pd.DataFrame(index=factors, columns=celltypes).fillna(0)
for ii in df.index:
    factor = df.loc[ii,'Factor']
    celltype = df.loc[ii,'Cell_line']
    df_out.loc[factor,celltype]+=1

# == remove sparse elements
num_thre = 10
df_out = df_out.loc[df_out.astype(bool).sum(axis=1)>num_thre,
                    df_out.astype(bool).sum(axis=0)>num_thre]
# == sort columns by total number of data
df_out = df_out.loc[df_out.astype(bool).sum(axis=1).sort_values(ascending=False).index]
df_out = df_out[df_out.astype(bool).sum().sort_values(ascending=False).index]
df_out.to_csv(outdir+os.sep+'cistrome2019_data_Count.csv')
    
    


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

# == select factors with sufficient data
df_count['HeLa_all'] = df_count[['HeLa-S3', 'HeLa']].sum(axis=1)
selected_samples = ['K562', 'MCF-7', 'HepG2', 'A549','H1', 'HCT-116', 'HeLa_all']
selected_df = df_count[selected_samples]
selected_df = selected_df.astype(bool).sum(axis=1).sort_values(ascending=False)
selected_factors = selected_df[selected_df>=6].index
selected_factors = np.append(selected_factors,'MED1')

# == save the selected chip-seq data
selected_data_df = pd.DataFrame()
selected_data_id = pd.DataFrame()
selected_data_GSM = pd.DataFrame()
for factor in selected_factors[:]:
    for celltype in selected_samples[:]:
        df_sub = df[(df.Factor==factor)&(df.Cell_line==celltype)]
        if celltype=='HeLa_all':
            df_sub = df[(df.Factor==factor)&
                        ((df.Cell_line=='HeLa')|(df.Cell_line=='HeLa-S3'))]
        if df_sub.shape[0]>0:
            # kept_index = df_sub.dropna().sort_values(by='FRiP',ascending=False).index[0]
            kept_index = df_sub.sort_values(by='PeaksFoldChangeAbove10',ascending=False).index[0]
            selected_data_df = pd.concat([selected_data_df,df_sub.loc[[kept_index]]])
            celltype = celltype.split('_')[0]
            selected_data_id.loc[factor,celltype]=int(kept_index)
            gsmID = df_sub.loc[kept_index].GSMID
            selected_data_GSM.loc[factor,celltype]=gsmID

# == save the data
selected_data_df.to_csv('{}/selected_data_details.csv'.format(outdir))    
selected_data_id.to_csv('{}/selected_data_IDs.csv'.format(outdir))           
selected_data_GSM.to_csv('{}/selected_data_GSMID.csv'.format(outdir))           





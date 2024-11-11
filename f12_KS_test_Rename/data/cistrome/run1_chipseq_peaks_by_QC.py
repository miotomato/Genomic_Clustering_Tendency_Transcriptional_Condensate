import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect


# ==== keep high quality data
# data_dir = '/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new'
data_dir = '/Volumes/zanglab/zw5j/data/cistrome/cistrome_db_2019_new'
df_count = pd.read_csv('{}/human_factor_peak_num.txt'.format(data_dir),sep='\s+',header=None)
df_count.columns = ['total_peaks','file_name']
df_count = df_count.iloc[:df_count.shape[0]-1]
df_count['DCid'] = [int(os.path.basename(ii).split('_')[0]) for ii in df_count.file_name]
df_count.index = df_count.DCid
df_count = df_count['total_peaks']

df = pd.read_csv('{}/human_factor_full_QC.txt'.format(data_dir),sep='\t',index_col=0)
df = pd.concat([df,df_count],axis=1)

# remove potential un-accurate data
df = df[df.Factor!='T'].dropna()
df = df[(df.PeaksFoldChangeAbove10>200)&(df.FastQC>15)&
        (df.UniquelyMappedRatio>0.3)&(df.PBC>0.3)&
        (df.PeaksUnionDHSRatio>0.3)&(df.FRiP>0.005)&
        (df.total_peaks>2000)]
# df = df[(df.PeaksFoldChangeAbove10>200)&(df.total_peaks>2000)]

df.Factor = df.Factor.replace('NOTCH1','RBPJ')
df.Cell_line = df.Cell_line.replace('Jurkat E6-1','T-ALL')
df.Cell_line = df.Cell_line.replace('Jurkat','T-ALL')
df.Cell_line = df.Cell_line.replace('CUTLL1','T-ALL')
df.Cell_line = df.Cell_line.replace('HeLa-S3','HeLa')
df.to_csv('cistrome2019_selected_QC.csv')           


# ==== count the data per cell type per factor
celltypes = sorted(df.Cell_line.unique())
factors = sorted(df.Factor.unique())
df_out = pd.DataFrame(index=factors, columns=celltypes).fillna(0)
for ii in df.index:
    factor = df.loc[ii,'Factor']
    celltype = df.loc[ii,'Cell_line']
    df_out.loc[factor,celltype]+=1
# == repeatively remove sparse elements
num_thre = 1
df_out = df_out.loc[df_out.astype(bool).sum(axis=1)>=num_thre,
                    df_out.astype(bool).sum(axis=0)>=num_thre]
# == sort columns by total number of data
df_out = df_out.loc[df_out.astype(bool).sum(axis=1).sort_values(ascending=False).index]
df_out = df_out[df_out.astype(bool).sum().sort_values(ascending=False).index]
df_out.to_csv('cistrome2019_data_Count.csv')
print(df_out.shape,df_out.sum().sum())
 
# ==== check overlap with motif and SE
# motif_dir="/nv/vol190/zanglab/shared/Motif/sites/hg38_fimo_jarspar/results/"
motif_dir="/Volumes/zanglab/shared/Motif/sites/hg38_fimo_jarspar/results/"
motif_files = glob.glob('{}/*.bed'.format(motif_dir))
motifs = [os.path.basename(i).split("_")[0] for i in motif_files if not re.search('~',i)]
SE_files = glob.glob('../SE_hg38/*.bed')
SEs = [os.path.basename(i).split(".bed")[0] for i in SE_files ]

# ==== profiles with motif
df_out = df_out.loc[df_out.index.intersection(motifs)]
df_out = df_out.loc[df_out.astype(bool).sum(axis=1)>=num_thre,
                    df_out.astype(bool).sum(axis=0)>=num_thre]
df_out.to_csv('cistrome2019_data_Count_with_motif.csv')
print(df_out.shape,df_out.sum().sum())

# ==== profiles with motif and SEs
df_out = df_out.loc[df_out.index.intersection(motifs),
                    df_out.columns.intersection(SEs)]
df_out = df_out.loc[df_out.astype(bool).sum(axis=1)>=num_thre,
                    df_out.astype(bool).sum(axis=0)>=num_thre]
df_out.to_csv('cistrome2019_data_Count_with_SE_motif.csv')
print(df_out.shape,df_out.sum().sum())




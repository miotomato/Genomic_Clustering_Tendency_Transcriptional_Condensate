import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect


pardir='/nv/vol190/zanglab/zw5j'
# pardir='/Volumes/zanglab/zw5j'

tf_file = '{}/work2017/tf_prediction/tf_DHS_analysis/info_analysis/3490_analysis/tf_repeats.txt'.format(pardir)
udhs_file='{}/data/unionDHS/hg38_unionDHS_fc4_50merge.bed'.format(pardir)

tf_df = pd.read_csv(tf_file)
tf_col = tf_df.columns[0]
udhs_df = pd.read_csv(udhs_file,header=None,sep='\t')

outfile = open('tf_on_UDHS.bed','w')
for ii in udhs_df.index[:]:
    chrome = udhs_df.loc[ii,0]
    start = udhs_df.loc[ii,1]
    end = udhs_df.loc[ii,2] 
    tf_info = tf_df.loc[ii,tf_col].split('\t') 
    dnaseID  = tf_info[0]
    total_binding  = tf_info[1]
    tf_names  = ','.join(tf_info[2:][::3])
    outfile.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(chrome,start,end,dnaseID,total_binding,tf_names))
outfile.close()



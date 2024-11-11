import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *
from scipy import stats
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#matplotlib.rcParams['font.size']=16
#matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
#matplotlib.rcParams["font.family"] = "sans-serif"
#import seaborn as sns
#sns.set(font_scale=2)
#sns.set_style("whitegrid", {'axes.grid' : False})

import re,bisect
#plus = re.compile('\+')
#minus = re.compile('\-')



# ==== main

# ATAC identifier
atac_id_info = pd.read_csv('TCGA_identifier_mapping.txt',sep='\t',index_col=3)
caseIDs = sorted(set(atac_id_info.index))

sig_file = 'mynorm_TCGA-ATAC_PanCan_Log2_QuantileNorm_Counts_plus5.txt'
out_sig_file = 'mynorm_TCGA-ATAC_PanCan_Log2_QuantileNorm_Counts_plus5.caseID.avg.txt'
with open(sig_file) as sig_inf:
    sig_df = pd.read_csv(sig_inf,sep='\t')

# averaged info for each caseID 
caseID_sig_df = pd.DataFrame(index = sig_df.index) # each row is a binding position
for caseID in caseIDs[:]:
    caseID_reps = atac_id_info.loc[caseID][['bam_prefix']].values
    if len(caseID_reps.shape)>1:
        caseID_reps = caseID_reps[:,0]
    caseID_reps = ['_'.join(i.split('-')) for i in caseID_reps]
    # caseID with cancertype+patient-replicate-ID (e.g., ), and matched patient uniq ID (e.g., 00f0f7dd-71de-4e4f-b4b5-df860324f2e8)
    caseID_reps = sig_df.columns.intersection(caseID_reps)#;print(caseID,caseID_reps);exit()
    if len(caseID_reps)>0:
        cancertype_readed = caseID_reps[0].split('_')[0]
        # print(cancertype_readed)
        caseID_sig_df['{}_{}'.format(cancertype_readed,caseID)] = sig_df[caseID_reps].mean(axis=1)
    
for insert_col in ['#seqnames','start', 'end', 'name', 'score', 'annotation', 'GC'][::-1]:
    caseID_sig_df.insert(0,insert_col,sig_df[insert_col])
caseID_sig_df.to_csv(out_sig_file,sep='\t',index=False)
#exit()
        

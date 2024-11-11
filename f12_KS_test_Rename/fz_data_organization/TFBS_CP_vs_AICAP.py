import os,sys,argparse
# import fileinput,time
import glob
import re,bisect
import pandas as pd
import numpy as np
from scipy import stats
import glob


            
            
project_dir='/Volumes/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/f12_KS_test_Rename'
# project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/'

# ==== read TFMS CP
indir =  '{}//f1_TF_cluster_potential/f2_cor_CP_SE_AICAP/f8_per_CT_cor_AICAP_cor_CP_with_motif_only/_csv'.format(project_dir)
infiles = glob.glob('{}/*csv'.format(indir))
infiles.sort()
i=0
writer = pd.ExcelWriter('data/TFBS_CP_vs_AICAP.xlsx')
for infile in infiles:
    cellType = re.split('\.',os.path.basename(infile))[0]
    df = pd.read_csv(infile,index_col=0)
    df['-log2 AICAP'] = -1*df['log2 AICAP']
    df = df[['TFBS CP', '-log2 AICAP']]
    df = df.dropna()
    df = df[df['-log2 AICAP']>0]
    # print(cellType,df.shape)
    df = df.sort_values(by='-log2 AICAP',ascending=False)
    if df.shape[0]>=3 and cellType !='None':
        df.to_excel(writer,cellType)
        print(cellType,df.shape,i)
        i=i+1
        
writer.close()
            
    
    
    
    
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
# indir =  '{}/f1_TF_cluster_potential/f2_cor_CP_SE_AICAP/f9_per_CT_TFBS_CP_cor_zscore_CP_with_motif_SE/TFBS_CP/'.format(project_dir)
indir =  '{}/f1_TF_cluster_potential/f2_cor_CP_SE_AICAP/f5_per_CT_TFBS_CP_cor_SE/TFBS_CP/'.format(project_dir)
infiles = glob.glob('{}/*csv'.format(indir))
infiles.sort()

writer = pd.ExcelWriter('data/TFBS_CP_vs_CP_zscore.xlsx')
for infile in infiles:
    cellType = re.split('_|\.',os.path.basename(infile))[-2]
    print(cellType)
    df = pd.read_csv(infile,index_col=0)
    df = df[['TFBS CP', 'Z-scored CP', 'avg rank','log2 Odds Ratio']]
    df.columns = ['TFBS CP', 'CP Z-score', 'Average rank','SE enrichment log2 Odds Ratio']
    df.to_excel(writer,cellType)
    
writer.close()
            
    
    
    
    
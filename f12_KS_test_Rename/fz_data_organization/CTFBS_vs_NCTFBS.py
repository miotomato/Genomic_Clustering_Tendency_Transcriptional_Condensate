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
indir =  '{}//f2_TCGA_clinical/fz_combi_figs/f2_ATAC_HiC_compr_TFBS_heatmap'.format(project_dir)
infiles = glob.glob('{}/*ExtendMerge_50KB.csv'.format(indir))
infiles.sort()

writer = pd.ExcelWriter('data/CTFBS_vs_NCTFBS_chromatin_activity.xlsx')
for infile in infiles:
    cellType = re.split('_',os.path.basename(infile))[1]
    print(cellType)
    df = pd.read_csv(infile,index_col=0)
    df.columns = ['SE enrichment Odds Ratio', 'SE enrichment pvalue',
       'Differential ATAC-seq T-statistics', 'Differential ATAC-seq pvalue',
       'ATAC-seq RP T-statistics', 'ATAC-seq RP pvalue', 'Hi-C CI statistics',
       'Hi-C CI pvalue']
    df.to_excel(writer,cellType)        
writer.close()
            
    
    
    
    
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


outdir = 'f2_diffATAC_sig_by_clinical'
os.makedirs(outdir,exist_ok=True)
    
data_dir='/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF'
# data_dir='/Volumes/zanglab/zw5j/work2017/T_ALL_CTCF'
clinical_date = '2019-10-09'    
cancertypes = ['BRCA','COAD','LUAD','PRAD']

# ATAC identifier
atac_id_info = pd.read_csv('{}/15_TCGA_patient_data/ATAC_seq/data/TCGA_identifier_mapping.txt'.format(data_dir),sep='\t',index_col=3)
atac_id_info['bam_prefix_new'] = ['_'.join(i.split('-')) for i in atac_id_info.bam_prefix]
caseIDs = sorted(set(atac_id_info.index))

# get the differential ATAC signal
for cancertype in cancertypes[:3]:
    
    # read the clinical info
    clinical_df = pd.read_csv('f1_clinical_data_JSON_{}/{}_clinical_info.csv'.format(clinical_date,cancertype),index_col=0)
    
    # read the ids that have long days to last follow up
    alive_thre = np.percentile(clinical_df.days_to_last_follow_up.dropna(),q=50)
    alive_ids = clinical_df[clinical_df.days_to_last_follow_up>alive_thre].index
    alive_ids = atac_id_info.index.intersection(alive_ids)
    alive_bam_index = atac_id_info.loc[alive_ids].bam_prefix_new
    
    # read the ids that have short days to death
    dead_thre = np.percentile(clinical_df.days_to_death.dropna(),q=50)
    dead_ids = clinical_df[clinical_df.days_to_death < dead_thre].index
    dead_ids = atac_id_info.index.intersection(dead_ids)
    dead_bam_index = atac_id_info.loc[dead_ids].bam_prefix_new
    print('# ids of alive/dead',len(alive_ids),len(dead_ids))
    print('# bam of alive/dead',len(alive_bam_index),len(dead_bam_index))
    
    # read atac sig
    atac_file = '{}/15_TCGA_patient_data/ATAC_seq/data/Count_matrices/TCGA-ATAC_Cancer_Type-specific_Count_Matrices_log2norm_counts/{}_log2norm.txt'.format(data_dir,cancertype)
    atac_df = pd.read_csv(atac_file,sep='\t')
    s,p = stats.ttest_ind(atac_df[dead_bam_index],atac_df[alive_bam_index],axis=1)
    
    atac_sig = atac_df[['seqnames', 'start', 'end']]
    atac_sig['stats'] = s
    atac_sig.to_csv(outdir+os.sep+'{}_diffATAC_DoverA.bed'.format(cancertype),sep='\t',index=False,header=None)
    
    atac_sig['stats'] = -1*s
    atac_sig.to_csv(outdir+os.sep+'{}_diffATAC_AoverD.bed'.format(cancertype),sep='\t',index=False,header=None)
    

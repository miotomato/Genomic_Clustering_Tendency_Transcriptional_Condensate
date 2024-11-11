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

rivanna_dir = '/Volumes/zanglab/zw5j/'
rivanna_dir = '/nv/vol190/zanglab/zw5j/'
bigwig_dir = '{}/work2017/T_ALL_CTCF/15_TCGA_patient_data/ATAC_seq/data/bigWigs/'.format(rivanna_dir)

atac_id_info = pd.read_csv('../TCGA_identifier_mapping.txt',sep='\t',index_col=3)
# caseIDs = sorted(set(atac_id_info.index))
filtered_df = pd.read_excel('../TCGA-ATAC_clustered_samples.xlsx',index_col=1)   
name_match = pd.read_excel('..//TCGA-ATAC_SE_cancerType_match.xlsx',index_col=0)   
name_match = name_match.dropna()



for cancertype in ['BRCA','COAD','CESC','LIHC','PRAD'][:]:        
    
    outdir = cancertype
    os.makedirs(outdir,exist_ok=True)
    
    # clinical_df = pd.read_csv('..//clinical.project-TCGA-{}.2022-03-20.csv'.format(cancertype),index_col=1)
    filtered_ids = filtered_df[filtered_df.cohort==cancertype].index
    for filtered_id in filtered_ids[:]:
        # == read the submitter_id and bam prefix 
        caseID_reps = atac_id_info.loc[filtered_id][['bam_prefix']].values
        submitter_id = filtered_df.loc[filtered_id].submitter_id
        # clinical_time = clinical_df.loc[submitter_id,'time'].astype(int)
        # clinical_status = clinical_df.loc[submitter_id,'status'].astype(int)
        # == in case there are multiple replicates
        if len(caseID_reps.shape)>1:
            caseID_reps = caseID_reps[:,0]
        # == copy each replicates
        for ii in np.arange(caseID_reps.shape[0]):
            caseID_rep = caseID_reps[ii]
            # flag = caseID_rep.split('-')[-2]
            caseID_rep = '_'.join(caseID_rep.split('-'))
            bigwig_file = '{}/{}_bigWigs/{}.insertions.bw'.format(bigwig_dir,cancertype,caseID_rep)
            # cp_file = '{}_status{}_time{}_{}_rep{}.insertions.bw'.format(cancertype,clinical_status,clinical_time,submitter_id,ii+1)
            cp_file = '{}_{}_rep{}.insertions.bw'.format(cancertype,submitter_id,ii+1)
            commandLine = 'cp {} {}/{}'.format(bigwig_file,outdir, cp_file)
            os.system(commandLine)
            print(commandLine)
        

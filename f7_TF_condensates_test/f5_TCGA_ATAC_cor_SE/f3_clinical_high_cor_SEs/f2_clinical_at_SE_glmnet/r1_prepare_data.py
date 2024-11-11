import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import json
#from GenomeData import *
import re,bisect
 


def return_survival_df(df):

    df['time_max'] = df[['days_to_death','days_to_last_follow_up']].fillna(0).astype(float).max(axis=1)   
    df['time_sum'] = df[['days_to_death','days_to_last_follow_up']].fillna(0).astype(float).sum(axis=1)   
#     df = df[df['time_max']==df['time_sum']]
    df['time']=df['time_max']
    df = df[df['time']>0]
    
    df['status']= df.loc[df.index]['vital_status'] 
    df.loc[df['status']=='Dead','status']=1
    df.loc[df['status']=='Alive','status']=0
    return df

    
# ==== main  
outdir = 'f1_data'
os.makedirs(outdir,exist_ok=True)
name_match = pd.read_excel('../../f1_TCGA_ATAC_at_SE_rep/TCGA_ATAC_SE_cancerType_match.xlsx',index_col=0)   
name_match = name_match.dropna()
se_dir = '../../f2_TCGA_ATAC_at_SE_cor_survival/f1_ATAC_overlap_SE_caseID'
clinical_dir = '../f1_clinical_at_each_SE/f2_caseID_each_SE_vs_clinical/'
 

# for cancertype in name_match.index[:]:
for cancertype in ['BRCA','COAD']:
    cancertype_SE = name_match.loc[cancertype].SE
    cancertype_SE_rename = name_match.loc[cancertype].SE_rename
    
    # ==== read TCGA peaks overlap with SE data
    overlap_file = '{}/{}_ATAC_overlap_{}_SE_caseID.bed'.format(se_dir,cancertype,cancertype_SE_rename)
    with open(overlap_file) as sig_inf:
        sig_df = pd.read_csv(sig_inf,sep='\t',index_col=3)
    sig_df = sig_df.iloc[:,6:]
    # patient ID by TCGA data
    sig_case_id = [i.split('_')[1] for i in sig_df.columns]
    sig_df.columns = sig_case_id
    sig_df = np.transpose(sig_df)
    
    
    # ==== read clinical data
    clinical_df = pd.read_csv('{}/{}_clinical_info.csv'.format(clinical_dir,cancertype),index_col=0)
    clinical_df =  return_survival_df(clinical_df)
    clinical_df = clinical_df[['time','status']]

    # ==== prepare data for Glmnet
    shared_case = clinical_df.index.intersection(sig_df.index)
    sig_df.loc[shared_case].to_csv('{}/{}_atacseq_sig.csv'.format(outdir,cancertype))
    clinical_df.loc[shared_case].to_csv('{}/{}_clinical.csv'.format(outdir,cancertype))
         

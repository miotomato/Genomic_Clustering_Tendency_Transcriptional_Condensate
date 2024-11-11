import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import json
#from GenomeData import *
import re,bisect




def return_survival_df(region_df,region,cut_off=0.5):

    df = pd.DataFrame(index = region_df.index)
    df['time_max'] = region_df[['days_to_death','days_to_last_follow_up']].fillna(0).astype(float).max(axis=1)   
    df['time_sum'] = region_df[['days_to_death','days_to_last_follow_up']].fillna(0).astype(float).sum(axis=1)   
#     df = df[df['time_max']==df['time_sum']]
    df['time']=df['time_max']
    df = df[df['time']>0]
    
    # ==== get the clinical info, using top/bottom 50/25% as cutoff
    # ==== keep at least 5 patients
    a_index = df.index[:max(int(len(df.index)*cut_off),0)]
    b_index = df.index[-1*max(int(len(df.index)*cut_off),0):]
    df = df.loc[a_index.union(b_index)]
    # print('treat:\t',len(a_index),'\nctrl:\t',len(b_index))
    df.loc[a_index,'group']='treat'   
    df.loc[b_index,'group']='ctrl'  
    df['status']= region_df.loc[df.index]['vital_status'] 
    df['score']= region_df.loc[df.index][region] 
    df.loc[df['status']=='Dead','status']=1
    df.loc[df['status']=='Alive','status']=0
    return df

    



    
# ==== main  
# indir = 'f1_ATAC_overlap_SE_caseID'
outdir = 'f1_SE_openness_add_clinical'
os.makedirs(outdir,exist_ok=True)
# os.makedirs(outdir+os.sep+'fig',exist_ok=True)

name_match = pd.read_excel('../../data/TCGA/TCGA-ATAC_SE_cancerType_match.xlsx',index_col=0)   
name_match = name_match.dropna()
atac_seq_dir = '../../f3_clinical_hicor_SE_overlapped_peaks/f1_clinical_at_each_peak/f1_ATAC_overlap_SE_caseID'  
clinical_dir = '../../f3_clinical_hicor_SE_overlapped_peaks/f1_clinical_at_each_peak/f2_caseID_each_SE_vs_clinical' 


for cancertype in ['BRCA','COAD'][:]:
    cancertype_SE = name_match.loc[cancertype].SE
    cancertype_SE_rename = name_match.loc[cancertype].SE_rename
    atac_df = pd.read_csv('{}/{}_ATAC_overlap_{}_SE_caseID.bed'.format(atac_seq_dir,cancertype,cancertype_SE_rename),sep='\t')
    # ==== read clinical data
    clinical_df = pd.read_csv('{}/{}_clinical_info.csv'.format(clinical_dir,cancertype),index_col=0)
    clinical_df['time'] = clinical_df[['days_to_death','days_to_last_follow_up']].fillna(0).astype(float).max(axis=1)   
    clinical_df = clinical_df[clinical_df['time']>0]
    clinical_df['status']= clinical_df['vital_status'] 
    clinical_df.loc[clinical_df['status']=='Dead','status']=1
    clinical_df.loc[clinical_df['status']=='Alive','status']=0

    # == read the SE info
    se_df = pd.read_csv('../data/{}_SE_with_overlapped_peaks.bed'.format(cancertype),sep='\t',header=None)    
    se_bed_columns = ['SE_chr','SE_start','SE_end']
    se_df.columns = np.append(se_bed_columns,atac_df.columns)
    se_df.index = ['_'.join(i.astype(str)) for i in se_df[se_bed_columns].values]
    
    # average ATAC-seq signal at each SE and combine with clnical info
    out_df = pd.DataFrame(index = clinical_df.index)
    for ii in se_df.index.unique()[:]:
        if len(se_df.loc[ii].shape)==1 :
            continue
        se_df_tmp = se_df.loc[ii].iloc[:,10:]
        se_df_tmp.columns = [i.split('_')[1] for i in se_df_tmp.columns]
        se_df_tmp = se_df_tmp[clinical_df.index].mean()
        out_df[ii] = se_df_tmp
    out_df['all_SE'] = out_df.mean(axis=1)
    out_df = pd.concat([clinical_df,out_df],axis=1)
    out_df.to_csv(outdir+os.sep+'{}_SE_openness_add_clinical.csv'.format(cancertype))
        
         

import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import json
#from GenomeData import *

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
# from stats_pvalues import irwin_hall_cdf



def irwin_hall_return_pvalue(df):
    # for each sample/col, get the rank (smaller is better) of each item (row)
    # calculate the rank sum of all items across all samples
    # df_rank_sum = df.rank(ascending=False).sum(axis=1)/df.shape[0]
    df_rank_sum = df.rank(ascending=False).mean(axis=1)
    df['rank_sum'] = df_rank_sum
#     df['irwin_hall_pvalue']=[irwin_hall_cdf(i,df.shape[1]) for i in df_rank_sum]
    return df
    
    
    
    
# main  
indir = 'f1_ATAC_overlap_SE_caseID'
outdir = 'f2_caseID_ranksum_vs_clinical'
os.makedirs(outdir,exist_ok=True)
name_match = pd.read_excel('../f1_TCGA_ATAC_at_SE_rep/TCGA_ATAC_SE_cancerType_match.xlsx',index_col=0)   
name_match = name_match.dropna()
  

for cancertype in name_match.index[:]:
    cancertype_SE = name_match.loc[cancertype].SE
    cancertype_SE_rename = name_match.loc[cancertype].SE_rename
    se_file = '../data/SE_hg38/{}.bed'.format(cancertype_SE)
    # ==== get TCGA peaks overlap with SE data
    overlap_file = '{}/{}_ATAC_overlap_{}_SE_caseID.bed'.format(indir,cancertype,cancertype_SE_rename)
    with open(overlap_file) as sig_inf:
        sig_df = pd.read_csv(sig_inf,sep='\t')
    
    # ==== rank caseID by TCGA ATAC-seq signal  
    sig_df = sig_df.iloc[:,7:]
    sig_df = np.transpose(sig_df)#;print(sig_df)
    sig_df = irwin_hall_return_pvalue(sig_df) # rank patient samples for each region/col
    sig_df = sig_df.sort_values(by=['rank_sum'],ascending=True)
    sig_df.to_csv(outdir+os.sep+'{}_irwin_hall.csv'.format(cancertype))
    
    # cancer_type_index = [i for i in sig_df.index if re.search(cancertype_match_names[name_type],i)]
    # sig_df_p = irwin_hall_return_pvalue(sig_df.loc[cancer_type_index])
    # sig_df_p = sig_df_p.sort_values(by=['rank_sum'],ascending=True)
    # sig_df_p.to_csv(outdir+os.sep+'{}_{}_irwin_hall_cancertype_ID.csv'.format(name_type,binding_type))

    # ==== read clinical data
    clinical_data = '2022-01-05'  
    clinical_file = '..//data_TCGA/clinical.cases_selection.{}.{}.json'.format(cancertype,clinical_data)
    with open(clinical_file) as clinical_inf:
        clinical_list = json.load(clinical_inf)
    clinical_case_id=[i['case_id'] for i in clinical_list]

    # ==== combi ranksum and clinical info
    sig_df = sig_df.iloc[:,-1:] # keep the rank sum col
    sig_case_id = [i.split('_')[1] for i in sig_df.index] # remove the cancer type prefix
        
    # search by clinical_list rather than sig_df
    # combine clinical info with openness info
    out_df = pd.DataFrame()
    for clinical_case in clinical_list:
        case_id = clinical_case['case_id']
        if case_id in sig_case_id:
            out_df.loc[case_id,'rank_sum'] = sig_df.loc['{}_{}'.format(cancertype,case_id)]['rank_sum']
            # out_df.loc[case_id,'tumor_stage'] = clinical_case['diagnoses'][0]['tumor_stage']
            out_df.loc[case_id,'age_at_diagnosis'] = clinical_case['diagnoses'][0]['age_at_diagnosis']
            out_df.loc[case_id,'vital_status'] = clinical_case['demographic']['vital_status'] 
            if out_df.loc[case_id,'vital_status']=='Alive':
                out_df.loc[case_id,'days_to_last_follow_up'] = clinical_case['diagnoses'][0]['days_to_last_follow_up']
            elif out_df.loc[case_id,'vital_status']=='Dead':
                out_df.loc[case_id,'days_to_death'] = clinical_case['demographic']['days_to_death']

    out_df = out_df.sort_values(by=['rank_sum'],ascending=True)
    out_df.index = ['{}_{}'.format(cancertype,i) for i in out_df.index]
    out_df.to_csv(outdir+os.sep+'{}_irwin_hall_with_clinical_info.csv'.format(cancertype))
         

import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import json
#from GenomeData import *
import re,bisect
# from lifelines.statistics import logrank_test
# from lifelines import KaplanMeierFitter



# ==== main  

for cancertype in ['BRCA','COAD'][:]:    
    # ==== read clinical data
    clinical_data = '2022-03-20'  
    clinical_file = 'clinical.project-TCGA-{}.{}.json'.format(cancertype,clinical_data)
    outfile = 'clinical.project-TCGA-{}.{}.csv'.format(cancertype,clinical_data)
    
    with open(clinical_file) as clinical_inf:
        clinical_list = json.load(clinical_inf)
    # patient ID by clinical data
    clinical_case_id=[i['case_id'] for i in clinical_list]

    
    # ==== read clinical info for each case ID
    clinical_df = pd.DataFrame()
    for clinical_case in clinical_list:
        case_id = clinical_case['case_id']

        try:
            clinical_df.loc[case_id,'submitter_id'] = clinical_case['demographic']['submitter_id'].split('_demographic')[0]
            clinical_df.loc[case_id,'vital_status'] = clinical_case['demographic']['vital_status'] 
            if clinical_df.loc[case_id,'vital_status']=='Alive':
                clinical_df.loc[case_id,'days_to_last_follow_up'] = clinical_case['diagnoses'][0]['days_to_last_follow_up']
            elif clinical_df.loc[case_id,'vital_status']=='Dead':
                clinical_df.loc[case_id,'days_to_death'] = clinical_case['demographic']['days_to_death']
        except:
            print('\n\n=======\n',case_id)


    clinical_df['time'] = clinical_df[['days_to_death','days_to_last_follow_up']].fillna(0).astype(float).max(axis=1)   
    clinical_df.loc[clinical_df['vital_status']=='Dead','status']=1
    clinical_df.loc[clinical_df['vital_status']=='Alive','status']=0

    clinical_df.index.name='case_id'
    clinical_df.to_csv(outfile)

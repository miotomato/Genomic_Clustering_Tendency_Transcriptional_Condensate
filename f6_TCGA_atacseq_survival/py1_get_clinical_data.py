import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *

import re,bisect
#plus = re.compile('\+')
#minus = re.compile('\-')
import json
    

clinical_date = '2019-10-09'    
outdir = 'f1_clinical_data_JSON_{}'.format(clinical_date)
os.makedirs(outdir,exist_ok=True)
data_dir='/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF'
data_dir='/Volumes/zanglab/zw5j/work2017/T_ALL_CTCF'

# selected cancer-specific ATAC info (overlapped with cancer specific CTCFs)
atac_sig_file_dir = 'f2_ATAC_sample_irwin_hall'
cancertypes = ['BRCA','COAD','LUAD','PRAD']

for cancer_type in cancertypes[:]:
    # this is for clinical info for each cancer type
    clinical_file = '{}/15_TCGA_patient_data/clinical_survival/data_TCGA_clinical/clinical.project-TCGA-{}.{}.json'.format(data_dir,cancer_type,clinical_date)
    with open(clinical_file) as clinical_inf:
        clinical_list = json.load(clinical_inf)
    # clinical_case_id=[i['case_id'] for i in clinical_list]
    
    # read the clinical info
    out_df = pd.DataFrame()
    for clinical_case in clinical_list[:]:
        case_id = clinical_case['case_id']
        # if 'demographic' in clinical_case.keys() and 'diagnoses' in clinical_case.keys():
        try:
            # out_df.loc[case_id,'tumor_stage'] = clinical_case['diagnoses'][0]['tumor_stage']
            # out_df.loc[case_id,'age_at_diagnosis'] = clinical_case['diagnoses'][0]['age_at_diagnosis']
            out_df.loc[case_id,'vital_status'] = clinical_case['demographic']['vital_status'] 
            if out_df.loc[case_id,'vital_status']=='Alive':
                out_df.loc[case_id,'days_to_last_follow_up'] = clinical_case['diagnoses'][0]['days_to_last_follow_up']
            elif out_df.loc[case_id,'vital_status']=='Dead':
                out_df.loc[case_id,'days_to_death'] = clinical_case['demographic']['days_to_death']
        except:
            pass
    
    out_df.to_csv(outdir+os.sep+'{}_clinical_info.csv'.format(cancer_type))

        

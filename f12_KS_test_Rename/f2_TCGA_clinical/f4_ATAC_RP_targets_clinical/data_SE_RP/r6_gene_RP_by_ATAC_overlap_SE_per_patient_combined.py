import os,sys,argparse,glob,re
import numpy as np
import pandas as pd
# import find_overlap_keep_info_NOT_sep_strand_asimport



# == main 
indir = 'ATAC_overlap_SE_per_patient_sig_RP'
outdir = 'ATAC_overlap_SE_per_patient_sig_RP_combined'
os.makedirs(outdir,exist_ok=True)


project_dir='/Volumes/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/'
# project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/'

filtered_file = '{}/f7_TF_condensates_test/f6_revised_TCGA_ATAC_cor_SE//data/TCGA/TCGA-ATAC_clustered_samples.xlsx'.format(project_dir)
filtered_df = pd.read_excel(filtered_file,index_col=0)   

# read matched names between TCGA ATAC and SE data
name_match = pd.read_excel('{}/f8_TF_condensates_V2/f3_clinical_outcome/data/TCGA-ATAC_SE_cancerType_match.xlsx'.format(project_dir),index_col=0)   
name_match = name_match.dropna()

for cancertype in ['BRCA','COAD'][:]:
    cancertype_SE = name_match.loc[cancertype].SE
    cancertype_SE_rename = name_match.loc[cancertype].SE_rename
    filtered_id = filtered_df[filtered_df.cohort==cancertype].case_id
    
    out_df = pd.DataFrame()
    for patient in filtered_id[:]:
        infile = '{}/{}_ATAC_overlap_{}_SE_{}_RP.tsv'.format(indir,cancertype,cancertype_SE,patient)
        df = pd.read_csv(infile,header=None,index_col=0,sep='\t')
        out_df = pd.concat([out_df,df],axis=1)
    
    out_df.columns = filtered_id
    out_df = out_df.loc[out_df.mean(axis=1)!=0]
    out_df.to_csv('{}/{}_ATAC_overlap_{}_SE_RP_combined.tsv'.format(outdir,cancertype,cancertype_SE))










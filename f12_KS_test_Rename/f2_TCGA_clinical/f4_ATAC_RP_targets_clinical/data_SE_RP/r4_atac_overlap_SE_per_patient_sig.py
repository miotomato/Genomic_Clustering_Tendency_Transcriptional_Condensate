import os,sys,argparse,glob,re
import numpy as np
import pandas as pd
# import find_overlap_keep_info_NOT_sep_strand_asimport



# == main 
# indir = 'ATAC_overlap_SE_overlap_TFMS'
outdir = 'ATAC_overlap_SE_per_patient_sig'
os.makedirs(outdir,exist_ok=True)


project_dir='/Volumes/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/'
# project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/f7_TF_condensates_test'
atac_overlap_SE_dir = '{}/f7_TF_condensates_test/f6_revised_TCGA_ATAC_cor_SE/f3_clinical_hicor_SE_overlapped_peaks/f1_clinical_at_each_peak/f1_ATAC_overlap_SE_caseID'.format(project_dir)
# atac_file = '{}/f6_revised_TCGA_ATAC_cor_SE/data/TCGA/tcga_atac.bed'.format(project_dir)

filtered_file = '{}/f7_TF_condensates_test/f6_revised_TCGA_ATAC_cor_SE//data/TCGA/TCGA-ATAC_clustered_samples.xlsx'.format(project_dir)
filtered_df = pd.read_excel(filtered_file,index_col=0)   

# read matched names between TCGA ATAC and SE data
name_match = pd.read_excel('{}/f8_TF_condensates_V2/f3_clinical_outcome/data/TCGA-ATAC_SE_cancerType_match.xlsx'.format(project_dir),index_col=0)   
name_match = name_match.dropna()
for cancertype in ['BRCA','COAD'][:]:
    cancertype_SE = name_match.loc[cancertype].SE
    cancertype_SE_rename = name_match.loc[cancertype].SE_rename
    atac_overlap_SE_file='{}/{}_ATAC_overlap_{}_SE_caseID.bed'.format(atac_overlap_SE_dir,cancertype,cancertype_SE)
    atac_overlap_SE_df = pd.read_csv(atac_overlap_SE_file,sep='\t',index_col=3)
    filtered_id = filtered_df[filtered_df.cohort==cancertype].case_id
    
    patient_col = atac_overlap_SE_df.columns[6:]
    patient_col = [i.split('_')[1] for i in patient_col]
    atac_overlap_SE_df.columns = np.append(atac_overlap_SE_df.columns[:6],patient_col)
                
    for patient in filtered_id[:]:
        kept_columns = ['#seqnames', 'start', 'end', patient, 'annotation']
        bed_df = atac_overlap_SE_df[kept_columns]
        bed_df.insert(3,'name',bed_df.index)
        bed_df.to_csv('{}/{}_ATAC_overlap_{}_SE_{}.bed'.format(outdir,cancertype,cancertype_SE,patient),sep='\t',index=False,header=None)





  


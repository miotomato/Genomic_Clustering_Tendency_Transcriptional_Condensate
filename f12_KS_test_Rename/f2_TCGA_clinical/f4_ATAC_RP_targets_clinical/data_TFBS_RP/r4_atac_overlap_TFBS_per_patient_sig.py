import os,sys,argparse,glob,re
import numpy as np
import pandas as pd
# import find_overlap_keep_info_NOT_sep_strand_asimport



# == main 
indir = 'f1_ATAC_overlap_TFBS_caseID'
outdir = 'ATAC_overlap_TFBS_per_patient_sig'
os.makedirs(outdir,exist_ok=True)


project_dir='/Volumes/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/'
project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/'
# atac_overlap_SE_dir = '{}/f7_TF_condensates_test/f6_revised_TCGA_ATAC_cor_SE/f3_clinical_hicor_SE_overlapped_peaks/f1_clinical_at_each_peak/f1_ATAC_overlap_SE_caseID'.format(project_dir)
# atac_file = '{}/f6_revised_TCGA_ATAC_cor_SE/data/TCGA/tcga_atac.bed'.format(project_dir)

# get patient samples
filtered_file = '{}/f7_TF_condensates_test/f6_revised_TCGA_ATAC_cor_SE//data/TCGA/TCGA-ATAC_clustered_samples.xlsx'.format(project_dir)
filtered_df = pd.read_excel(filtered_file,index_col=0)   

# read matched names between TCGA ATAC and SE data
name_match = pd.read_excel('{}/f8_TF_condensates_V2/f3_clinical_outcome/data/TCGA-ATAC_SE_cancerType_match.xlsx'.format(project_dir),index_col=0)   
name_match = name_match.dropna()
tcga_file = '{}/f7_TF_condensates_test/f6_revised_TCGA_ATAC_cor_SE//data/TCGA/mynorm_TCGA-ATAC_PanCan_Log2_QuantileNorm_Counts_plus5.caseID.avg.txt'.format(project_dir)

# get top3 factors
# c_tfbs_dir = '{}/f12_KS_test_Rename/f1_TF_cluster_potential/f3_clustered_TFBS/f5_atac_overlap_coBinding_TFBS'.format(project_dir)
selected_factors = {}
tfbs_cp_dir = '{}/f12_KS_test_Rename/f1_TF_cluster_potential/f2_cor_CP_SE_AICAP/f9_per_CT_TFBS_CP_cor_zscore_CP_with_motif_SE/TFBS_CP/'.format(project_dir)
for ct in ['MCF-7','HCT-116'][:]:
    df = pd.read_csv('{}/_CP_TFBS_all_vs_TFMS_{}.csv'.format(tfbs_cp_dir,ct),index_col=0)
    selected_factors['{} top_TFBSCP'.format(ct)] = df['TFBS CP rank'].sort_values().iloc[:3].index
    selected_factors['{} top_zscored_TFBSCP'.format(ct)] = df['avg rank'].sort_values().iloc[:3].index

# get atac-seq sig per patient
for cancertype in ['BRCA','COAD'][:]:
    ct = name_match.loc[cancertype].SE
    for treat_flag in ['percentile_T','percentile_T_ExtendMerge'][:]:
        for factorType in ['top_TFBSCP','top_zscored_TFBSCP'][:]:
            subdir = '{}_{}'.format(ct,factorType)
            os.makedirs(outdir+os.sep+subdir,exist_ok=True)            
            factors = selected_factors['{} {}'.format(ct,factorType)]            

            xticklabels = ['ALL','Union',
                           '{}'.format('-'.join([i for i in factors[:1]])),
                           '{}'.format('-'.join([i for i in factors[1:2]])),
                           '{}'.format('-'.join([i for i in factors[2:3]])),
                           '{}'.format('-'.join([i for i in factors[:2]])),
                           '{}'.format('-'.join([i for i in factors[:3]])),
                           ]

            for ii in np.arange(len(xticklabels))[:]:
                print(ct,treat_flag,xticklabels[ii])
                # tfbs_file = '{}/{}/{}_{}_{}.bed'.format(c_tfbs_dir,subdir,treat_flag,ct,xticklabels[ii])
                atac_overlap_file = '{}/{}/ATAC_overlap_{}_{}_{}.bed'.format(indir,subdir,treat_flag,ct,xticklabels[ii])
                if xticklabels[ii]=='ALL':
                     atac_overlap_file = tcga_file
                    
                # atac_overlap_file='{}/{}_ATAC_overlap_{}_SE_caseID.bed'.format(atac_overlap_SE_dir,cancertype,cancertype_SE)
                atac_overlap_df = pd.read_csv(atac_overlap_file,sep='\t',index_col=3)
                filtered_id = filtered_df[filtered_df.cohort==cancertype].case_id
                
                patient_col = atac_overlap_df.columns[6:]
                patient_col = [i.split('_')[1] for i in patient_col]
                atac_overlap_df.columns = np.append(atac_overlap_df.columns[:6],patient_col)
                            
                for patient in filtered_id[:]:
                    kept_columns = ['#seqnames', 'start', 'end', patient, 'annotation']
                    bed_df = atac_overlap_df[kept_columns]
                    bed_df.insert(3,'name',bed_df.index)
                    bed_df.to_csv('{}/{}/ATAC_overlap_{}_{}_{}_{}.bed'.format(outdir,subdir,treat_flag,ct,xticklabels[ii],patient),sep='\t',index=False,header=None)




  


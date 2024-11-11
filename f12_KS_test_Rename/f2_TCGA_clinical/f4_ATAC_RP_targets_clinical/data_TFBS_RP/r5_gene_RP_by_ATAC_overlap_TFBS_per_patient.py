import os,sys,argparse,glob,re
import numpy as np
import pandas as pd
# import find_overlap_keep_info_NOT_sep_strand_asimport



# == main 
indir = 'ATAC_overlap_TFBS_per_patient_sig'
outdir = 'ATAC_overlap_TFBS_per_patient_sig_RP'
os.makedirs(outdir,exist_ok=True)


project_dir='/Volumes/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/'
project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/'

filtered_file = '{}/f7_TF_condensates_test/f6_revised_TCGA_ATAC_cor_SE//data/TCGA/TCGA-ATAC_clustered_samples.xlsx'.format(project_dir)
filtered_df = pd.read_excel(filtered_file,index_col=0)   

# read matched names between TCGA ATAC and SE data
name_match = pd.read_excel('{}/f8_TF_condensates_V2/f3_clinical_outcome/data/TCGA-ATAC_SE_cancerType_match.xlsx'.format(project_dir),index_col=0)   
name_match = name_match.dropna()

# get top3 factors
selected_factors = {}
tfbs_cp_dir = '{}/f12_KS_test_Rename/f1_TF_cluster_potential/f2_cor_CP_SE_AICAP/f9_per_CT_TFBS_CP_cor_zscore_CP_with_motif_SE/TFBS_CP/'.format(project_dir)
for ct in ['MCF-7','HCT-116'][:]:
    df = pd.read_csv('{}/_CP_TFBS_all_vs_TFMS_{}.csv'.format(tfbs_cp_dir,ct),index_col=0)
    selected_factors['{} top_TFBSCP'.format(ct)] = df['TFBS CP rank'].sort_values().iloc[:3].index
    selected_factors['{} top_zscored_TFBSCP'.format(ct)] = df['avg rank'].sort_values().iloc[:3].index


for cancertype in ['BRCA','COAD'][:]:
    ct = name_match.loc[cancertype].SE
    filtered_id = filtered_df[filtered_df.cohort==cancertype].case_id

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
                # print(ct,treat_flag,xticklabels[ii])
                for patient in filtered_id[:]:
                    bed_file = '{}/{}/ATAC_overlap_{}_{}_{}_{}.bed'.format(indir,subdir,treat_flag,ct,xticklabels[ii],patient)
                    out_file = '{}/{}/ATAC_overlap_{}_{}_{}_{}_RP.tsv'.format(outdir,subdir,treat_flag,ct,xticklabels[ii],patient)
            
                    if os.path.isfile(bed_file):
                        commandLine = 'python2 get-regulatory-potential-on-genes_peak_level.py \
                            -b {} -s hg38 -g hg38_unique_geneSymbol.ucsc -o {}\n'.format(bed_file,out_file,)
                        os.system(commandLine)
                        print(commandLine)










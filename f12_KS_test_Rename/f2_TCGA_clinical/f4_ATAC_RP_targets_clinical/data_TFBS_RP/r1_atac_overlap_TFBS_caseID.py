import os,sys,argparse,glob
import numpy as np
import pandas as pd
# import find_overlap_keep_info_NOT_sep_strand_asimport




# == main 
outdir = 'f1_ATAC_overlap_TFBS_caseID'
os.makedirs(outdir,exist_ok=True)

# read matched names between TCGA ATAC and SE data
project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR'
# project_dir='/Volumes/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR'

name_match = pd.read_excel('{}/f7_TF_condensates_test/f6_revised_TCGA_ATAC_cor_SE//data/TCGA/TCGA-ATAC_SE_cancerType_match.xlsx'.format(project_dir),index_col=0)   
name_match = name_match.dropna()
pyfile="find_overlap_keep_info_NOT_sep_strand_revised.py"
tcga_file = '{}/f7_TF_condensates_test/f6_revised_TCGA_ATAC_cor_SE//data/TCGA/mynorm_TCGA-ATAC_PanCan_Log2_QuantileNorm_Counts_plus5.caseID.avg.txt'.format(project_dir)
c_tfbs_dir = '{}/f12_KS_test_Rename/f1_TF_cluster_potential/f3_clustered_TFBS/f5_atac_overlap_coBinding_TFBS'.format(project_dir)

# get top3 factors
selected_factors = {}
tfbs_cp_dir = '{}/f12_KS_test_Rename/f1_TF_cluster_potential/f2_cor_CP_SE_AICAP/f9_per_CT_TFBS_CP_cor_zscore_CP_with_motif_SE/TFBS_CP/'.format(project_dir)
for ct in ['MCF-7','HCT-116'][:]:
    df = pd.read_csv('{}/_CP_TFBS_all_vs_TFMS_{}.csv'.format(tfbs_cp_dir,ct),index_col=0)
    selected_factors['{} top_TFBSCP'.format(ct)] = df['TFBS CP rank'].sort_values().iloc[:3].index
    selected_factors['{} top_zscored_TFBSCP'.format(ct)] = df['avg rank'].sort_values().iloc[:3].index


# get overlap per C-TFBS 
for cancertype in ['BRCA','COAD'][:]:
    ct = name_match.loc[cancertype].SE
    for treat_flag in ['percentile_T','percentile_T_ExtendMerge'][:]:
        for factorType in ['top_TFBSCP','top_zscored_TFBSCP'][:]:
            subdir = '{}_{}'.format(ct,factorType)
            os.makedirs(outdir+os.sep+subdir,exist_ok=True)            
            factors = selected_factors['{} {}'.format(ct,factorType)]            

            xticklabels = ['Union',
                           '{}'.format('-'.join([i for i in factors[:1]])),
                           '{}'.format('-'.join([i for i in factors[1:2]])),
                           '{}'.format('-'.join([i for i in factors[2:3]])),
                           '{}'.format('-'.join([i for i in factors[:2]])),
                           '{}'.format('-'.join([i for i in factors[:3]])),
                           ]

            for ii in np.arange(len(xticklabels))[:]:
                tfbs_file = '{}/{}/{}_{}_{}.bed'.format(c_tfbs_dir,subdir,treat_flag,ct,xticklabels[ii])

                # se_file = '../../data/SE_hg38/{}.bed'.format(cancertype_SE)
                # == get TCGA peaks overlap/NOT-overlap with SE data
                overlap_file = '{}/{}/ATAC_overlap_{}_{}_{}.bed'.format(outdir,subdir,treat_flag,ct,xticklabels[ii])
                print('python {} -a {} -b {} -s hg38 -p {} -e2 0'.format(pyfile,tcga_file,tfbs_file,overlap_file))
                os.system('python {} -a {} -b {} -s hg38 -p {} -e2 0'.format(pyfile,tcga_file,tfbs_file,overlap_file))


  


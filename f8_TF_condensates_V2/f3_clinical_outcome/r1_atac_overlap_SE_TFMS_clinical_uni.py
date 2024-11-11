import os,sys,argparse,glob,re
import numpy as np
import pandas as pd
# import find_overlap_keep_info_NOT_sep_strand_asimport



# == main 
indir = 'data/ATAC_overlap_SE_overlap_TFMS'
outdir = 'f1_ATAC_overlap_SE_TFMS_clinical_uni'
os.makedirs(outdir+os.sep+'_csv',exist_ok=True)

project_dir='/Volumes/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/f7_TF_condensates_test'
# project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/f7_TF_condensates_test'
atac_overlap_SE_dir = '{}/f6_revised_TCGA_ATAC_cor_SE/f3_clinical_hicor_SE_overlapped_peaks/f1_clinical_at_each_peak/f1_ATAC_overlap_SE_caseID'.format(project_dir)
clinical_dir = '{}/f6_revised_TCGA_ATAC_cor_SE/f3_clinical_hicor_SE_overlapped_peaks/f1_clinical_at_each_peak/f2_caseID_each_SE_vs_clinical'.format(project_dir)
# atac_file = '{}/f6_revised_TCGA_ATAC_cor_SE/data/TCGA/tcga_atac.bed'.format(project_dir)

motif_file_dir='/nv/vol190/zanglab/shared/Motif/sites/hg38_fimo_jarspar/results/'
motif_file_dir='/Volumes/zanglab/shared/Motif/sites/hg38_fimo_jarspar/results/'
motif_infiles = glob.glob('{}/*'.format(motif_file_dir))
motif_infiles = [i for i in motif_infiles if not re.search('~',i)]
motif_infiles.sort()


# read matched names between TCGA ATAC and SE data
name_match = pd.read_excel('data/TCGA-ATAC_SE_cancerType_match.xlsx',index_col=0)   
name_match = name_match.dropna()
for cancertype in ['BRCA','COAD'][:]:
    cancertype_SE = name_match.loc[cancertype].SE
    cancertype_SE_rename = name_match.loc[cancertype].SE_rename
    clinical_file = '{}/{}_logrank_info.csv'.format(clinical_dir,cancertype)
    clinical_df = pd.read_csv(clinical_file,index_col=0)
    # atac_overlap_SE_file='{}/{}_ATAC_overlap_{}_SE_caseID.bed'.format(atac_overlap_SE_dir,cancertype,cancertype_SE)
    out_df = pd.DataFrame()
    for motif_file in motif_infiles[:]:
        motif_name = os.path.basename(motif_file).split('.bed')[0]
        overlapped_bed = '{}/{}_ATAC_overlap_{}_SE_overlap_{}.bed'.format(indir,cancertype,cancertype_SE,motif_name)
        try:
            overlapped_df = pd.read_csv(overlapped_bed,sep='\t',index_col=3)
            df = clinical_df.loc[overlapped_df.index]
            df.to_csv('{}/_csv/{}_overlap_SE_{}_ranksum.csv'.format(outdir,cancertype,motif_name))
            df_p = df[df['log rank p']<0.05]
            df_pt = df[(df['treat time']<df['ctrl time'])&(df['log rank p']<0.05)]
            out_df.loc[motif_name,'total'] = df.shape[0]
            out_df.loc[motif_name,'#P<0.05'] = df_p.shape[0]
            out_df.loc[motif_name,'%P<0.05'] = np.round(df_p.shape[0]/df.shape[0],2)
            out_df.loc[motif_name,'#TreatTime<CtrlTime'] = df_pt.shape[0]
            out_df.loc[motif_name,'%TreatTime<CtrlTime'] = np.round(df_pt.shape[0]/df.shape[0],2)
        except:
            print(motif_file)
    out_df.to_csv('{}/{}_clinical_summary.csv'.format(outdir,cancertype))
        
    
            

  


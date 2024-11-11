import os,sys,argparse,glob,re
import numpy as np
import pandas as pd
# import find_overlap_keep_info_NOT_sep_strand_asimport



# == main 
indir = 'ATAC_overlap_SE_overlap_TFMS'
outdir = 'ATAC_overlap_SE_overlap_TFMS_avg_sig'
os.makedirs(outdir,exist_ok=True)


project_dir='/Volumes/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/f7_TF_condensates_test'
# project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/f7_TF_condensates_test'
atac_overlap_SE_dir = '{}/f6_revised_TCGA_ATAC_cor_SE/f3_clinical_hicor_SE_overlapped_peaks/f1_clinical_at_each_peak/f1_ATAC_overlap_SE_caseID'.format(project_dir)
# atac_file = '{}/f6_revised_TCGA_ATAC_cor_SE/data/TCGA/tcga_atac.bed'.format(project_dir)

motif_file_dir='/nv/vol190/zanglab/shared/Motif/sites/hg38_fimo_jarspar/results/'
motif_file_dir='/Volumes/zanglab/shared/Motif/sites/hg38_fimo_jarspar/results/'
motif_infiles = glob.glob('{}/*'.format(motif_file_dir))
motif_infiles = [i for i in motif_infiles if not re.search('~',i)]
motif_infiles.sort()

filtered_df = pd.read_excel('{}/f6_revised_TCGA_ATAC_cor_SE//data/TCGA/TCGA-ATAC_clustered_samples.xlsx'.format(project_dir),index_col=0)   


# read matched names between TCGA ATAC and SE data
name_match = pd.read_excel('TCGA-ATAC_SE_cancerType_match.xlsx',index_col=0)   
name_match = name_match.dropna()
for cancertype in ['BRCA','COAD'][:]:
    cancertype_SE = name_match.loc[cancertype].SE
    cancertype_SE_rename = name_match.loc[cancertype].SE_rename
    atac_overlap_SE_file='{}/{}_ATAC_overlap_{}_SE_caseID.bed'.format(atac_overlap_SE_dir,cancertype,cancertype_SE)
    atac_overlap_SE_df = pd.read_csv(atac_overlap_SE_file,sep='\t',index_col=3)
    filtered_id = filtered_df[filtered_df.cohort==cancertype].case_id
    
    # read bed overlapp w/ each factor motif
    for motif_file in motif_infiles[:]:
        motif_name = os.path.basename(motif_file).split('.bed')[0]
        try:
            overlapped_bed = '{}/{}_ATAC_overlap_{}_SE_overlap_{}.bed'.format(indir,cancertype,cancertype_SE,motif_name)
            overlapped_df = pd.read_csv(overlapped_bed,sep='\t',index_col=3,header=None)
            overlapped_df.columns = atac_overlap_SE_df.columns
            
            bed_df = overlapped_df.iloc[:,:6]
            sig_df = overlapped_df.iloc[:,6:]
            sig_df.columns = [i.split('_')[1] for i in sig_df.columns]
            bed_df['avg-score'] = sig_df[filtered_id].mean(axis=1)
            bed_df.insert(3,'name',bed_df.index)
            bed_df = bed_df[['#seqnames', 'start', 'end', 'name', 'avg-score', 'annotation']]
            bed_df.to_csv('{}/{}_ATAC_overlap_{}_SE_overlap_{}.bed'.format(outdir,cancertype,cancertype_SE,motif_name),sep='\t',index=False,header=None)
        except:
            print(motif_name)

  


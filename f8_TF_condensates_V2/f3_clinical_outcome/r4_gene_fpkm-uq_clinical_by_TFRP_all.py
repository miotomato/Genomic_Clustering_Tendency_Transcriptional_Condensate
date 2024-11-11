import os,sys,argparse,glob,re
import numpy as np
import pandas as pd
# import find_overlap_keep_info_NOT_sep_strand_asimport
from collections import Counter
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=14
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams["mathtext.rm"] = "Arial"
import seaborn as sns
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks")
from scipy import stats




def add_logrank_info(out_df,df,motif_name):

    df_p = df[df['log rank p']<0.05]
    df_pt = df[(df['treat time']<df['ctrl time'])&(df['log rank p']<0.05)]
    out_df.loc[motif_name,'total'] = df.shape[0]
    out_df.loc[motif_name,'#P<0.05'] = df_p.shape[0]
    out_df.loc[motif_name,'%P<0.05'] = np.round(df_p.shape[0]/df.shape[0],2)
    out_df.loc[motif_name,'#TreatTime<CtrlTime&P<0.05'] = df_pt.shape[0]
    out_df.loc[motif_name,'%TreatTime<CtrlTime&P<0.05'] = np.round(df_pt.shape[0]/df.shape[0],2)
    return out_df



# ==== main 

project_dir='/Volumes/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/f7_TF_condensates_test'
# project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/f7_TF_condensates_test'
clinical_dir = '{}/f6_revised_TCGA_ATAC_cor_SE/f5_clinical_per_gene/f1_clinical_per_gene/'.format(project_dir)


outdir = 'f4_gene_fpkm-uq_clinical_by_TFRP_all'
os.makedirs(outdir,exist_ok=True)
rp_dir = 'data/ATAC_overlap_SE_overlap_TFMS_avg_sig_RP'

# ==== read matched names between TCGA ATAC and SE data
name_match = pd.read_excel('data/TCGA-ATAC_SE_cancerType_match.xlsx',index_col=0)   
name_match = name_match.dropna()


motif_file_dir='/nv/vol190/zanglab/shared/Motif/sites/hg38_fimo_jarspar/results/'
motif_file_dir='/Volumes/zanglab/shared/Motif/sites/hg38_fimo_jarspar/results/'
motif_infiles = glob.glob('{}/*'.format(motif_file_dir))
motif_infiles = [i for i in motif_infiles if not re.search('~',i)]
motif_infiles.sort()


for data_type in ['SE','CP'][:1]:
    for cancertype in ['BRCA','COAD'][:]:
        cancertype_SE = name_match.loc[cancertype].SE
        cancertype_SE_rename = name_match.loc[cancertype].SE_rename
        # atac_overlap_SE_file='{}/{}_ATAC_overlap_{}_SE_caseID.bed'.format(atac_overlap_SE_dir,cancertype,cancertype_SE)
        # box_vals,xticklabels = [],[]
        clinical_file = '{}/{}_logrank_info.csv'.format(clinical_dir,cancertype)
        clinical_df = pd.read_csv(clinical_file,index_col=0)
        
        # == save the log rank P info of each TF
        out_df = pd.DataFrame()
        for motif_file in motif_infiles[:]:
            motif_name = os.path.basename(motif_file).split('.bed')[0]
            out_file = '{}/{}_ATAC_overlap_{}_SE_overlap_{}_RP.tsv'.format(rp_dir,cancertype,cancertype_SE,motif_name)
            
            try:
                rp_df = pd.read_csv(out_file,sep='\t',index_col=0,header=None)
                rp_df.columns=['RP']
                rp_df = rp_df.sort_values(by='RP',ascending=False)                
                selected_df = clinical_df.loc[rp_df[rp_df.RP>1].index].dropna()
                out_df = add_logrank_info(out_df,selected_df,motif_name)
            
            except:
                print(data_type,cancertype,motif_name)
                
        out_df.to_csv('{}/{}_TF_targets_RP_GT1.csv'.format(outdir,cancertype))
        
    



            

  


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
sns.set(font_scale=1.1)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks")
from scipy import stats





def add_logrank_info(out_df,df,motif_name,flag):

    df_p = df[df['log rank p']<p_thre]
    # if cp_fig:
    #     for gene in df_p.index[:10]:
    #         fig_file = '{}/fig/{}_{}_survival.pdf'.format(clinical_dir,cancertype,gene)
    #         os.system('cp {} {}/_fig/'.format(fig_file,outdir))
    df_pt = df[(df['treat time']<df['ctrl time'])&(df['log rank p']<p_thre)]
    out_df.loc[motif_name,'{} total'.format(flag)] = df.shape[0]
    out_df.loc[motif_name,'{} #P<{}'.format(flag,p_thre)] = df_p.shape[0]
    out_df.loc[motif_name,'{} %P<{}'.format(flag,p_thre)] = np.round(df_p.shape[0]/df.shape[0],2)
    out_df.loc[motif_name,'{} #TreatTime<CtrlTime&P<{}'.format(flag,p_thre)] = df_pt.shape[0]
    out_df.loc[motif_name,'{} %TreatTime<CtrlTime&P<{}'.format(flag,p_thre)] = np.round(df_pt.shape[0]/df.shape[0],2)
    return out_df




# ==== main 

project_dir='/Volumes/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/'
# project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR'
# atac_overlap_SE_dir = '{}/f7_TF_condensates_test/f6_revised_TCGA_ATAC_cor_SE/f3_clinical_hicor_SE_overlapped_peaks/f1_clinical_at_each_peak/f1_ATAC_overlap_SE_caseID'.format(project_dir)
clinical_per_peak_dir = '{}/f7_TF_condensates_test/f6_revised_TCGA_ATAC_cor_SE/f3_clinical_hicor_SE_overlapped_peaks/f1_clinical_at_each_peak/f2_caseID_each_SE_vs_clinical/'.format(project_dir)
clinical_per_gene_dir = '{}/f7_TF_condensates_test/f6_revised_TCGA_ATAC_cor_SE/f5_clinical_per_gene/f1_clinical_per_gene/'.format(project_dir)
# atac_file = '{}/f7_TF_condensates_test/f6_revised_TCGA_ATAC_cor_SE/data/TCGA/tcga_atac.bed'.format(project_dir)
# atac_df = pd.read_csv(atac_file,sep='\t',header=None,index_col=3)
atac_overlap_SE_TFMS_dir = '{}/f8_TF_condensates_V2/f3_clinical_outcome/data/ATAC_overlap_SE_overlap_TFMS_avg_sig/'.format(project_dir)
rp_dir = '{}/f8_TF_condensates_V2/f3_clinical_outcome/data/ATAC_overlap_SE_overlap_TFMS_avg_sig_RP'.format(project_dir)

# ==== read matched names between TCGA ATAC and SE data
name_match = pd.read_excel('../../data/TCGA/TCGA-ATAC_SE_cancerType_match.xlsx',index_col=0)   
name_match = name_match.dropna()
# ==== read TFBS CP info
tfbs_dir = '../../f1_TF_cluster_potential/f2_cor_CP_SE_AICAP/f3_TFBS_CP_heatmap/_csv/'
tfbs_cp_type = 'CP_TFBS_vs_TFMS'
# cp_type = 'TFBS_CP'
tf_thre = 6
rp_thre = 1
p_thre = 0.05


# outdir = 'f2_clinical_gene_vs_peak_gene_treatLTcontrol'
outdir = 'f2_clinical_gene_vs_peak_gene'
os.makedirs(outdir,exist_ok=True)


for cancertype in ['BRCA','COAD'][:]:
    cancertype_SE = name_match.loc[cancertype].SE
    cancertype_SE_rename = name_match.loc[cancertype].SE_rename
    
    # == TFBS CP
    tfbs_cp_s = pd.read_csv('{}/data_fisher_{}_CP_RankSum_statistics.csv'.format(tfbs_dir,tfbs_cp_type),index_col=0)
    tfbs_se_s = pd.read_csv('{}/data_fisher_{}_SE_statistics.csv'.format(tfbs_dir,tfbs_cp_type),index_col=0)

    
    # == save the log rank P info of each gene/peak
    out_df = pd.DataFrame()
    
    # == logrank pvalue of SE overlapped peaks 
    clinical_per_peak_file = '{}/{}_logrank_info.csv'.format(clinical_per_peak_dir,cancertype)
    clinical_per_peak_df = pd.read_csv(clinical_per_peak_file,index_col=0)
    out_df = add_logrank_info(out_df,clinical_per_peak_df,'All','SE-overlapped-peaks')
    
    # == logrank pvalue of all genes 
    clinical_per_gene_file = '{}/{}_logrank_info.csv'.format(clinical_per_gene_dir,cancertype)
    clinical_per_gene_df = pd.read_csv(clinical_per_gene_file,index_col=0)
    # clinical_per_gene_df = clinical_per_gene_df[(clinical_per_gene_df['treat time']<clinical_per_gene_df['ctrl time'])]
    out_df = add_logrank_info(out_df,clinical_per_gene_df,'All','Target-genes')
    
    # ==== target-gene and peak- clinical info for each TF
    for motif_name in tfbs_cp_s.index[:]:
        
        # ==== add target gene clinical info
        rp_files = glob.glob('{}/{}_ATAC_overlap_{}_SE_overlap_{}_*.tsv'.format(rp_dir,cancertype,cancertype_SE,motif_name))
        rp_files = [i for i in rp_files if not re.search('~',i)]
        assert len(rp_files)==1
        rp_df = pd.read_csv(rp_files[0],sep='\t',index_col=0,header=None)
        rp_df.columns=['RP']
        # ==== get genes with RP>1, read the gene clinical info
        df = clinical_per_gene_df.loc[rp_df[rp_df.RP>rp_thre].index].dropna();
        out_df = add_logrank_info(out_df,df,motif_name,'Target-genes')
        
        
        # ==== add SE-overlapped-peak clinical info
        overlapped_beds = glob.glob('{}/{}_ATAC_overlap_{}_SE_overlap_{}_*.bed'.format(atac_overlap_SE_TFMS_dir,cancertype,cancertype_SE,motif_name))
        overlapped_beds = [i for i in overlapped_beds if not re.search('~',i)]
        assert len(overlapped_beds)==1
        overlapped_df = pd.read_csv(overlapped_beds[0],sep='\t',index_col=3,header=None)
        df = clinical_per_peak_df.loc[overlapped_df.index]
        out_df = add_logrank_info(out_df,df,motif_name,'SE-overlapped-peaks')
        
    # ==== add TFBS CP info   
    out_df['avg %P<0.05'] = (out_df['SE-overlapped-peaks %P<0.05']+out_df['Target-genes %P<0.05'])/2
    out_df = pd.concat([out_df,tfbs_cp_s[cancertype_SE].rename('TFBS CP')],axis=1)
    out_df = pd.concat([out_df,tfbs_se_s[cancertype_SE].rename('SE enrichment log2 OR')],axis=1)
    out_df.to_csv('{}/{}_clinical_Gene_vs_Peaks.csv'.format(outdir,cancertype))
                
                

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





def add_logrank_info(out_df,df,motif_name,cp_fig=False):

    df_p = df[df['log rank p']<p_thre]
    # if cp_fig:
    #     for gene in df_p.index[:10]:
    #         fig_file = '{}/fig/{}_{}_survival.pdf'.format(clinical_dir,cancertype,gene)
    #         os.system('cp {} {}/_fig/'.format(fig_file,outdir))
    df_pt = df[(df['treat time']<df['ctrl time'])&(df['log rank p']<p_thre)]
    out_df.loc[motif_name,'total'] = df.shape[0]
    out_df.loc[motif_name,'#P<{}'.format(p_thre)] = df_p.shape[0]
    out_df.loc[motif_name,'%P<{}'.format(p_thre)] = np.round(df_p.shape[0]/df.shape[0],2)
    out_df.loc[motif_name,'#TreatTime<CtrlTime&P<{}'.format(p_thre)] = df_pt.shape[0]
    out_df.loc[motif_name,'%TreatTime<CtrlTime&P<{}'.format(p_thre)] = np.round(df_pt.shape[0]/df.shape[0],2)
    return out_df




# ==== main 
indir = 'f4_clinical_per_gene_by_TFBS_RP'
outdir = 'f5_TF_targets_clinical_by_gene_TFBS_RP'
os.makedirs(outdir,exist_ok=True)
# os.makedirs(outdir+,exist_ok=True)

project_dir='/Volumes/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/'
# project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR'

# ==== read matched names between TCGA ATAC and SE data
name_match = pd.read_excel('{}/f8_TF_condensates_V2/f3_clinical_outcome/data/TCGA-ATAC_SE_cancerType_match.xlsx'.format(project_dir),index_col=0)   
name_match = name_match.dropna()

# get top3 factors
selected_factors = {}
tfbs_cp_dir = '{}/f12_KS_test_Rename/f1_TF_cluster_potential/f2_cor_CP_SE_AICAP/f9_per_CT_TFBS_CP_cor_zscore_CP_with_motif_SE/TFBS_CP/'.format(project_dir)
for ct in ['MCF-7','HCT-116'][:]:
    df = pd.read_csv('{}/_CP_TFBS_all_vs_TFMS_{}.csv'.format(tfbs_cp_dir,ct),index_col=0)
    selected_factors['{} top_TFBSCP'.format(ct)] = df['TFBS CP rank'].sort_values().iloc[:3].index
    selected_factors['{} top_zscored_TFBSCP'.format(ct)] = df['avg rank'].sort_values().iloc[:3].index


p_thres = [0.1,0.05]
for p_thre in p_thres:
    for cancertype in ['BRCA','COAD'][:1]:
        # ==== read clinical data
        ct = name_match.loc[cancertype].SE
        # cancertype_SE_rename = name_match.loc[cancertype].SE_rename
        for treat_flag in ['percentile_T','percentile_T_ExtendMerge'][1:]:
            for factorType in ['top_TFBSCP','top_zscored_TFBSCP'][:]:
                subdir = '{}_{}'.format(ct,factorType)
                os.makedirs(outdir+os.sep+subdir+os.sep+'fig',exist_ok=True)            
                factors = selected_factors['{} {}'.format(ct,factorType)]            
    
                xticklabels = ['ALL','Union',
                               '{}'.format('-'.join([i for i in factors[:1]])),
                               '{}'.format('-'.join([i for i in factors[1:2]])),
                               '{}'.format('-'.join([i for i in factors[2:3]])),
                               '{}'.format('-'.join([i for i in factors[:2]])),
                               '{}'.format('-'.join([i for i in factors[:3]])),
                               ]
    
                out_df = pd.DataFrame()
                for ii in np.arange(len(xticklabels))[:]:
                    basename = '{}_{}_{}'.format(treat_flag,ct,xticklabels[ii])           
                    clinical_file = '{}/{}/{}_logrank_info.csv'.format(indir,subdir,basename)
                    clinical_df = pd.read_csv(clinical_file,index_col=0)
                    # clinical_df = clinical_df[(clinical_df['treat time']<clinical_df['ctrl time'])]
                    out_df = add_logrank_info(out_df,clinical_df,xticklabels[ii])                
    
                out_df = out_df.drop(['Union'])
                check_col = '#P<{}'.format(p_thre)
                # check_col = '#TreatTime<CtrlTime&P<{}'.format(p_thre)

                plt.figure(figsize=(3,2.6))
                positions = np.arange(out_df.shape[0])                
                t_all = out_df.loc['ALL','total']
                p_all = out_df.loc['ALL',check_col]

                for position in positions:
                    total = out_df.iloc[position]['total']
                    a = out_df.iloc[position][check_col]
                    s,p = stats.fisher_exact([[a,total-a],[p_all-a,t_all-p_all-total+a]])
                    p_label='{:.1e}'.format(p)
                    out_df.loc[out_df.index[position],'{} fisher s'.format(check_col)] = s
                    out_df.loc[out_df.index[position],'{} fisher p'.format(check_col)] = p
                    # if p_label[-2]=='0':
                    #     p_label = p_label[:-2]+p_label[-1]
                    if p<0.05:
                        # plt.text(position,100*a/total, p_label , ha='center', va='bottom', color='k',fontsize=10)
                        plt.text(position,92*a/total, '*' , ha='center', va='bottom', color='tab:red',fontsize=20)
    
                    g0 = plt.bar(position,100*a/total,bottom=0,width = .68, 
                                 lw=0,color = 'tab:purple',alpha=.6)
                            
                plt.title('{}'.format(cancertype))
                plt.axes().set_xticks(positions) 
                plt.axes().set_xticklabels(out_df.index,rotation=45, ha='right',fontsize=13,color='k')
                plt.ylabel('% ATAC-seq target genes \n w/ logrank P<{}'.format(p_thre),fontsize=13)
                figname = '{}/{}/{}_{}_p{}.pdf'.format(outdir,subdir,treat_flag,ct,p_thre)
                plt.savefig(figname,bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
                plt.show()
                plt.close()

                out_df.to_csv('{}/{}/{}_{}_p{}_data.csv'.format(outdir,subdir,treat_flag,ct,p_thre))
                
            
            # copy fig
            SE_overlapped_promoter_df = pd.read_csv('data_gene_overlap_SE/promoter_on_MCF7.bed',sep='\t',header=None,index_col=3)
            clinical_df = clinical_df.sort_values(by=['log rank p'],ascending=True)
            # target_genes = clinical_df[clinical_df['log rank p']<0.1].index.intersection(SE_overlapped_promoter_df.index)
            target_genes = clinical_df[(clinical_df['log rank p']<0.1)&(clinical_df['treat high RP avg']>20)].index
            for ii in target_genes:
                pdf_file = '{}/{}/fig/{}_{}_{}_{}_survival.pdf'.format(indir,subdir,treat_flag,ct,xticklabels[-1],ii)
                commandLine = 'cp {} {}/{}/fig'.format(pdf_file,outdir,subdir)
                os.system(commandLine)


                
            
        
        

            

  


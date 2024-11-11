import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import json
#from GenomeData import *
import re,bisect
from lifelines.statistics import logrank_test
from lifelines import KaplanMeierFitter

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=14
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=1.1)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style('ticks')

import warnings
warnings.filterwarnings("ignore")



def fdr_adj_p(pvalues,p_index):
    df = pvalues.to_frame()
    n,k = len(pvalues),len(pvalues)      
    minimum = 1    
    for index in df.sort_values(by=p_index,ascending=False).index:
        pvalue = df.loc[index,p_index]
        fdr = n*pvalue/k  
        minimum = min(minimum,fdr)
        df.loc[index,'fdr'] = minimum
        k=k-1     
    return df['fdr']


    
def survival_for_two(df,treat,ctrl,legends,title,figname):
    
    # select the time and status info for treat and control group
    ix = df['group'] == treat
    t1 = df.loc[ix]['time']#;print(t1.shape)
    e1 = df.loc[ix]['status'] 
    t2 = df.loc[~ix]['time']#;print(t2.shape)
    e2 = df.loc[~ix]['status']
    
    results = logrank_test(t1,t2,event_observed_A = e1,event_observed_B = e2)
    pvalue = results.p_value#;print('pvalue:\t{}'.format(pvalue))
    
    # if 0:
    if pvalue<0.05:
        # survival curves
        plt.figure(figsize=(2.6,2.6))
        ax = plt.subplot(111)
        
        kmf_control = KaplanMeierFitter()
        g1 = kmf_control.fit(t1, e1).plot(ax=ax,show_censors=True,\
                                          label='more open',\
                            censor_styles={'ms': 12, 'marker': '+'},ci_show=False,color='red',ls='-')
        
        kmf_exp = KaplanMeierFitter()
        g2 = kmf_exp.fit(t2, e2).plot(ax=ax,show_censors=True,\
                                      label='less open',\
                        censor_styles={'ms': 12, 'marker': '+'},ci_show=False,color='k',ls='--')
        
        handles, labels = ax.get_legend_handles_labels()
        lg = ax.legend(handles, legends,loc='lower left',borderaxespad=0,handletextpad=.2,labelspacing=.2,
                       handlelength=1,frameon=False)
        if pvalue<1:
              plt.axes().text(df['time'].max()*0.45,0.45,'p={:.1e}'.format(pvalue),fontsize=16,ha='center')
        plt.ylim([0.0,1.0])
        plt.title('{}'.format(title),fontsize=12)
        plt.xlabel('Days',fontsize=14)
        plt.ylabel('Survival probability',fontsize=14)
        plt.savefig(figname,bbox_inches='tight',pad_inches=.1,dpi=600,transparent=True)
        # plt.show()
        plt.close()
    return results




    
# ==== main  
# indir = 'f1_ATAC_overlap_SE_caseID'
outdir = 'f4_clinical_per_gene_by_TFBS_RP'
os.makedirs(outdir,exist_ok=True)
# os.makedirs(outdir+os.sep+'fig',exist_ok=True)

project_dir='/Volumes/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/'
project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR'

name_match = pd.read_excel('{}/f8_TF_condensates_V2/f3_clinical_outcome/data/TCGA-ATAC_SE_cancerType_match.xlsx'.format(project_dir),index_col=0)   
name_match = name_match.dropna()
clinical_dir = '{}/f7_TF_condensates_test/f6_revised_TCGA_ATAC_cor_SE/f3_clinical_hicor_SE_overlapped_peaks/f1_clinical_at_each_peak/f2_caseID_each_SE_vs_clinical'.format(project_dir) 
tss_df = pd.read_csv('{}/../../data//geneID_annotation/hg38/hg38_4k_promoter_geneID.bed'.format(project_dir),sep='\t',index_col=3)

# get top3 factors
selected_factors = {}
tfbs_cp_dir = '{}/f12_KS_test_Rename/f1_TF_cluster_potential/f2_cor_CP_SE_AICAP/f9_per_CT_TFBS_CP_cor_zscore_CP_with_motif_SE/TFBS_CP/'.format(project_dir)
for ct in ['MCF-7','HCT-116'][:]:
    df = pd.read_csv('{}/_CP_TFBS_all_vs_TFMS_{}.csv'.format(tfbs_cp_dir,ct),index_col=0)
    selected_factors['{} top_TFBSCP'.format(ct)] = df['TFBS CP rank'].sort_values().iloc[:3].index
    selected_factors['{} top_zscored_TFBSCP'.format(ct)] = df['avg rank'].sort_values().iloc[:3].index


for cancertype in ['BRCA','COAD'][:]:
    # ==== read clinical data
    ct = name_match.loc[cancertype].SE
    # cancertype_SE_rename = name_match.loc[cancertype].SE_rename
    for treat_flag in ['percentile_T','percentile_T_ExtendMerge'][:]:
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

            for ii in np.arange(len(xticklabels))[:]:
                print(ct,treat_flag,xticklabels[ii])
                basename = '{}_{}_{}'.format(treat_flag,ct,xticklabels[ii])           
                rp_file = 'data_TFBS_RP/ATAC_overlap_TFBS_per_patient_sig_RP_combined/{}/ATAC_overlap_{}_RP_combined.tsv'.format(subdir,basename)
                rp_df = pd.read_csv(rp_file,index_col=0,)
                rp_df = rp_df.loc[rp_df.index.intersection(tss_df.index)]
                rp_df = rp_df[~rp_df.index.duplicated()]
              
                clinical_df = pd.read_csv('{}/{}_clinical_info.csv'.format(clinical_dir,cancertype),index_col=0)
                clinical_df['time'] = clinical_df[['days_to_death','days_to_last_follow_up']].fillna(0).astype(float).max(axis=1)   
                clinical_df = clinical_df[clinical_df['time']>0]
                clinical_df['status']= clinical_df['vital_status'] 
                clinical_df.loc[clinical_df['status']=='Dead','status']=1
                clinical_df.loc[clinical_df['status']=='Alive','status']=0
                clinical_df.to_csv('{}/{}/{}_clinical_info.csv'.format(outdir,subdir,basename))
                
                # patient samples per cancertype
                kept_col = clinical_df.index.intersection(rp_df.columns)
                clinical_df = clinical_df.loc[kept_col]
                rp_df_tmp = rp_df[kept_col]
                rp_df_tmp = np.transpose(rp_df_tmp)
                
                # # ==== save the data
                out_df = pd.concat([clinical_df,rp_df_tmp],axis=1)
                out_df.to_csv(outdir+os.sep+'{}/{}_RP_add_clinical.csv'.format(subdir,basename))
            
            
                
                # == survival per gene
                patients_each_group = int(len(kept_col)*0.5)
                region_outdf = pd.DataFrame()
                for region in rp_df_tmp.columns[:]:
                    df_tmp = out_df[['time','status',region]]
                    df_tmp = df_tmp.sort_values(by=region,ascending=False)
                    
                    # separate the patients into two groups
                    treat_index = df_tmp.index[:patients_each_group]
                    ctrl_index = df_tmp.index[-1*patients_each_group:]
                    df_tmp = df_tmp.loc[treat_index.union(ctrl_index)]
                    df_tmp.loc[treat_index,'group']='treat'   
                    df_tmp.loc[ctrl_index,'group']='ctrl'  
                    # ==== plot the figs
                    figname=outdir+os.sep+'{}/fig/{}_{}_survival.pdf'.format(subdir,basename,region)
                    legends = ['high RP','low RP']
                    results = survival_for_two(df_tmp,'treat','ctrl',legends,region,figname)
                    if results.p_value<0.05:
                        df_tmp.to_csv(outdir+os.sep+'{}/fig/{}_{}_survival.csv'.format(subdir,basename,region))
            
                    # save info
                    ix = df_tmp['group'] == 'treat'
                    region_outdf.loc[region,'treat high RP avg'] = df_tmp.loc[ix,region].values.mean().round(2)
                    region_outdf.loc[region,'ctrl low RP avg'] = df_tmp.loc[~ix,region].values.mean().round(2)
                    region_outdf.loc[region,'treat time'] = df_tmp.loc[ix,'time'].values.mean().round(2)
                    region_outdf.loc[region,'ctrl time'] = df_tmp.loc[~ix,'time'].values.mean().round(2)
                    # region_outdf.loc[region,'log rank s'] = results.test_statistic  
                    region_outdf.loc[region,'log rank p'] = results.p_value
                region_outdf['fdr'] = fdr_adj_p(region_outdf['log rank p'],'log rank p')
                region_outdf = region_outdf.sort_values(by='log rank p',ascending=True)
                region_outdf.to_csv('{}/{}/{}_logrank_info.csv'.format(outdir,subdir,basename))
            
                     
                        

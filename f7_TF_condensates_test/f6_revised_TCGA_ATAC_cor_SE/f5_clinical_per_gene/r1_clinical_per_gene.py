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
        
        handles, labels = ax.get_legend_handles_labels();print(labels)
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
outdir = 'f1_clinical_per_gene'
os.makedirs(outdir+os.sep+'fig',exist_ok=True)
# os.makedirs(outdir+os.sep+'fig',exist_ok=True)


clinical_dir = '../f3_clinical_hicor_SE_overlapped_peaks/f1_clinical_at_each_peak/f2_caseID_each_SE_vs_clinical' 
tss_df = pd.read_csv('../../../../../data//geneID_annotation/hg38/hg38_4k_promoter_geneID.bed',sep='\t',index_col=3)
expr_df = pd.read_csv('../data/TCGA/TCGA-ATAC_clustered_samples.htseq_fpkm-uq.tsv',sep='\t',index_col=0)
expr_df = expr_df.loc[expr_df.index.intersection(tss_df.index)]
expr_df = expr_df[~expr_df.index.duplicated()]

for cancertype in ['BRCA','COAD'][:]:
    # ==== read clinical data
    clinical_df = pd.read_csv('{}/{}_clinical_info.csv'.format(clinical_dir,cancertype),index_col=0)
    clinical_df['time'] = clinical_df[['days_to_death','days_to_last_follow_up']].fillna(0).astype(float).max(axis=1)   
    clinical_df = clinical_df[clinical_df['time']>0]
    clinical_df['status']= clinical_df['vital_status'] 
    clinical_df.loc[clinical_df['status']=='Dead','status']=1
    clinical_df.loc[clinical_df['status']=='Alive','status']=0
    
    # patient samples per cancertype
    kept_col = clinical_df.index.intersection(expr_df.columns)
    clinical_df = clinical_df.loc[kept_col]
    expr_df_tmp = expr_df[kept_col]
    expr_df_tmp = np.transpose(expr_df_tmp)
    
    # ==== save the data
    out_df = pd.concat([clinical_df,expr_df_tmp],axis=1)
    out_df.to_csv(outdir+os.sep+'{}_gene_fpkm-uq_add_clinical.csv'.format(cancertype))


    # == survival per gene
    patients_each_group = int(len(kept_col)*0.5)
    region_outdf = pd.DataFrame()
    for region in expr_df_tmp.columns[:]:
        df_tmp = out_df[['time','status',region]]
        df_tmp = df_tmp.sort_values(by=region,ascending=False)
        
        # separate the patients into two groups
        treat_index = df_tmp.index[:patients_each_group]
        ctrl_index = df_tmp.index[-1*patients_each_group:]
        df_tmp = df_tmp.loc[treat_index.union(ctrl_index)]
        df_tmp.loc[treat_index,'group']='treat'   
        df_tmp.loc[ctrl_index,'group']='ctrl'  
        # ==== plot the figs
        figname=outdir+os.sep+'/fig/{}_{}_survival.pdf'.format(cancertype,region)
        legends = ['high fpkm','low fpkm']
        results = survival_for_two(df_tmp,'treat','ctrl',legends,region,figname)
        if results.p_value<0.05 :
            df_tmp.to_csv(outdir+os.sep+'fig/{}_{}_survival.csv'.format(cancertype,region))

        # save info
        ix = df_tmp['group'] == 'treat'
        region_outdf.loc[region,'treat high fpkm avg'] = df_tmp.loc[ix,region].values.mean().round(2)
        region_outdf.loc[region,'ctrl low fpkm avg'] = df_tmp.loc[~ix,region].values.mean().round(2)
        region_outdf.loc[region,'treat time'] = df_tmp.loc[ix,'time'].values.mean().round(2)
        region_outdf.loc[region,'ctrl time'] = df_tmp.loc[~ix,'time'].values.mean().round(2)
        # region_outdf.loc[region,'log rank s'] = results.test_statistic  
        region_outdf.loc[region,'log rank p'] = results.p_value
    region_outdf['fdr'] = fdr_adj_p(region_outdf['log rank p'],'log rank p')
    region_outdf = region_outdf.sort_values(by='log rank p',ascending=True)
    region_outdf.to_csv('{}/{}_logrank_info.csv'.format(outdir,cancertype))

         
            

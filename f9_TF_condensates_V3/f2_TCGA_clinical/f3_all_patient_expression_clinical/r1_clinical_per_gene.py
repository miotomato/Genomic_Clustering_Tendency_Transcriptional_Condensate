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
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=12
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=1)
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


    
def survival_for_two(df,treat,ctrl,legends,title,figname,pvalue_thre):
    
    # select the time and status info for treat and control group
    ix = df['group'] == treat
    t1 = df.loc[ix]['time']#;print(t1.shape)
    e1 = df.loc[ix]['status'] 
    t2 = df.loc[~ix]['time']#;print(t2.shape)
    e2 = df.loc[~ix]['status']
    
    results = logrank_test(t1,t2,event_observed_A = e1,event_observed_B = e2)
    pvalue = results.p_value#;print('pvalue:\t{}'.format(pvalue))
    
    # if 0:
    if pvalue<pvalue_thre:
        # survival curves
        plt.figure(figsize=(2.2,2.2))
        ax = plt.subplot(111)
        
        kmf_control = KaplanMeierFitter()
        g1 = kmf_control.fit(t1, e1).plot(ax=ax,show_censors=True,\
                                          label=legends[0],\
                            censor_styles={'ms': 8, 'marker': '+'},ci_show=False,color='red',ls='-',lw=1)
        
        kmf_exp = KaplanMeierFitter()
        g2 = kmf_exp.fit(t2, e2).plot(ax=ax,show_censors=True,\
                                      label=legends[1],\
                        censor_styles={'ms': 8, 'marker': '+'},ci_show=False,color='k',ls='--',lw=1)
        
        handles, labels = ax.get_legend_handles_labels()
        # print(labels)
        lg = ax.legend(handles, legends,loc='lower left',borderaxespad=0,handletextpad=.2,labelspacing=.2,
                       handlelength=1,frameon=False,fontsize=11)
        if pvalue<1:
              plt.axes().text(df['time'].max()*0.75,0.9,'p={:.1e}'.format(pvalue),
                              fontsize=12,ha='center')
        plt.ylim([0.0,1.0])
        plt.title('{}'.format(title),fontsize=12)
        plt.xlabel('Days',fontsize=12)
        plt.ylabel('Survival probability',fontsize=12)
        plt.savefig(figname,bbox_inches='tight',pad_inches=.1,dpi=600,transparent=True)
        # plt.show()
        plt.close()
    return results




    
# ==== main  
# indir = 'f1_ATAC_overlap_SE_caseID'
outdir = 'f1_clinical_per_gene'
os.makedirs(outdir+os.sep+'fig',exist_ok=True)
# os.makedirs(outdir+os.sep+'fig',exist_ok=True)
tcga_data_dir = '../../data/TCGA/'


for cancertype in ['BRCA','COAD'][:]:
    
    # expr_df = pd.read_csv('{}/GDC-PANCAN.htseq_fpkm-uq.{}.head100.tsv'.format(tcga_data_dir,cancertype),sep='\t',index_col=0)
    expr_df = pd.read_csv('{}/GDC-PANCAN.htseq_fpkm-uq.{}.tsv'.format(tcga_data_dir,cancertype),sep='\t',index_col=0)    
    clinical_df = pd.read_csv('{}/clinical.project-TCGA-{}.2022-03-20.csv'.format(tcga_data_dir,cancertype),index_col=1)
    clinical_df = clinical_df.loc[expr_df.columns][['time','status']]
    
    # == survival per gene
    pvalue_thre = 0.001
    patients_each_group = int(clinical_df.shape[0]*0.5)
    
    region_outdf = pd.DataFrame()
    for region in expr_df.index[:]:
        # print(region)
        df_tmp = pd.concat([clinical_df,expr_df.loc[region]],axis=1)
        df_tmp = df_tmp.sort_values(by=region,ascending=False)
        
        # separate the patients into two groups
        treat_index = df_tmp.index[:patients_each_group]
        ctrl_index = df_tmp.index[-1*patients_each_group:]
        df_tmp = df_tmp.loc[treat_index.union(ctrl_index)]
        df_tmp.loc[treat_index,'group']='treat'   
        df_tmp.loc[ctrl_index,'group']='ctrl'  
        # ==== plot the figs
        figname=outdir+os.sep+'/fig/{}_{}_survival.pdf'.format(cancertype,region)
        legends = ['high fpkm-uq','low fpkm-uq']
        results = survival_for_two(df_tmp,'treat','ctrl',legends,region,figname,pvalue_thre)
        if results.p_value<pvalue_thre:
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

         
            

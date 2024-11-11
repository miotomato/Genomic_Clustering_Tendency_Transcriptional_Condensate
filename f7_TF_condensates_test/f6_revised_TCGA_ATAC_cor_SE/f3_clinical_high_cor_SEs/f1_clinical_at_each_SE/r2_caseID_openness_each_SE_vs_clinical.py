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
matplotlib.rcParams['font.size']=16
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=1.2)
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


def return_survival_df(region_df,region,cut_off=0.5):

    df = pd.DataFrame(index = region_df.index)
    df['time_max'] = region_df[['days_to_death','days_to_last_follow_up']].fillna(0).astype(float).max(axis=1)   
    df['time_sum'] = region_df[['days_to_death','days_to_last_follow_up']].fillna(0).astype(float).sum(axis=1)   
#     df = df[df['time_max']==df['time_sum']]
    df['time']=df['time_max']
    df = df[df['time']>0]
    
    # ==== get the clinical info, using top/bottom 50/25% as cutoff
    # ==== keep at least 5 patients
    a_index = df.index[:max(int(len(df.index)*cut_off),0)]
    b_index = df.index[-1*max(int(len(df.index)*cut_off),0):]
    df = df.loc[a_index.union(b_index)]
    # print('treat:\t',len(a_index),'\nctrl:\t',len(b_index))
    df.loc[a_index,'group']='treat'   
    df.loc[b_index,'group']='ctrl'  
    df['status']= region_df.loc[df.index]['vital_status'] 
    df['score']= region_df.loc[df.index][region] 
    df.loc[df['status']=='Dead','status']=1
    df.loc[df['status']=='Alive','status']=0
    return df

    

    
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
    if pvalue<0.01:
        # survival curves
        plt.figure(figsize=(2.6,2.6))
        ax = plt.subplot(111)
        
        kmf_control = KaplanMeierFitter()
        #g1 = kmf_control.fit(t1, e1, label=legends[0]).plot(ax=ax,show_censors=True,\
        g1 = kmf_control.fit(t1, e1).plot(ax=ax,show_censors=True,\
                                          label='more open',\
                            censor_styles={'ms': 12, 'marker': '+'},ci_show=False,color='red',ls='-')
        
        kmf_exp = KaplanMeierFitter()
        #g2 = kmf_exp.fit(t2, e2, label=legends[1]).plot(ax=ax,show_censors=True,\
        g2 = kmf_exp.fit(t2, e2).plot(ax=ax,show_censors=True,\
                                      label='less open',\
                        censor_styles={'ms': 12, 'marker': '+'},ci_show=False,color='k',ls='--')
        
        handles, labels = ax.get_legend_handles_labels();print(labels)
        lg = ax.legend(handles, legends,loc='lower left',borderaxespad=0,handletextpad=.2,labelspacing=.2,
                       handlelength=1,frameon=False)
        if pvalue<1:
              plt.axes().text(df['time'].max()*0.45,0.45,'p={:.1e}'.format(pvalue),fontsize=16,ha='center')
        plt.ylim([0.0,1.0])
    #     plt.xlim([0,max_val*1])
        plt.title('{}'.format(title),fontsize=16)
        plt.xlabel('Days',fontsize=16)
        # plt.ylabel('Survival probability \n of {} patients'.format(df.shape[0]),fontsize=18)
        plt.ylabel('Survival probability',fontsize=16)
        plt.savefig(figname,bbox_inches='tight',pad_inches=.1,dpi=600,transparent=True)
        # plt.show()
        plt.close()
    return results





    
# ==== main  
indir = 'f1_ATAC_overlap_SE_caseID'
outdir = 'f2_caseID_each_SE_vs_clinical'
os.makedirs(outdir,exist_ok=True)
os.makedirs(outdir+os.sep+'fig',exist_ok=True)
name_match = pd.read_excel('../../data/TCGA/TCGA-ATAC_SE_cancerType_match.xlsx',index_col=0)   
name_match = name_match.dropna()
filtered_df = pd.read_excel('../../data/TCGA/TCGA-ATAC_clustered_samples.xlsx',index_col=0)   


for cancertype in ['BRCA','COAD']:
    cancertype_SE = name_match.loc[cancertype].SE
    cancertype_SE_rename = name_match.loc[cancertype].SE_rename
    filtered_id = filtered_df[filtered_df.cohort==cancertype].case_id
    # ==== read TCGA peaks overlap with SE data
    overlap_file = '{}/{}_ATAC_overlap_{}_SE_caseID.bed'.format(indir,cancertype,cancertype_SE_rename)
    with open(overlap_file) as sig_inf:
        sig_df = pd.read_csv(sig_inf,sep='\t',index_col=3)
    sig_df = sig_df.iloc[:,6:]
    # patient ID by TCGA data
    sig_case_id = [i.split('_')[1] for i in sig_df.columns]
    sig_df.columns = sig_case_id
    # == remove non-basal BRCA
    sig_case_id = [i for i in sig_case_id if i in filtered_id.values]
    
    
    # ==== read clinical data
    clinical_data = '2022-03-20'  
    clinical_file = '../../data/TCGA//clinical.project-TCGA-{}.{}.json'.format(cancertype,clinical_data)
    with open(clinical_file) as clinical_inf:
        clinical_list = json.load(clinical_inf)
    # patient ID by clinical data
    clinical_case_id=[i['case_id'] for i in clinical_list]


    # ==== read clinical info for each case ID
    clinical_df = pd.DataFrame()
    for clinical_case in clinical_list:
        case_id = clinical_case['case_id']
        if case_id in sig_case_id:
            try:
                # clinical_df.loc[case_id,'age_at_diagnosis'] = clinical_case['diagnoses'][0]['age_at_diagnosis']
                clinical_df.loc[case_id,'vital_status'] = clinical_case['demographic']['vital_status'] 
                if clinical_df.loc[case_id,'vital_status']=='Alive':
                    clinical_df.loc[case_id,'days_to_last_follow_up'] = clinical_case['diagnoses'][0]['days_to_last_follow_up']
                elif clinical_df.loc[case_id,'vital_status']=='Dead':
                    clinical_df.loc[case_id,'days_to_death'] = clinical_case['demographic']['days_to_death']
            except:
                print('\n\n=======\n',case_id)
    clinical_df.to_csv(outdir+os.sep+'{}_clinical_info.csv'.format(cancertype))


    # ==== for each ATAC-seq peak, check openness cor clinical results
    region_outdf = pd.DataFrame()
    for region in sig_df.index[:]:
        sig_score = sig_df.loc[region][clinical_df.index]
        region_df = pd.concat([clinical_df,sig_score],axis=1)
        region_df = region_df.sort_values(by=[region],ascending=False)
        df = return_survival_df(region_df,region,0.5)  
        figname=outdir+os.sep+'/fig/{}_{}_survival_50th.pdf'.format(cancertype,region)
        legends = ['more open','less open']
        results = survival_for_two(df,'treat','ctrl',legends,region,figname)
        if results.p_value<0.01:
            df.to_csv(outdir+os.sep+'fig/{}_{}_survival_50th.csv'.format(cancertype,region))
        
        # save info
        ix = df['group'] == 'treat'
        region_outdf.loc[region,'treat more open avg atac-score'] = df.loc[ix,'score'].values.mean().round(2)
        region_outdf.loc[region,'ctrl less open avg atac-score'] = df.loc[~ix,'score'].values.mean().round(2)
        region_outdf.loc[region,'treat time'] = df.loc[ix,'time'].values.mean().round(2)
        region_outdf.loc[region,'ctrl time'] = df.loc[~ix,'time'].values.mean().round(2)
        # region_outdf.loc[region,'log rank s'] = results.test_statistic  
        region_outdf.loc[region,'log rank p'] = results.p_value
    region_outdf['fdr'] = fdr_adj_p(region_outdf['log rank p'],'log rank p')
    region_outdf = region_outdf.sort_values(by='log rank p',ascending=True)
    region_outdf.to_csv('{}/{}_logrank_info.csv'.format(outdir,cancertype))

         

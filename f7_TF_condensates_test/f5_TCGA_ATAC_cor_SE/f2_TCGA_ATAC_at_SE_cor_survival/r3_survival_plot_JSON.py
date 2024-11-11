import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *

import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=18
import seaborn as sns
sns.set(font_scale=1.2)

sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style('ticks')
import re,bisect
from lifelines.statistics import logrank_test
from lifelines import KaplanMeierFitter
matplotlib.rcParams["font.sans-serif"] = ["Arial"]



def survival_for_two(df,treat,ctrl,legends,title,figname):
    
    # select the time and status info for treat and control group
    ix = df['group'] == treat
    t1 = df.loc[ix]['time'];print(t1.shape)
    e1 = df.loc[ix]['status'] 
    t2 = df.loc[~ix]['time'];print(t2.shape)
    e2 = df.loc[~ix]['status']
    
    results = logrank_test(t1,t2,event_observed_A = e1,event_observed_B = e2)
    pvalue = results.p_value;print('pvalue:\t{}'.format(pvalue))
    
    # survival curves
    plt.figure(figsize=(2.6,2.6))
    ax = plt.subplot(111)
    
    kmf_control = KaplanMeierFitter()
    #g1 = kmf_control.fit(t1, e1, label=legends[0]).plot(ax=ax,show_censors=True,\
    g1 = kmf_control.fit(t1, e1).plot(ax=ax,show_censors=True,\
                                      label='more open',\
                        censor_styles={'ms': 12, 'marker': '+'},ci_show=False,c='red',ls='-')
    
    kmf_exp = KaplanMeierFitter()
    #g2 = kmf_exp.fit(t2, e2, label=legends[1]).plot(ax=ax,show_censors=True,\
    g2 = kmf_exp.fit(t2, e2).plot(ax=ax,show_censors=True,\
                                  label='less open',\
                    censor_styles={'ms': 12, 'marker': '+'},ci_show=False,c='k',ls='--')
    
    handles, labels = ax.get_legend_handles_labels();print(labels)
    lg = ax.legend(handles, legends,loc='lower left',borderaxespad=0,handletextpad=.2,labelspacing=.2,
                   handlelength=1,frameon=False)
    if pvalue<1:
          plt.axes().text(df['time'].max()*0.45,0.45,'p={:.1e}'.format(pvalue),fontsize=16,ha='center')
    plt.ylim([0.0,1.0])
#     plt.xlim([0,max_val*1])
    plt.title('{}'.format(title),fontsize=18)
    plt.xlabel('Days',fontsize=18)
    # plt.ylabel('Survival probability \n of {} patients'.format(df.shape[0]),fontsize=18)
    plt.ylabel('Survival probability',fontsize=18)
    plt.savefig(figname,bbox_inches='tight',pad_inches=.1,dpi=600,transparent=True)
    plt.show()
    plt.close()
    return results




def return_survival_df(clinical_df,cut_off=0.5):

    df = pd.DataFrame(index = clinical_df.index)
    df['time_max'] = clinical_df[['days_to_death','days_to_last_follow_up']].fillna(0).astype(float).max(axis=1)   
    df['time_sum'] = clinical_df[['days_to_death','days_to_last_follow_up']].fillna(0).astype(float).sum(axis=1)   
#     df = df[df['time_max']==df['time_sum']]
    df['time']=df['time_max']
    df = df[df['time']>0]
    
    # ==== get the clinical info, using top/bottom 50/25% as cutoff
    # ==== keep at least 5 patients
    a_index = df.index[:max(int(len(df.index)*cut_off),0)]
    b_index = df.index[-1*max(int(len(df.index)*cut_off),0):]
    df = df.loc[a_index.union(b_index)]
    print('treat:\t',len(a_index),'\nctrl:\t',len(b_index))
    df.loc[a_index,'group']='treat'   
    df.loc[b_index,'group']='ctrl'  
    df['status']= clinical_df.loc[df.index]['vital_status'] 
    df.loc[df['status']=='Dead','status']=1
    df.loc[df['status']=='Alive','status']=0
    return df


# ==== main 
indir = 'f2_caseID_ranksum_vs_clinical'
outdir = 'f3_survival_figs'
os.makedirs(outdir,exist_ok=True)

name_match = pd.read_excel('../f1_TCGA_ATAC_at_SE_rep/TCGA_ATAC_SE_cancerType_match.xlsx',index_col=0)   
name_match = name_match.dropna()
  
for cancertype in name_match.index[:]:
    if cancertype in ['CESC','PRAD']:
        continue
    cancertype_SE = name_match.loc[cancertype].SE
    cancertype_SE_rename = name_match.loc[cancertype].SE_rename
    clinical_file = indir+os.sep+'{}_irwin_hall_with_clinical_info.csv'.format(cancertype)
    clinical_df = pd.read_csv(clinical_file,sep=',',index_col=0) 
    treat='treat'
    ctrl='ctrl'
    print('\n\n=====\n',cancertype)

    ### separate patient into two groups evenly, rank by mean different score
    # ==== cutoff: top/bottom 50%
    df = return_survival_df(clinical_df,0.5)
    df.to_csv(outdir+os.sep+'{}_survival_50th.csv'.format(cancertype))
    title = '{} \n  {} patients in total'.format(cancertype,df.shape[0])
    figname=outdir+os.sep+'{}_survival_50th.pdf'.format(cancertype)
    legends = ['top 50%','bottom 50%']
    survival_for_two(df,treat,ctrl,legends,title,figname)


    # ==== cutoff: top/bottom 25%
    df = return_survival_df(clinical_df,.25)
    df.to_csv(outdir+os.sep+'{}_survival_25th.csv'.format(cancertype))
    # title = '{} \n total {} patients'.format(cancertype,df.shape[0])
    figname=outdir+os.sep+'{}_survival_25th.pdf'.format(cancertype)
    legends = ['top 25%','bottom 25%']
    survival_for_two(df,treat,ctrl,legends,title,figname)




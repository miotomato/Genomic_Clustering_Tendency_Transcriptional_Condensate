import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import json
#from GenomeData import *
import re,bisect

import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=16
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style('ticks')

from lifelines.statistics import logrank_test
from lifelines import KaplanMeierFitter

    

    
def survival_for_two(df,treat,ctrl,figname,cancertype,labmda_cutoff):
    
    # select the time and status info for treat and control group
    ix = df['group'] == treat
    t1 = df.loc[ix]['time']#;print(t1.shape)
    e1 = df.loc[ix]['status'] 
    t2 = df.loc[~ix]['time']#;print(t2.shape)
    e2 = df.loc[~ix]['status']
    
    results = logrank_test(t1,t2,event_observed_A = e1,event_observed_B = e2)
    pvalue = results.p_value#;print('pvalue:\t{}'.format(pvalue))
    
    if 1:
    # if pvalue<0.01:
        # survival curves
        plt.figure(figsize=(2.6,2.6))
        ax = plt.subplot(111)
        
        kmf_control = KaplanMeierFitter()
        g1 = kmf_control.fit(t2, e2).plot(ax=ax,show_censors=True,\
                                          label='more open',\
                            censor_styles={'ms': 12, 'marker': '+'},ci_show=False,color='red',ls='-')
        
        kmf_exp = KaplanMeierFitter()
        g2 = kmf_exp.fit(t1, e1).plot(ax=ax,show_censors=True,\
                                      label='less open',\
                        censor_styles={'ms': 12, 'marker': '+'},ci_show=False,color='k',ls='--')
        
        handles, labels = ax.get_legend_handles_labels();print(labels)
        lg = ax.legend(handles, labels,loc='lower left',borderaxespad=0,handletextpad=.2,labelspacing=.2,
                       handlelength=1,frameon=False)
        if pvalue<1:
              plt.axes().text(df['time'].max()*0.45,0.45,'p={:.1e}'.format(pvalue),fontsize=16,ha='center')
        plt.ylim([0.0,1.0])
    #     plt.xlim([0,max_val*1])
        plt.title('{}'.format(cancertype),fontsize=16)
        plt.xlabel('Days',fontsize=16)
        # plt.ylabel('Survival probability \n of {} patients'.format(df.shape[0]),fontsize=18)
        plt.ylabel('Survival probability',fontsize=16)
        plt.savefig(figname,bbox_inches='tight',pad_inches=.1,dpi=600,transparent=True)
        plt.show()
        plt.close()
    return results





    
# ==== main  
outdir = 'f3_survival_figs'
os.makedirs(outdir,exist_ok=True)
name_match = pd.read_excel('../../f1_TCGA_ATAC_at_SE_rep/TCGA_ATAC_SE_cancerType_match.xlsx',index_col=0)   
name_match = name_match.dropna()
  

# for cancertype in name_match.index[:1]:
for cancertype in ['BRCA','COAD']:
    atac_data = pd.read_csv('f1_data/{}_atacseq_sig.csv'.format(cancertype),index_col=0)
    clinical_data = pd.read_csv('f1_data/{}_clinical.csv'.format(cancertype),index_col=0)

    # save the logrank info
    summary_score=pd.DataFrame()
    for labmda_cutoff in [0.01,0.03,0.05,0.07,0.09,0.11,0.13][:]:
        coef = pd.read_csv('f2_figs/{}_coef_{}.csv'.format(cancertype,labmda_cutoff),header=None,index_col=0)
        coef = coef[coef!=0].dropna()
        coef.columns = ['coef']
        # coef.to_csv('{}/{}_coef.csv'.format(outdir,labmda_cutoff))
        
        inner_product = np.inner(atac_data[coef.index],coef['coef'])
        inner_df = pd.DataFrame(inner_product)
        inner_df.index = atac_data.index
        c1_lower = inner_df.sort_values(by=[0],ascending=True).index[:int(.5*inner_df.shape[0])]
        c2_higher = inner_df.sort_values(by=[0],ascending=True).index[-1*int(.5*inner_df.shape[0]):]
        # plot the clinical figs
        survival_df = pd.DataFrame(index = np.append(c1_lower,c2_higher))#;print(survival_df)
        survival_df.loc[set(c1_lower).intersection(clinical_data.index),'group']='c1'
        survival_df.loc[set(c2_higher).intersection(clinical_data.index),'group']='c2'
        survival_df['status']= clinical_data.loc[survival_df.index]['status']         
        survival_df['time']= clinical_data.loc[survival_df.index]['time']
        survival_df = survival_df.dropna(axis=0,how='any')
        avg_time_c1 = survival_df[survival_df['group']=='c1']['time'].mean()
        avg_time_c2 = survival_df[survival_df['group']=='c2']['time'].mean()
        figname='{}/{}_{}_clinical.pdf'.format(outdir,cancertype,labmda_cutoff)
        results = survival_for_two(survival_df,'c1','c2',figname,cancertype,labmda_cutoff)
        survival_df.to_csv('{}/{}_{}_clinical_info.csv'.format(outdir,cancertype,labmda_cutoff))
        # data for clinical analysis
        summary_score.loc[labmda_cutoff,'logrank_pvalue'] = results.p_value
        summary_score.loc[labmda_cutoff,'logrank_statistics'] = results.test_statistic
        summary_score.loc[labmda_cutoff,'logrank_c1_low_days'] = avg_time_c1
        summary_score.loc[labmda_cutoff,'logrank_c2_high_days'] = avg_time_c2
    
    summary_score.to_csv('{}/{}_summary_score.csv'.format(outdir,cancertype))        
        



       

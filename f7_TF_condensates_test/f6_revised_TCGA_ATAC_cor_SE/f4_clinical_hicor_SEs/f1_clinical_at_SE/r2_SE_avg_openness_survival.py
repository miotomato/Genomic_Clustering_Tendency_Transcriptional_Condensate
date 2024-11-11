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
    if pvalue<0.05 or re.search('SE',title):
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
        plt.title('{}'.format(title),fontsize=12)
        plt.xlabel('Days',fontsize=14)
        # plt.ylabel('Survival probability \n of {} patients'.format(df.shape[0]),fontsize=18)
        plt.ylabel('Survival probability',fontsize=14)
        plt.savefig(figname,bbox_inches='tight',pad_inches=.1,dpi=600,transparent=True)
        # plt.show()
        plt.close()
    return results





    
# ==== main  
outdir = 'f2_SE_avg_openness_survival'
os.makedirs(outdir,exist_ok=True)
os.makedirs(outdir+os.sep+'fig',exist_ok=True)


for cancertype in ['BRCA','COAD']:
    df = pd.read_csv('f1_SE_openness_add_clinical/{}_SE_openness_add_clinical.csv'.format(cancertype),index_col=0)
    se_cols = df.columns[5:]
    patients_each_group = int(len(df.index)*0.5)
    print(df.shape)
    
    # ==== survival analysis for each SE
    region_outdf = pd.DataFrame()
    for region in se_cols:
        df_tmp = df[['time','status',region]]
        df_tmp = df_tmp.sort_values(by=region,ascending=False)
        
        # separate the patients into two groups
        treat_index = df_tmp.index[:patients_each_group]
        ctrl_index = df_tmp.index[-1*patients_each_group:]
        df_tmp = df_tmp.loc[treat_index.union(ctrl_index)]
        df_tmp.loc[treat_index,'group']='treat'   
        df_tmp.loc[ctrl_index,'group']='ctrl'  
        # ==== plot the figs
        figname=outdir+os.sep+'/fig/{}_{}_survival.pdf'.format(cancertype,region)
        legends = ['more open','less open']
        results = survival_for_two(df_tmp,'treat','ctrl',legends,region,figname)
        if results.p_value<0.05 or re.search('SE',region):
            df_tmp.to_csv(outdir+os.sep+'fig/{}_{}_survival.csv'.format(cancertype,region))

        # save info
        ix = df_tmp['group'] == 'treat'
        region_outdf.loc[region,'treat more open avg atac-score'] = df_tmp.loc[ix,region].values.mean().round(2)
        region_outdf.loc[region,'ctrl less open avg atac-score'] = df_tmp.loc[~ix,region].values.mean().round(2)
        region_outdf.loc[region,'treat time'] = df_tmp.loc[ix,'time'].values.mean().round(2)
        region_outdf.loc[region,'ctrl time'] = df_tmp.loc[~ix,'time'].values.mean().round(2)
        # region_outdf.loc[region,'log rank s'] = results.test_statistic  
        region_outdf.loc[region,'log rank p'] = results.p_value
    region_outdf['fdr'] = fdr_adj_p(region_outdf['log rank p'],'log rank p')
    region_outdf = region_outdf.sort_values(by='log rank p',ascending=True)
    region_outdf.to_csv('{}/{}_logrank_info.csv'.format(outdir,cancertype))

         

import os,sys,argparse
# import fileinput,time
import glob
import re,bisect
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=9
import seaborn as sns
sns.set(font_scale=.9)
sns.set_style("whitegrid", {'axes.grid' : False})
import scipy
import scipy.optimize
sns.set_style("ticks")
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams["mathtext.rm"] = "Arial"
from scipy.interpolate import interpn
from scipy.stats import gaussian_kde
from matplotlib import gridspec


def box_plot(pos_len,vals,mark_p=False):
    positions = np.arange(pos_len)
    box_step = int(len(vals)/len(positions))
    box_vals = []
    for ii in positions:
        box_val = vals[ii*box_step:(ii+1)*box_step]
        box_vals.append(box_val)
        s,p = scipy.stats.ttest_1samp(box_val,0)#;print(s,p)
        p_val = '*' if p<.05 else ''
        if mark_p == True:
            ax.text(ii, np.mean(box_val), p_val, fontsize=27,c='red',ha='center',va='center')
    g = plt.boxplot(box_vals,positions=positions,widths = .6,patch_artist=True,\
                boxprops=dict(color='k',facecolor='w',fill=None,lw=1),\
                medianprops=dict(color='grey'),showfliers=False)  

    


# == check the TF motifs
outdir = 'f2_TFMS_CP_cor_SE_heatmap'
os.makedirs(outdir,exist_ok=True)
test_types = ['KS_statistics','Rank_Sum_statistics','T_test_statistics']


for test_type in test_types[:]:
    df = pd.read_csv('data/TFMS_CP_SE_enrich.csv',index_col=0)
    if test_type == 'KS_statistics':
        df = df.sort_values(by='log10-dis ks_2samp-s signed',ascending=False)
    elif test_type == 'Rank_Sum_statistics':
        df = df.sort_values(by='log10-dis Wilcoxon-rank-sum-s',ascending=False)
    else:
        df = df.sort_values(by='log10-dis T-test-s',ascending=False)
    df = df.iloc[:,:]
    value_dic = {'Motif_num':df['#TFMS'],
                 'Motif_length':df['len-of-TFMS'],
                 'T_test_statistics':df['log10-dis T-test-s'],
                 'Rank_Sum_statistics':df['log10-dis Wilcoxon-rank-sum-s'],
                 'KS_statistics':df['log10-dis ks_2samp-s signed'],
                 'P_LT_005_num':df['#p<0.05'],
                 'P_LT_005_percentage':df['#p<0.05']/df['#TFMS'],
                 'P_LT_005_log2_OR_at_SE':np.log2(df['TFMS-p<0.05 enrich-at-SE-fisher-exact-s']),
                 }
    
    
    plt.figure(figsize=(7,3))
    width_ratios = [9,.1]
    height_ratios = [.0,.1,.3,1,1,1]
    # gs = gridspec.GridSpec(1,3,width_ratios = width_ratios, hspace=.15)
    gs = gridspec.GridSpec(6,2,width_ratios = width_ratios,
                           height_ratios = height_ratios, 
                           wspace=.1,hspace=.1)


    # color bar
    cmap = plt.cm.bwr
    vmax=1
    vmin=-1
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    
    ax = plt.subplot(gs[1:3,1]) 
    cbar = matplotlib.colorbar.ColorbarBase(ax,cmap=cmap,norm=norm,orientation='vertical')
    cbar.set_ticks([vmin,0,vmax])
    cbar.set_label('log2 OR',labelpad=-27,
                   fontsize=9,va='center')
    ax.tick_params(axis='y',direction='out', length=2, width=1, colors='black')   
    
    # ==== heatmap of SE enrichment
    ax = plt.subplot(gs[2,0])
    vals = value_dic['P_LT_005_log2_OR_at_SE']
    vals = np.transpose(vals.to_frame())
    g = sns.heatmap(vals,cmap=cmap,cbar=False \
                    ,vmax=vmax,vmin=vmin\
                    ,xticklabels=False,yticklabels=False\
                    ,ax=ax)
    
    # ==== box of SE enrichment
    ax = plt.subplot(gs[3,0])
    vals = value_dic['P_LT_005_log2_OR_at_SE']
    box_plot(10,vals,mark_p=True)
    ax.axhline(y=0,color='k',lw=1.2,ls='--')
    ax.set_xticks([])
    ax.set_ylabel('log2 OR\n',ha='center')
    
    
    # ==== box of SE enrichment
    # ax = plt.subplot(gs[2,0])
    # vals = value_dic['P_LT_005_percentage']
    # box_plot(10,vals,mark_p=False)
    # # ax.axhline(y=0,color='k',lw=1.2,ls='--')
    # ax.set_xticks([])
    # ax.set_ylabel('log2 OR')
    
    # ==== % of TFMS with p<0.05
    ax = plt.subplot(gs[4,0])
    vals = 100*value_dic['P_LT_005_percentage']
    # g = ax.scatter(np.arange(len(vals)),vals,s=3,c='k')
    g = ax.bar(np.arange(len(vals)),vals,lw=0,width=1,color='k')
    ax.axhline(y=5,color='r',lw=1.2,ls='--')
    ax.set_xlim([-.5,len(vals)-.5])
    ax.set_ylim([0,55])
    ax.set_xticks([])
    ax.set_ylabel('% p<0.05')
    
    
    # ==== plot the CP
    ax = plt.subplot(gs[5,0])
    vals = value_dic[test_type]
    g = ax.scatter(np.arange(len(vals)),vals,s=3,c='k')
    ax.axhline(y=0,color='k',lw=1.2,ls='--')
    ax.set_xlim([-.5,len(vals)-.5])
    ax.set_ylabel('TFMS CP ',ha='center')
    
    
    # plt.xlabel('Rank')
    plt.xlabel('TF rank')
    plt.savefig('{}/CP_by_{}.pdf'.format(outdir,test_type),bbox_inches='tight',pad_inches=0.02,transparent=True)
    plt.show()
    plt.close()

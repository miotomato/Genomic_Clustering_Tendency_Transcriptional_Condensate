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


def box_plot(pos_len,vals,or_cutoff,mark_p=False):
    positions = np.arange(pos_len)
    box_step = int(len(vals)/len(positions))
    box_vals = []
    for ii in positions:
        box_val = vals[ii*box_step:(ii+1)*box_step]
        # box_val = np.abs(box_val)
        box_vals.append(box_val)
        s,p = scipy.stats.ttest_1samp(box_val,np.log2(or_cutoff),alternative='greater');print(len(box_val),s,p)
        p_val = '*' if p<.05 else ''
        if mark_p == True:
            # ax.text(ii, np.mean(box_val), p_val, fontsize=27,c='red',ha='center',va='center')
            ax.text(ii, np.percentile(vals,96), p_val, fontsize=27,c='red',ha='center',va='center')
    g = plt.boxplot(box_vals,positions=positions,widths = .6,patch_artist=True,\
                boxprops=dict(color='k',facecolor='w',fill=None,lw=1),\
                medianprops=dict(color='grey'),showfliers=False)  

    


# == check the TF motifs
outdir = 'f3_TFMS_CP_cor_SE_heatmap'
os.makedirs(outdir,exist_ok=True)

df1 = pd.read_csv('../f1_TF_cluster_potential/f2_cor_CP_SE_AICAP/data//TFMS_CP_SE_enrich.csv',index_col=0)
df1 = df1[['enrich-at-SE-fisher-exact-s']]
df2 = pd.read_csv('f2_gamma_alpha_vs_KS/data_gamma_fit.csv',index_col=0)
df = pd.concat([df1,df2],axis=1,join='inner')
df = df.sort_values(by='log10-dis ks_2samp-s signed',ascending=False)
df.to_csv('{}/TFMS_CP_SeEnrich_GammaFit.csv'.format(outdir))

test_types = ['TFMS_CP','Gamma_alpha','function_alpha']

for test_type in test_types[:]:
    
    if test_type == 'TFMS_CP':
        df = df.sort_values(by='log10-dis ks_2samp-s signed',ascending=False)
    elif test_type == 'Gamma_alpha':
        df = df.sort_values(by='alpha',ascending=True)
    elif test_type == 'Function_alpha':
        df = df.sort_values(by='f(alpha)',ascending=False)
    df = df.iloc[:,:]
    value_dic = {
        # 'Motif_num':df['#TFMS'],
        # 'Motif_length':df['len-of-TFMS'],
        'TFMS_CP':df['log10-dis ks_2samp-s signed'],
        'Gamma_alpha':df['alpha'],
        'function_alpha':df['f(alpha)'],
        'log2_OR_at_SE':np.log2(df['enrich-at-SE-fisher-exact-s']),
        # 'P_LT_005_num':df['#p<0.05'],
        # 'P_LT_005_percentage':df['#p<0.05']/df['#TFMS'],
        # 'P_LT_005_log2_OR_at_SE':np.log2(df['TFMS-p<0.05 enrich-at-SE-fisher-exact-s']),
        # 'P_LT_001_num':df['#p<0.01'],
        # 'P_LT_001_percentage':df['#p<0.01']/df['#TFMS'],
        # 'P_LT_001_log2_OR_at_SE':np.log2(df['TFMS-p<0.01 enrich-at-SE-fisher-exact-s']),
        }
    
    
    plt.figure(figsize=(6,4))
    width_ratios = [9,.2]
    height_ratios = [1,1,1]
    # gs = gridspec.GridSpec(1,3,width_ratios = width_ratios, hspace=.15)
    gs = gridspec.GridSpec(3,2,width_ratios = width_ratios,
                           height_ratios = height_ratios, 
                           wspace=.1,hspace=.1)



    # ==== plot the CP
    ax = plt.subplot(gs[0,0])
    vals = value_dic[test_type]
    g = ax.scatter(np.arange(len(vals)),vals,s=3,c='k')
    ax.axhline(y=0,color='k',lw=1.2,ls='--')
    ax.set_xlim([-.5,len(vals)-.5])
    ax.xaxis.tick_top()
    ax.set_ylabel(test_type,ha='center')
    if test_type == 'Gamma_alpha':
        ax.set_ylabel('Gamma $k$',ha='center')
    # ax.set_title(test_type,va='bottom')


    # color bar
    cmap = plt.cm.bwr
    vmax= 1
    vmin=-1*vmax
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    
    ax = plt.subplot(gs[1,1]) 
    cbar = matplotlib.colorbar.ColorbarBase(ax,cmap=cmap,norm=norm,orientation='vertical')
    cbar.set_ticks([vmin,0,vmax])
    cbar.set_label('log2 OR',labelpad=-30,
                   fontsize=9,va='center')
    ax.tick_params(axis='y',direction='out', length=2, width=1, colors='black')   
    
    # ==== heatmap of SE enrichment
    ax = plt.subplot(gs[1,0])
    vals = value_dic['log2_OR_at_SE']
    vals = np.transpose(vals.to_frame())
    g = sns.heatmap(vals,cmap=cmap,cbar=False \
                    ,vmax=vmax,vmin=vmin\
                    ,xticklabels=False,yticklabels=False\
                    ,ax=ax)
    ax.set_ylabel('Enrichment at \n SE (log2 OR)\n',ha='center')


    # ==== box of SE enrichment
    ax = plt.subplot(gs[2,0])
    vals = value_dic['log2_OR_at_SE']
    or_cutoff = 1.2
    box_plot(20,vals,or_cutoff,mark_p=True)
    ax.axhline(y=np.log2(or_cutoff),color='k',lw=1.2,ls='--')
    ax.set_xticks([])
    ax.set_ylabel('log2 OR',ha='center')
        
    
    # plt.xlabel('Rank')
    # plt.xlabel('TF rank')
    # plt.title(test_type)
    plt.savefig('{}/CP_by_{}.pdf'.format(outdir,test_type),bbox_inches='tight',pad_inches=0.02,transparent=True)
    plt.show()
    plt.close()

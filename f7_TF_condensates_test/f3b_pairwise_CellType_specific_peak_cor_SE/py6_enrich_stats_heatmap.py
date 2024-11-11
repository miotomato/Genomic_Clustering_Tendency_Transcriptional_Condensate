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
matplotlib.rcParams['font.size']=12
import seaborn as sns
sns.set(font_scale=1.1)
sns.set_style("whitegrid", {'axes.grid' : False})
import scipy
import scipy.optimize
sns.set_style("ticks")
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
from matplotlib import gridspec


# == check the TF motifs
indir = 'f5_cellType_specific_peaks_cor_SE_over_AllPeaks'
outdir = 'f6_stats_heatmap'
os.makedirs(outdir,exist_ok=True)

list_df = pd.read_csv('f2_selected_data/pairwise_sharedTF_count_id.csv',index_col=0)
for celltype_pairs in list_df.index[:]:

    df_s = pd.read_csv('{}/{}/fisher_s.csv'.format(indir,celltype_pairs),index_col=0)
    df_p = pd.read_csv('{}/{}/fisher_p.csv'.format(indir,celltype_pairs),index_col=0)
    df_p = df_p.replace(0,1e-299)
    
    kept_index = df_p[(df_p<0.001).sum(axis=1)>1].index
    df_s = df_s.loc[kept_index]
    df_p = df_p.loc[kept_index]
    
    # log2 of odds ratio
    df_s = np.log2(df_s)
    # df_s = df_s.fillna(0)
    # rank col and row by sum
    ranked_row = df_s.sum(axis=1).sort_values(ascending=True).index
    df_s = df_s.loc[ranked_row]
    # ranked_col = df_s.sum(axis=0).sort_values(ascending=False).index
    # df_s = df_s[ranked_col]
    
    
    # == plot the heatmap
    plt.figure(figsize=(3,3.3+np.sqrt(df_s.shape[0])))
    width_ratios = [3,.4,.1,.4]
    gs = gridspec.GridSpec(4,4,width_ratios=width_ratios,wspace=.3,hspace=.5)
    
    # heatmap
    pal = plt.cm.bwr
    vmin=-2
    vmax=2
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    color_map = matplotlib.cm.ScalarMappable(norm=norm, cmap=pal)
    scatter_scale = .5
    ax = plt.subplot(gs[:,0])   
    for ii in np.arange(len(df_s.index)):
        for jj in np.arange(len(df_s.columns)):
            odds_ratio = df_s.loc[df_s.index[ii],df_s.columns[jj]]
            pvalue = df_p.loc[df_s.index[ii],df_s.columns[jj]]
            if not np.isnan(odds_ratio):
                color = color_map.to_rgba(odds_ratio)
                size = -1*np.log10(pvalue)
                ax.scatter(jj,ii,s = scatter_scale*size, color = color)
                # ax.text(jj,ii,jj,fontsize=12)
            else:
                ax.scatter(jj,ii,marker = 'x',s = scatter_scale*30, color = 'k')
    ax.set_yticks(np.arange(len(df_s.index)))
    ax.set_yticklabels(df_s.index,fontsize=13)
    ax.set_xticks(np.arange(len(df_s.columns)))
    ax.set_xticklabels(df_s.columns,rotation=60,fontsize=13)
    ax.set_xlim([-.5,len(df_s.columns)-.5])
    ax.set_ylim([-.5,len(df_s.index)-.5])
    ax.set_title(celltype_pairs,fontsize=13)
    
    # color bar
    ax = plt.subplot(gs[0,3]) 
    cbar = matplotlib.colorbar.ColorbarBase(ax,cmap=pal,norm=norm,orientation='vertical')
    cbar.set_ticks([vmin,0,vmax])
    cbar.set_label('log2 (odds ratio)',labelpad=-50,fontsize=13)
    # pvalue legend
    ax = plt.subplot(gs[1,3])   
    ax.scatter(0,.15,s = -1*scatter_scale*np.log10(1e-10), color = 'k')
    ax.scatter(0,.5,s = -1*scatter_scale*np.log10(1e-50), color = 'k')
    ax.scatter(0,.9,s = -1*scatter_scale*np.log10(1e-100), color = 'k')
    ax.text(0,.0,'1e-10',ha='center')
    ax.text(0,.3,'1e-50',ha='center')
    ax.text(0,.67,'1e-100',ha='center')
    ax.set_xlim([-1,1])
    ax.set_ylim([.0,1.05])
    ax.set_title('P-value',ha='center',fontsize=13)
    ax.axis('off')
    
    
    # heatmap of AICAP
    # read TF AICAP
    idrfile = '../f2_TF_motif_IDR_SE_mutation/data/13059_2021_2456_MOESM2_ESM.xlsx'
    idrdf = pd.read_excel(idrfile,sheet_name='condition 2(1,6-HD-2)',index_col=0)
    vmin=-2;vmax=2
    ax = plt.subplot(gs[:,1]) 
    plot_data = idrdf.loc[df_s.index[::-1],['AICAP']]
    plot_data = np.log2(plot_data)
    g = sns.heatmap(plot_data,cmap=plt.cm.PiYG,cbar=False\
                    ,vmax=vmax,vmin=vmin\
                    ,ax=ax\
                    ,xticklabels=True,yticklabels=False)
    ax.set_xticklabels(['AICAP'],rotation=60)

    # color bar
    ax = plt.subplot(gs[2,3]) 
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    cbar = matplotlib.colorbar.ColorbarBase(ax,cmap=plt.cm.PiYG,norm=norm,orientation='vertical')
    cbar.set_ticks([-2,0,2])
    cbar.set_label('log2 AICAP',labelpad=-50,fontsize=13)
    
    
    plt.savefig('{}/{}_heatmap.pdf'.format(outdir,celltype_pairs),bbox_inches='tight',pad_inches=0.05,transparent=True)
    plt.show()
    plt.close()
    







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
matplotlib.rcParams['font.size']=14
import seaborn as sns
sns.set(font_scale=1.)
sns.set_style("whitegrid", {'axes.grid' : False})
import scipy
import scipy.optimize
sns.set_style("ticks")
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
from matplotlib import gridspec



outdir = 'f1_CP_heatmap'
os.makedirs(outdir+os.sep+'_csv',exist_ok=True)

indirs = ['CP_TFBS_vs_TFMS','CP_TFBS_with_motif_vs_TFMS','CP_TFBS_without_motif_vs_TFMS']
group_info = {
'CP': {'data_median_fisher':
       ['log10-distance CP t-stats',
        'log10-distance CP pvalue Fisher-TF-Celltype-combined'],
       'data_fisher_fisher':
       ['log10-distance CP t-stats',
        'log10-distance CP pvalue Fisher-data-combined Fisher-TF-Celltype-combined'],
    },
'SE': {'data_median_fisher':
       ['fisher_exact_s',
        'fisher_exact_p Fisher-TF-Celltype-combined'],
       'data_fisher_fisher': 
       ['fisher_exact_s',
        'fisher_exact_p Fisher-data-combined Fisher-TF-Celltype-combined'],
    }}

    
or_thre = {'CP':9,'SE':6}  
label_name = {'CP':'log2 T-statistics','SE':'log2 Odds Ratio'}

index_cutoff = 18
celltype_cutoff = 6


for indir in indirs[:]:
    # outdir = 'f1_CP_heatmap'
    # os.makedirs(outdir,exist_ok=True)
    infile = '../f2_TFMS_TFBS_CP/{}/per_TF_per_Celltype_{}.csv'.format(indir,indir)
    df = pd.read_csv(infile,index_col = 0)
    # === data for CP and SE enrichment
    for or_key in group_info.keys():
        for fisher_P_key in group_info[or_key].keys():
            s_col = group_info[or_key][fisher_P_key][0]
            p_col = group_info[or_key][fisher_P_key][1]
            df_s = pd.DataFrame()
            df_p = pd.DataFrame()
            basename = '{}_{}_{}'.format(or_key,fisher_P_key,indir)
            # print(p_col,basename)
            for ii in df.index:
                factor = ii.split('_')[0]
                celltype = ii.split('_')[1]
                df_s.loc[factor,celltype] = df.loc[ii,s_col]
                df_p.loc[factor,celltype] = df.loc[ii,p_col]
            df_s.to_csv(outdir+os.sep+'_csv/{}_OR.csv'.format(basename))
            df_p.to_csv(outdir+os.sep+'_csv/{}_P.csv'.format(basename))
            
            # ==== select the TF/cell-type with sufficient data
            selected_index = df_s.notnull().sum(axis=1).sort_values(ascending=False).index[:index_cutoff][::-1]
            selected_col = df_s.notnull().sum(axis=0).sort_values(ascending=False).index[:celltype_cutoff]
            df_s = df_s.loc[selected_index,selected_col]
            df_s = np.log2(df_s)
            df_p = df_p.loc[selected_index,selected_col]
            df_p = df_p.replace(0,1e-299)
            
            # == plot the heatmap
            plt.figure(figsize=(4,7))
            width_ratios = [4,.0,0,.5]
            gs = gridspec.GridSpec(3,4,width_ratios=width_ratios,wspace=.15,hspace=.5)
            pal = plt.cm.bwr
            vmax = or_thre[or_key]
            vmin = -1*vmax
            norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
            color_map = matplotlib.cm.ScalarMappable(norm=norm, cmap=pal)
            scatter_scale = 1
            ax = plt.subplot(gs[:,0])   
            for ii in np.arange(len(df_s.index)):
                for jj in np.arange(len(df_s.columns)):
                    odds_ratio = df_s.loc[df_s.index[ii],df_s.columns[jj]]
                    pvalue = df_p.loc[df_s.index[ii],df_s.columns[jj]]
                    if not np.isnan(odds_ratio):
                        color = color_map.to_rgba(odds_ratio)
                        size = -1*np.log10(pvalue)#;print(size)
                        ax.scatter(jj,ii,s = scatter_scale*size, color = color)
                        # ax.text(jj,ii,df_s.index[ii],fontsize=12)
                    else:
                        ax.scatter(jj,ii,marker = 'x',s = scatter_scale*30, color = 'k')
            
            ax.set_yticks(np.arange(len(df_s.index)))
            ax.set_yticklabels(df_s.index,fontsize=13)
            ax.set_xticks(np.arange(len(df_s.columns)))
            ax.set_xticklabels(df_s.columns,rotation=60,fontsize=13)
            ax.set_xlim([-.5,len(df_s.columns)-.5])
            ax.set_title('{}\n{}\n{}'.format(or_key,fisher_P_key,indir),fontsize=13,va='bottom')
            
            # color bar
            ax = plt.subplot(gs[0,3]) 
            cbar = matplotlib.colorbar.ColorbarBase(ax,cmap=pal,norm=norm,orientation='vertical')
            cbar.set_ticks([vmin,0,vmax])
            cbar.set_label(label_name[or_key],labelpad=-60,fontsize=12)
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
            # idrfile = '../f2_TF_motif_IDR_SE_mutation/data/13059_2021_2456_MOESM2_ESM.xlsx'
            # idrdf = pd.read_excel(idrfile,sheet_name='condition 2(1,6-HD-2)',index_col=0)
            # vmin=-2;vmax=2
            # ax = plt.subplot(gs[:,1]) 
            # plot_data = idrdf.loc[df_s.index[::-1],['AICAP']]
            # plot_data = np.log2(plot_data)
            # g = sns.heatmap(plot_data,cmap=plt.cm.PiYG,cbar=False\
            #                 ,vmax=vmax,vmin=vmin\
            #                 ,ax=ax\
            #                 ,xticklabels=True,yticklabels=False)
            # ax.set_xticklabels(['AICAP'],rotation=60)
            # # color bar
            # ax = plt.subplot(gs[2,3]) 
            # norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
            # cbar = matplotlib.colorbar.ColorbarBase(ax,cmap=plt.cm.PiYG,norm=norm,orientation='vertical')
            # # cbar.set_ticks([-1,0,2])
            # cbar.set_label('log2 AICAP',labelpad=-55,fontsize=13)
            
            
            plt.savefig('{}/{}.pdf'.format(outdir,basename),bbox_inches='tight',pad_inches=0.02,transparent=True)
            plt.show()
            plt.close()
            
            






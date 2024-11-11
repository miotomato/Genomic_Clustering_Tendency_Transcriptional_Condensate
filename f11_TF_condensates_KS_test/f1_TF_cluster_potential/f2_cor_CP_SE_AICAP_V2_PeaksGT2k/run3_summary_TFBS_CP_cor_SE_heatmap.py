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



outdir = 'f3_TFBS_CP_heatmap'
os.makedirs(outdir+os.sep+'_csv',exist_ok=True)
indirs = ['CP_TFBS_all_vs_TFMS','CP_TFBS_overlap_motif_vs_TFMS','CP_TFBS_NOT_overlap_motif_vs_TFMS']

group_info = {
'CP_KStest': {'data_median':
       ['log10-dis ks_2samp-s signed',
        'log10-dis ks_2samp-p signed'],
       'data_fisher':
       ['log10-dis ks_2samp-s signed',
        'log10-dis ks_2samp-p signed Fisher-TF-Celltype-combined'],
    },
'CP_Ttest': {'data_median':
       ['log10-dis T-test-s',
        'log10-dis T-test-p'],
       'data_fisher':
       ['log10-dis T-test-s',
        'log10-dis T-test-p Fisher-TF-Celltype-combined'],
    },
'CP_RankSum': {'data_median':
       ['log10-dis Wilcoxon-rank-sum-s',
        'log10-dis Wilcoxon-rank-sum-p'],
       'data_fisher':
       ['log10-dis Wilcoxon-rank-sum-s',
        'log10-dis Wilcoxon-rank-sum-p Fisher-TF-Celltype-combined'],
    },
'SE': {'data_median':
       ['fisher_exact_s',
        'fisher_exact_p'],
       'data_fisher': 
       ['fisher_exact_s',
        'fisher_exact_p Fisher-TF-Celltype-combined'],
    }}

    
or_thre = {'CP_KStest':.3,
           'CP_Ttest':60,
           'CP_RankSum':60,
           'SE':4}  
label_name = {'CP_KStest':'TFBS CP',
              'CP_Ttest':'TFBS CP',
              'CP_RankSum':'TFBS CP',
              'SE':'log2 OR'}
cmp_dict   = {'CP_KStest':plt.cm.PiYG_r, 
              'CP_Ttest':plt.cm.PiYG_r,
              'CP_RankSum':plt.cm.PiYG_r,
              'SE':plt.cm.bwr}

index_cutoff = 20
column_cutoff = 6


for indir in indirs[:]:
    # outdir = 'f1_CP_heatmap'
    # os.makedirs(outdir,exist_ok=True)
    infile = '../f1_TFMS_TFBS_CP_V2_PeaksGT2k/{}/per_TF_per_Celltype_{}.csv'.format(indir,indir)
    df = pd.read_csv(infile,index_col = 0)
    
    # === data for CP and SE enrichment
    for or_key in list(group_info.keys())[:]:
        for fisher_P_key in list(group_info[or_key].keys())[1:]:
            s_col = group_info[or_key][fisher_P_key][0]
            p_col = group_info[or_key][fisher_P_key][1]
            df_s = pd.DataFrame()
            df_p = pd.DataFrame()
            basename = '{}_{}_{}'.format(fisher_P_key,indir,or_key,)
            # print(p_col,basename)
            for ii in df.index:
                factor = ii.split('_')[1]
                celltype = ii.split('_')[0]
                df_s.loc[factor,celltype] = df.loc[ii,s_col]
                df_p.loc[factor,celltype] = df.loc[ii,p_col]
            df_s.to_csv(outdir+os.sep+'_csv/{}_statistics.csv'.format(basename))
            df_p.to_csv(outdir+os.sep+'_csv/{}_P.csv'.format(basename))
            # ==== z-score of statistics
            df_s_z = scipy.stats.zscore(df_s,nan_policy='omit',axis=1)
            df_s_z = pd.DataFrame(df_s_z,index=df_s.index,columns=df_s.columns)
            df_s_z.to_csv(outdir+os.sep+'_csv/{}_statistics_zscored.csv'.format(basename))
            # ==== select the TF/cell-type with sufficient data
            selected_index = df_s.notnull().sum(axis=1).sort_values(ascending=False).index[:index_cutoff][::-1]
            selected_col = df_s.notnull().sum(axis=0).sort_values(ascending=False).index[:column_cutoff]
            df_s = df_s.loc[selected_index,selected_col]
            # selected_index_sort = df_s.median(axis=1).sort_values(ascending=True).index
            # selected_col_sort = df_s.notnull().median(axis=0).sort_values(ascending=True).index
            # df_s = df_s.loc[selected_index_sort,selected_col_sort]
            if or_key == 'SE':
                df_s = np.log2(df_s)
            # df_p = df_p.loc[selected_index_sort,selected_col_sort]
            df_p = df_p.loc[selected_index,selected_col]
            df_p = df_p.replace(0,1e-299)
            
            # == plot the heatmap
            plt.figure(figsize=(3.6,7.5))
            width_ratios = [4,.0,0,.4]
            height_ratios = [1,1,.8]
            gs = gridspec.GridSpec(3,4,width_ratios=width_ratios,height_ratios = height_ratios,
                                   wspace=.15,hspace=.5)
            pal = cmp_dict[or_key]
            vmax = or_thre[or_key]
            vmin = -1*vmax 
            norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
            color_map = matplotlib.cm.ScalarMappable(norm=norm, cmap=pal)
            scatter_scale = .9
            ax = plt.subplot(gs[:,0])   
            for ii in np.arange(len(df_s.index)):
                for jj in np.arange(len(df_s.columns)):
                    statistics = df_s.loc[df_s.index[ii],df_s.columns[jj]]
                    pvalue = df_p.loc[df_s.index[ii],df_s.columns[jj]]
                    if not np.isnan(statistics):
                        color = color_map.to_rgba(statistics)
                        size = -1*np.log10(pvalue)#;print(size)
                        ax.scatter(jj,ii,s = scatter_scale*size, color = color)
                        # ax.text(jj,ii,df_s.index[ii],fontsize=12)
                    else:
                        ax.scatter(jj,ii,marker = 'x',s = scatter_scale*30, color = 'k')
            
            ax.set_yticks(np.arange(len(df_s.index)))
            ax.set_yticklabels(df_s.index,fontsize=13)
            ax.set_xticks(np.arange(len(df_s.columns)))
            ax.set_xticklabels(df_s.columns,rotation=60,fontsize=13)
            ax.set_xlim([-.6,len(df_s.columns)-.4])
            ax.set_ylim([-.6,len(df_s.index)-.4])
            ax.set_title('{}\n{}\n{}\n'.format(or_key,fisher_P_key,indir),fontsize=13,va='bottom')
            
            # color bar
            ax = plt.subplot(gs[0,3]) 
            cbar = matplotlib.colorbar.ColorbarBase(ax,cmap=pal,norm=norm,orientation='vertical')
            cbar.set_ticks([vmin,0,vmax])
            cbar.set_label(label_name[or_key],labelpad=-48 if or_key == 'SE' else -55,
                           fontsize=12,va='center')
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
            
            






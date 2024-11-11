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
matplotlib.rcParams['font.size']=16
import seaborn as sns
sns.set(font_scale=1.3)
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

import warnings
warnings.filterwarnings("ignore")


def box_plot(pos_len,df,mark_p=False):
    positions = np.arange(pos_len)
    box_step = int(len(vals)/len(positions))
    box_vals = []
    for ii in positions:
        box_val = vals[ii*box_step:(ii+1)*box_step]
        box_vals.append(box_val)
        s,p = scipy.stats.ttest_1samp(box_val,0)#;print(s,p)
        p_val = '*' if p<.05 else ''
        if mark_p == True:
            ax.text(ii, np.mean(box_val), p_val, fontsize=30,c='r',ha='center',va='center')
    g = plt.boxplot(box_vals,positions=positions,widths = .6,patch_artist=True,\
                boxprops=dict(color='k',facecolor='w',fill=None,lw=1),\
                medianprops=dict(color='grey'),showfliers=False)  




def scatter_plot(x,y,xlabel,ylabel,figname,title):

    # ==== scatter plot with density
    plt.figure(figsize=(3,3))
    plt.scatter(x,y,c='k',s=9)
    
    label_i=0
    label_indexes = x[x<0].sort_values(ascending=True)[:half_len].index
    # label_indexes = x.loc[label_indexes].sort_values(ascending=True).index
    for label_index in label_indexes:
        plt.scatter(x[label_index],y[label_index],
                    c=colors1[label_i],
                    s=25,label=label_index)
        label_i+=1
    
    label_i=0
    label_indexes = x[x>0].sort_values(ascending=False)[:half_len].index
    # colors2 = colors2[-5:]
    for label_index in label_indexes[::-1]:
        plt.scatter(x[label_index],y[label_index],
                    c=colors2[label_i],
                    s=25,label=label_index)
        label_i+=1
    
    plt.legend(bbox_to_anchor=[1,1.05],
           markerscale=1.2,fontsize=12,borderaxespad=0.2,labelspacing=.2,
           handletextpad=0.2,handlelength=1.,loc="upper left",frameon=False)
    # ==== linear regression
    # plt.axhline(y=0,color='grey',lw=.5,ls='--')
    plt.axvline(x=0,color='grey',lw=.5,ls='--')
    plt.xlabel(xlabel)
    # plt.ylabel(ylabel)
    plt.title(title,ha='center')
    plt.ylim([ylim_min,ylim_max])
    # plt.axes().set_yticks([0,25,50,75,100,125])
    plt.xlim([-10,10])
    plt.axes().set_xticklabels(['-10'.rjust(10),0,'10'.ljust(7)])
    # if not title=='GM12878':
        # plt.yticks([])
        # plt.ylabel('')
    
    plt.savefig('{}'.format(figname),bbox_inches='tight',pad_inches=0.02,transparent=True)
    plt.show()
    plt.close()
    # return r_value, p_value



def return_colors(pal,half_len,a,b):
    norm = matplotlib.colors.Normalize(vmin=a*half_len,vmax=b*half_len)
    color_map = matplotlib.cm.ScalarMappable(norm=norm, cmap=pal)
    colors = [color_map.to_rgba(i) for i in np.arange(half_len)]
    return colors




# ==== main
# outdir = 'f8bj_per_CT_cor_AICAP_cor_CP_sep_by_num_all_TFs_negLogAICAP'
# os.makedirs(outdir+os.sep+'_csv',exist_ok=True)

# ==== read AICAP data
idrfile = '../../data/public/13059_2021_2456_MOESM2_ESM.xlsx'
idrdf = pd.read_excel(idrfile,sheet_name='condition 2(1,6-HD-2)',index_col=0)

# ==== TFMS CP
# tfms = pd.read_csv('data/TFMS_CP_SE_enrich.csv',index_col=0)
# tfms.index = [i.split('_')[0] for i in tfms.index]

# ==== TFBS info
tfbs_dir = 'f3_TFBS_CP_heatmap/_csv/'
tfbs_cp_types = ['CP_TFBS_all_vs_TFMS',
                 'CP_TFBS_overlap_motif_vs_TFMS',
                 'CP_TFBS_NOT_overlap_motif_vs_TFMS']
tfbs_cp_type = tfbs_cp_types[0]


half_len = 10
# pal1 = plt.cm.RdPu_r
# pal1 = plt.cm.GnBu
pal1 = plt.cm.cool_r
# pal1 = plt.cm.summer
colors1 = return_colors(pal1,half_len,-.1,1.1)
# colors2 = return_colors(pal2,half_len,-.5,1.1)
colors = colors1
# colors = plt.cm.tab20
# cmap = plt.cm.PiYG
# vmax=11
# vmin=-1*vmax; 
# norm = matplotlib.colors.Normalize(vmin=vmin,vmax=vmax)
# color_map = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)

# ylim_min = -.5
# ylim_max = 130
ylim_min = -.05
ylim_max = .35


profileTypes = ['with_motif_SE','with_motif_only']
for profileType in profileTypes[:]:
    # tfbs_dir = 'f3_TFBS_CP_heatmap/_csv/'
    tfbs_dir = 'f3_TFBS_CP_heatmap_{}/_csv/'.format(profileType)
    outdir = 'f8bj_per_CT_cor_AICAP_cor_CP_sep_by_num_all_TFs_negLogAICAP_{}'.format(profileType)
    os.makedirs(outdir+os.sep+'_csv',exist_ok=True)
    tfbs_cp_s = pd.read_csv('{}/data_fisher_{}_CP_KStest_statistics.csv'.format(tfbs_dir,tfbs_cp_type),index_col=0)
    tfbs_se_s = pd.read_csv('{}/data_fisher_{}_SE_statistics.csv'.format(tfbs_dir,tfbs_cp_type),index_col=0)
    # tfms = tfms.loc[tfbs_cp_s.index]
    ct_df = pd.DataFrame()
    x_all,y_all,z_all = [],[],[]
    for ct in tfbs_cp_s.columns[:]:
        y = tfbs_cp_s[ct].dropna()
        if len(y)<1:
            continue
        x = idrdf.loc[y.index.intersection(idrdf.index)].AICAP
        # keep only those with AICAP <1
        x = x[x<1]
        if len(x.dropna())<1:
            continue
        # print(ct,len(x.dropna()))
        ct_df.loc[ct,'TF num']=len(x.dropna())
        # ==== save the TFBS CP SE AICAP info for each cell type
        df = pd.concat([y,tfbs_se_s[ct].loc[y.index],np.log2(x)],axis=1)
        df.columns = ['TFBS CP','log2 OR','log2 AICAP',]
        df.sort_values(by='TFBS CP',ascending=False).round(2).to_csv(outdir+os.sep+'_csv/{}.csv'.format(ct))
        
        # ==== and plot TFBS CP vs. AICAP for each cell type
        # x = x[x<1]
        x = -1*np.log2(x.dropna())
        y = y.loc[x.index]
        xlabel = 'log2 AICAP'
        ylabel = 'TFBS CP'
        # figname = outdir+os.sep+'{}_{}.pdf'.format(tfbs_cp_type,ct)
        # scatter_plot(x,y,xlabel,ylabel,figname,ct)
        x_all = x.append(x_all)
        y_all = y.append(y_all)
        # z_all = z.append(z_all)   
        
    # ==== data merged from all cell types
    df = pd.concat([x_all,y_all],axis=1)
    df.columns = ['log2 AICAP','TFBS CP']
    df = df.sort_values(by='TFBS CP',ascending=True)
    df.to_csv(outdir+os.sep+'all_{}.csv'.format(tfbs_cp_type))
    ct_df.to_csv(outdir+os.sep+'all_{}_CountCT.csv'.format(tfbs_cp_type))
    pd.Series(df.index.unique()).to_csv(outdir+os.sep+'all_{}_uniq.csv'.format(tfbs_cp_type),header=None)
    # uniq_file = outdir+os.sep+'all_{}.csv'.format(tfbs_cp_type))
    
    
    # ==== merged plot of AICAP vs. TFBS CP
    plt.figure(figsize=(4,4))
    ax = plt.axes()
    # df = df.dropna()
    x,y = df['log2 AICAP'],df['TFBS CP']
    ax.scatter(x,y,c='k',s=9)
    label_index = df['log2 AICAP'].drop_duplicates().sort_values(ascending=False)[:half_len].index
    # colors = ['r']*len(label_index)
    # colors = ['','','','','','','']
    for ii in np.arange(len(label_index)):
        ax.scatter(x[label_index[ii]],y[label_index[ii]],
                   # color = colors[ii],label=label_index[ii],s=9)
                   color = plt.cm.tab20(ii),label=label_index[ii],s=9)
    ax.legend(bbox_to_anchor=[1,1],
           markerscale=2,fontsize=14,borderaxespad=0.3,labelspacing=.3,
           handletextpad=0.3,handlelength=1.,loc="upper left",frameon=False)
    # slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)   
    # x_sort = np.sort(x)
    # ax.plot(x_sort,x_sort*slope+intercept,c = 'grey',ls='--',lw=.6)
    s,p = stats.spearmanr(x,y) ;print(s,p)   
    ax.text(.5,.15,'$r = {:.2f}$'.format(s),fontsize=14,transform=ax.transAxes)
    ax.text(.5,.05,'$p = {:.1e}$'.format(p),fontsize=14,transform=ax.transAxes)
    # s,p = stats.pearsonr(x,y)    
    # ax.text(.1,.05,r'$\rho = {:.2f}$'.format(s),fontsize=11,transform=ax.transAxes)
    # ax.axhline(y=0,color='k',lw=1.2,ls='--')
    ax.axvline(x=0,color='grey',lw=.8,ls='--')
    ax.set_xlabel('-log2 AICAP')
    ax.set_ylabel('TFBS CP')
    ax.set_xlim([0,10])
    # ax.set_xticklabels(['-10'.rjust(10),0,'10'.ljust(7)])
    # ax.set_ylim([-5,130])
    # ax.set_yticks([0,25,50,75,100,125])
    ax.set_ylim([ylim_min,0.3])
    # ax.set_yticks([0,25,50,75,100,125])
    # ax.set_title('All cell types')
    figname = outdir+os.sep+'all_{}_scatter.pdf'.format(tfbs_cp_type) 
    plt.savefig('{}'.format(figname),bbox_inches='tight',pad_inches=0.02,transparent=True)
    plt.show()
    plt.close()
    


    # ==== box plot
    plt.figure(figsize=(4,4))
    ax = plt.axes()
    # ax = plt.subplot(gs[0,1])
    box_step = [0,
                int(.25*df.shape[0]),
                int(.5*df.shape[0]),
                int(.75*df.shape[0]),
                int(1*df.shape[0])+1]
    yticklabels = ['First quartile',
                    'Second quartile',
                    'Third quartile',
                    'Fourth quartile']
    # box_step = [0, .1, .2, .4,]
    # yticklabels = ['[0,0.1)',
    #                 '[0.1,0.2)',
    #                 '≥ 0.2']
    # colors = ['k']*len(box_step)
    positions = np.arange(len(box_step)-1)
    box_vals = []
    for ii in positions:
        df_ii = df.iloc[box_step[ii]:box_step[ii+1]]
        # df[(df['TFBS CP']>=box_step[ii]) & (df['TFBS CP']<box_step[ii+1])]
        box_val = df_ii['log2 AICAP'].dropna().values
        # box_val = np.log2(box_val)
        box_vals.append(box_val)
        s,p = scipy.stats.ttest_ind(box_val,box_vals[0],alternative='greater');print(s,p)
        p_val = '*' if (p<.05 and s>0) else ''
        p_val = '{:.0e}'.format(p) if (p<.05 and s>0) else ''
        # ax.text(np.percentile(box_val,50), box_step[ii]+0.03,  p_val, fontsize=25,c='red',ha='center',va='center')
        ax.text(ii, 7.7,  p_val, fontsize=14,c='k',ha='center',va='center')

    box_color = 'gray'
    g = ax.boxplot(box_vals,
                    positions = positions,
                    # positions=np.array(box_step[:-1])+0.05,
                    widths = .3,patch_artist=True,\
                boxprops=dict(color=box_color,facecolor='w',fill=None,lw=1),\
                capprops=dict(color=box_color),\
                whiskerprops=dict(color=box_color),\
                medianprops=dict(color='red',lw=2),showfliers=False,vert=True)  
    g = ax.violinplot(box_vals,
                   positions = positions,showmedians=False,
                   # positions=np.array(box_step[:-1])+0.05,
                   widths = .6,vert=True,showextrema=False,
                   )
                # boxprops=dict(color='k',facecolor='w',fill=None,lw=1),\
                # medianprops=dict(color='grey'),showfliers=False,vert=False)  
    # ax.axvline(x=0,color='grey',lw=.8,ls='--')
    # ax.set_xlim([-10,10])
    # ax.set_xticklabels(['-10'.rjust(10),0,'10'.ljust(7)])
    # ax.set_ylim([ylim_min,ylim_max])
    ax.tick_params(axis='x',length=0 ,right=True,left=False,labelright=True,labelleft=False)
    # ax.set_yticklabels(yticklabels)
    ax.set_xticklabels([])
    ax.set_ylabel('-log2 AICAP')
    figname = outdir+os.sep+'all_{}_box.pdf'.format(tfbs_cp_type) 
    
    plt.savefig('{}'.format(figname),bbox_inches='tight',pad_inches=0.02,transparent=True)
    plt.show()
    plt.close()
    
    


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
sns.set(font_scale=1.1)
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
            ax.text(ii, np.mean(box_val), p_val, fontsize=30,c='red',ha='center',va='center')
    g = plt.boxplot(box_vals,positions=positions,widths = .6,patch_artist=True,\
                boxprops=dict(color='k',facecolor='w',fill=None,lw=1),\
                medianprops=dict(color='grey'),showfliers=False)  




def scatter_plot(x,y,z,xlabel,ylabel,figname,title):

    # ==== scatter plot with density
    plt.figure(figsize=(2.6,3.6))
    plt.scatter(x,y,c='lightgrey',s=9)
    
    vmin=-2; vmax=2
    z = z.dropna()
    z = np.log2(z).clip(vmin,vmax)
    for label_index in z.index:
        plt.scatter(x[label_index],y[label_index],c=z[label_index],
                    vmin=-2, vmax=2,
                    cmap=plt.cm.PiYG,
                    s=35,label=label_index)
        # label_i+=1
    # plt.legend(bbox_to_anchor=[.98,1.02],
    #        markerscale=1.2,fontsize=9,borderaxespad=0.2,labelspacing=.2,
    #        handletextpad=0.2,handlelength=1.,loc="upper left",frameon=False)
    # ==== linear regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)       
    x_sort = np.sort(x)
    plt.plot(x_sort,x_sort*slope+intercept,c = 'grey',ls='--',lw=.6)
    plt.text(.62,.88,'$R^2$ = {:.2f}'.format(r_value**2),fontsize=12,transform=plt.axes().transAxes)
    plt.axhline(y=0,color='k',lw=1.2,ls='--')
    plt.axvline(x=0,color='k',lw=1.2,ls='--')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.savefig('{}'.format(figname),bbox_inches='tight',pad_inches=0.02,transparent=True)
    plt.show()
    plt.close()
    return r_value, p_value


# ==== color panel
# colors = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple',
#           'tab:brown','tab:pink','tab:grey','tab:olive','tab:cyan','gold','lightgreen']

# half_len = 7
# norm = matplotlib.colors.Normalize(vmin=0,vmax=2*half_len)
# pal = plt.cm.rainbow_r
# color_map = matplotlib.cm.ScalarMappable(norm=norm, cmap=pal)
# colors = [color_map.to_rgba(i) for i in np.arange(2*half_len)]



outdir = 'f7_per_CT_cor_AICAP_sep_by_CP_SE'
os.makedirs(outdir,exist_ok=True)

# ==== read AICAP data
idrfile = 'data/13059_2021_2456_MOESM2_ESM.xlsx'
idrdf = pd.read_excel(idrfile,sheet_name='condition 2(1,6-HD-2)',index_col=0)

# ==== TFMS CP
tfms = pd.read_csv('data/TFMS_CP_SE_enrich.csv',index_col=0)
tfms.index = [i.split('_')[0] for i in tfms.index]

# ==== TFBS info
tfbs_dir = 'f3_TFBS_CP_heatmap/_csv/'
# or_key = ['CP','SE']
# fisher_P_key = ['data_median','data_fisher']
tfbs_cp_types = ['CP_TFBS_vs_TFMS',
                 'CP_TFBS_overlap_motif_vs_TFMS',
                 'CP_TFBS_NOT_overlap_motif_vs_TFMS']

cmap = plt.cm.PiYG

for tfbs_cp_type in tfbs_cp_types[:1]:
    # basename = '{}_{}_{}'.format(or_key,fisher_P_key,indir)
    tfbs_cp_s = pd.read_csv('{}/data_fisher_{}_CP_RankSum_statistics.csv'.format(tfbs_dir,tfbs_cp_type),index_col=0)
    tfbs_se_s = pd.read_csv('{}/data_fisher_{}_SE_statistics.csv'.format(tfbs_dir,tfbs_cp_type),index_col=0)
    # tfms = tfms.loc[tfbs_cp_s.index]
    
    x_all,y_all,z_all = [],[],[]
    for ct in tfbs_cp_s.columns[:]:
        y = tfbs_cp_s[ct].dropna()
        if len(y)<=3:
            continue
        x = tfbs_se_s[ct].loc[y.index]
        x = np.log2(x)
        z = idrdf.loc[y.index].AICAP
        if len(z.dropna())<=3:
            continue
        ylabel = 'TFBS CP'
        xlabel = 'TFBS SE enrichment (log2 OR)'
        figname = outdir+os.sep+'{}_{}.pdf'.format(tfbs_cp_type,ct)
        # r_value, p_value = scatter_plot(x,y,z,xlabel,ylabel,figname,ct)
        x_all = x.append(x_all)
        y_all = y.append(y_all)
        z_all = z.append(z_all)   
    
    df = pd.concat([x_all,y_all,z_all],axis=1)
    df.columns = ['TFBS SE enrichment (log2 OR)','TFBS CP','AICAP']
    df.to_csv(outdir+os.sep+'all_{}.csv'.format(tfbs_cp_type))
    
    
    # ==== add plot of AICAP
    plt.figure(figsize=(2.5,5))
    plt.figure(figsize=(2,2.5))
    height_ratios = [2.2,.1,]
    width_ratios = [.7,.3,.7]
    gs = gridspec.GridSpec(2,3,width_ratios = width_ratios,
                           height_ratios =height_ratios, hspace=1,wspace=.15)
    
    # ==== scatter plot
    vmin=-2; vmax=2
    # ax = plt.subplot(gs[0,0:2])
    # # df = df.dropna()
    # x,y,z = df['TFBS SE enrichment (log2 OR)'],df['TFBS CP'],df['AICAP']
    # ax.scatter(x,y,c='white',s=9)
    # z = z.dropna()
    # z = np.log2(z).clip(vmin,vmax)
    # for label_index in z.index:
    #     plt.scatter(x[label_index],y[label_index],c=z[label_index],
    #                 vmin=-2, vmax=2,
    #                 cmap=cmap,
    #                 s=11,label=label_index)
    # # plt.legend(bbox_to_anchor=[.98,1.02],
    # #        markerscale=1.2,fontsize=9,borderaxespad=0.2,labelspacing=.2,
    # #        handletextpad=0.2,handlelength=1.,loc="upper left",frameon=False)
    # # ==== linear regression
    # slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)       
    # x_sort = np.sort(x)
    # ax.plot(x_sort,x_sort*slope+intercept,c = 'grey',ls='--',lw=.6)
    # ax.text(2.5,120,'$R^2$={:.2f}'.format(r_value**2),fontsize=11,)
    # # ax.axhline(y=0,color='k',lw=1.2,ls='--')
    # # ax.axvline(x=0,color='k',lw=1.2,ls='--')
    # ax.set_xlabel('log2 OR')
    # ax.set_ylabel('TFBS CP')
    # ax.set_ylim([0,130])
    # ax.set_yticks([0,25,50,75,100,125])
    

    # ==== box plot
    ax = plt.subplot(gs[0,2])
    ax = plt.subplot(gs[:,:])
    # box_step = [0,20,40,60,80,200]
    box_step = [0,25,50,75,100,200]
    yticklabels = ['[0,25)',
                   '[25,50)',
                   '[50,75)',
                   '[75,100)',
                   '≥100']
    colors = ['k']*len(box_step)
    positions = np.arange(len(box_step)-1)
    box_vals = []
    for ii in positions:
        df_ii = df[(df['TFBS CP']>=box_step[ii]) & (df['TFBS CP']<box_step[ii+1])]
        box_val = df_ii['AICAP'].dropna().values
        box_val = np.log2(box_val)
        box_vals.append(box_val)
        s,p = scipy.stats.ttest_1samp(box_val,0);print(s,p)
        p_val = '*' if (p<.05 and s<0) else ''
        ax.text(np.median(box_val), box_step[ii]+8,  p_val, fontsize=25,c='red',ha='center',va='center')
    g = ax.boxplot(box_vals,positions=np.array(box_step[:-1])+12.5,widths = 15,patch_artist=True,\
                boxprops=dict(color='k',facecolor='w',fill=None,lw=1),\
                medianprops=dict(color='grey'),showfliers=False,vert=False)  
    ax.axvline(x=0,color='k',lw=1.2,ls='--')
    ax.set_xlim([-10,10])
    ax.set_ylim([-5,130])
    ax.tick_params(axis='y',length=3 ,right=True,left=False,labelright=True,labelleft=False)
    ax.set_yticklabels(yticklabels,fontsize=13)
    ax.set_xlabel('log2 AICAP')
    
    
    figname = outdir+os.sep+'all_{}.pdf'.format(tfbs_cp_type) 
    plt.savefig('{}'.format(figname),bbox_inches='tight',pad_inches=0.02,transparent=True)
    plt.show()
    plt.close()
    
    
    # ==== color bar
    # plt.figure(figsize=(.2,.8))
    plt.figure(figsize=(.8,.2))
    # ax = plt.subplot(gs[1,1]) 
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    cbar = matplotlib.colorbar.ColorbarBase(ax=plt.axes(),cmap=cmap,norm=norm,orientation='horizontal')
    cbar.set_label('log2 AICAP',labelpad=-45,fontsize=12)
    plt.tick_params(length=3,direction='out',labelsize='small')
    plt.savefig(outdir+os.sep+'colorbar.pdf',bbox_inches='tight',pad_inches=0.02,transparent=True)
    plt.show()
    plt.close()
    

    
    
    
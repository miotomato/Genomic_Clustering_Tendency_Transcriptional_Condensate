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
            ax.text(ii, np.mean(box_val), p_val, fontsize=30,c='r',ha='center',va='center')
    g = plt.boxplot(box_vals,positions=positions,widths = .6,patch_artist=True,\
                boxprops=dict(color='k',facecolor='w',fill=None,lw=1),\
                medianprops=dict(color='grey'),showfliers=False)  




def scatter_plot(x,y,xlabel,ylabel,figname,title):

    # ==== scatter plot with density
    plt.figure(figsize=(2,2.5))
    plt.scatter(x,y,c='k',s=9)
    
    # label the TFs
    # label_i=0
    # label_indexes = y.sort_values(ascending=False)[:half_len].index
    # label_indexes = label_indexes.append(x[x>0].sort_values(ascending=False)[:half_len].index)
    # # == if current index is < 2*half_len
    # current_index = set(label_indexes)
    # extra_index_num = 2*half_len - len(current_index) 
    # if extra_index_num>0:
    #     label_indexes = label_indexes.append(x[x<0].loc[y.index.difference(current_index)].sort_values(ascending=True)[:extra_index_num].index)
    # label_indexes = x.loc[set(label_indexes)].sort_values(ascending=True).index
    
    # for label_index in label_indexes:
    #     plt.scatter(x[label_index],y[label_index],
    #                 # c=colors[label_i],
    #                 c=color_map.to_rgba(x[label_index]),
    #                 s=25,label=label_index)
    #     label_i+=1
    
    label_i=0
    label_indexes = x[x<-1].sort_values(ascending=True)[:half_len].index
    label_indexes = y.loc[label_indexes].sort_values(ascending=False).index
    for label_index in label_indexes:
        plt.scatter(x[label_index],y[label_index],
                    c=colors1[label_i],
                    # c=color_map.to_rgba(x[label_index]),
                    s=25,label=label_index)
        label_i+=1
    
    label_i=0
    label_indexes = x[x>1].sort_values(ascending=False)[:half_len].index
    # colors2 = colors2[-5:]
    for label_index in label_indexes[::-1]:
        plt.scatter(x[label_index],y[label_index],
                    c=colors2[label_i],
                    # c=color_map.to_rgba(x[label_index]),
                    s=25,label=label_index)
        label_i+=1
    
    
    
    
    plt.legend(bbox_to_anchor=[.0,-.3],
           markerscale=1.2,fontsize=12,borderaxespad=0.2,labelspacing=.2,
           handletextpad=0.2,handlelength=1.,loc="upper left",frameon=False)
    # ==== linear regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)       
    # x_sort = np.sort(x)
    # plt.plot(x_sort,x_sort*slope+intercept,c = 'grey',ls='--',lw=.6)
    # plt.text(.62,.88,'$R^2$ = {:.2f}'.format(r_value**2),fontsize=12,transform=plt.axes().transAxes)
    plt.axhline(y=0,color='grey',lw=.5,ls='--')
    plt.axvline(x=0,color='grey',lw=.5,ls='--')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title,ha='center')
    plt.ylim([0,130])
    plt.xlim([-10,10])
    if not title=='GM12878':
        plt.yticks([])
        plt.ylabel('')
    
    plt.savefig('{}'.format(figname),bbox_inches='tight',pad_inches=0.02,transparent=True)
    plt.show()
    plt.close()
    return r_value, p_value



def return_colors(pal,half_len,a,b):
    norm = matplotlib.colors.Normalize(vmin=a*half_len,vmax=b*half_len)
    color_map = matplotlib.cm.ScalarMappable(norm=norm, cmap=pal)
    colors = [color_map.to_rgba(i) for i in np.arange(half_len)]
    return colors

# ==== color panel
half_len = 7
pal1 = plt.cm.RdPu_r
pal2 = plt.cm.GnBu
# pal1 = plt.cm.cool_r
# pal2 = plt.cm.summer_r
colors1 = return_colors(pal1,half_len,-.2,1.2)
colors2 = return_colors(pal2,half_len,-.5,1.1)
colors = colors1 + colors2 

norm = matplotlib.colors.Normalize(vmin=-1*half_len,vmax=1*half_len)
pal = plt.cm.PiYG
color_map = matplotlib.cm.ScalarMappable(norm=norm, cmap=pal)
# colors = [color_map.to_rgba(i) for i in np.arange(2*half_len)]


outdir = 'f8_per_CT_cor_AICAP_cor_CP'
os.makedirs(outdir,exist_ok=True)

# ==== read AICAP data
idrfile = 'data/13059_2021_2456_MOESM2_ESM.xlsx'
idrdf = pd.read_excel(idrfile,sheet_name='condition 2(1,6-HD-2)',index_col=0)

# ==== TFMS CP
tfms = pd.read_csv('data/TFMS_CP_SE_enrich.csv',index_col=0)
tfms.index = [i.split('_')[0] for i in tfms.index]

# ==== TFBS info
tfbs_dir = 'f3_TFBS_CP_heatmap/_csv/'
tfbs_cp_types = ['CP_TFBS_vs_TFMS',
                 'CP_TFBS_overlap_motif_vs_TFMS',
                 'CP_TFBS_NOT_overlap_motif_vs_TFMS']


for tfbs_cp_type in tfbs_cp_types[:1]:
    # basename = '{}_{}_{}'.format(or_key,fisher_P_key,indir)
    tfbs_cp_s = pd.read_csv('{}/data_fisher_{}_CP_RankSum_statistics.csv'.format(tfbs_dir,tfbs_cp_type),index_col=0)
    tfbs_se_s = pd.read_csv('{}/data_fisher_{}_SE_statistics.csv'.format(tfbs_dir,tfbs_cp_type),index_col=0)
    # tfms = tfms.loc[tfbs_cp_s.index]
    
    # x_all,y_all,z_all = [],[],[]
    for ct in tfbs_cp_s.columns[:]:
        y = tfbs_cp_s[ct].dropna()
        if len(y)<=3:
            continue
        x = idrdf.loc[y.index].AICAP
        x = np.log2(x.dropna())
        y = y.loc[x.index]
        if len(x.dropna())<=3:
            continue
        xlabel = 'log2 AICAP'
        ylabel = 'TFBS CP'
        figname = outdir+os.sep+'{}_{}.pdf'.format(tfbs_cp_type,ct)
        r_value, p_value = scatter_plot(x,y,xlabel,ylabel,figname,ct)
        # x_all = x.append(x_all)
        # y_all = y.append(y_all)
        # z_all = z.append(z_all)   
    
    # df = pd.concat([x_all,y_all,z_all],axis=1)
    # df.columns = ['TFBS SE enrichment (log2 OR)','TFBS CP','AICAP']
    # df.to_csv(outdir+os.sep+'all_{}.csv'.format(tfbs_cp_type))
    
    

    
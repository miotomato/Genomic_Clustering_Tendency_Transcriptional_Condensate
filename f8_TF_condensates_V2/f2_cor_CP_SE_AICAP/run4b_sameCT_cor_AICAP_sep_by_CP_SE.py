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
sns.set(font_scale=1.0)
sns.set_style("whitegrid", {'axes.grid' : False})
import scipy
import scipy.optimize
sns.set_style("ticks")
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams["mathtext.rm"] = "Arial"
from scipy.interpolate import interpn
from scipy.stats import gaussian_kde



def mark_pvalue(compr_pos,positions,box_vals):
    s,p = stats.ttest_ind(box_vals[compr_pos[0]],box_vals[compr_pos[1]],nan_policy='omit')
    y, h, col = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),95)*.95 ,1.05, 'k'
    y2 = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),5)*0.99
    x1,x2 = positions[compr_pos[0]],positions[compr_pos[1]]
    indicator = 'down' if s>0 else 'up'
    p_label='{} \n {:.1e}'.format(indicator,p)
    if p_label[-2]=='0':
        p_label = p_label[:-2]+p_label[-1]
    if p>=0.05:
        p_label = 'n.s.'
    if compr_pos[2] == 't':
        plt.plot([x1*1.03, x1*1.03, x2*0.97, x2*0.97], [y, y*h, y*h, y], lw=1, c=col)
        plt.text((x1+x2)*.5, y*h, p_label, ha='center', va='bottom', color=col,fontsize=12)
    else:
        plt.plot([x1*1.03, x1*1.03, x2*0.97, x2*0.97], [y2, y2*.91, y2*.91, y2], lw=1, c=col)
        plt.text((x1+x2)*.5, y2*.95, p_label, ha='center', va='top', color=col,fontsize=12)



def scatter_plot(x,y,z,xlabel,ylabel,figname,title,tfbs_per_data):

    # ==== scatter plot with density
    plt.figure(figsize=(2.6,2.6))
    x_linear,y_linear = [],[]
    for tf in y.index:
        per_data_index = [i for i in tfbs_per_data.index if re.search('_{}_{}'.format(tf,title),i)]
        x_new = tfbs_per_data.loc[per_data_index]['fisher_exact_s']
        x_new = np.log2(x_new).values
        y_new = tfbs_per_data.loc[per_data_index]['log10-distance CP t-stats'].values
        plt.scatter(x_new,y_new,color='lightgrey',s=8)
        x_linear = np.append(x_linear,x_new)
        y_linear = np.append(y_linear,y_new)
    # plt.scatter(x,y,c='lightgrey',s=9)
    
    vmin=-2; vmax=2
    z = z.dropna()
    z = np.log2(z)#.clip(vmin,vmax)
    df_out = pd.DataFrame()
    for label_index in z.index:
        per_data_index = [i for i in tfbs_per_data.index if re.search('_{}_{}'.format(label_index,title),i)]
        x_new = tfbs_per_data.loc[per_data_index]['fisher_exact_s']
        x_new = np.log2(x_new).values
        y_new = tfbs_per_data.loc[per_data_index]['log10-distance CP t-stats'].values
        plt.scatter(x_new,y_new,c=[z[label_index]]*len(x_new),
                    vmin=-2, vmax=2,
                    cmap=plt.cm.PiYG,
                    s=35,label=label_index)
        df_tmp = pd.DataFrame(data = {'TFBS_index':tfbs_per_data.loc[per_data_index].index,
                                      'SE enrich (log2 OR)':x_new,
                                      'TFBS CP':y_new,
                                      'log2 AICAP':[z[label_index]]*len(x_new),},
                              index = [label_index]*len(x_new))
        df_out = pd.concat([df_out,df_tmp])
        # label_i+=1
    # plt.legend(bbox_to_anchor=[.98,1.02],
    #        markerscale=1.2,fontsize=9,borderaxespad=0.2,labelspacing=.2,
    #        handletextpad=0.2,handlelength=1.,loc="upper left",frameon=False)
    # ==== linear regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(x_linear, y_linear)       
    x_sort = np.sort(x_linear)
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
    return r_value, p_value,df_out


# ==== color panel
# colors = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple',
#           'tab:brown','tab:pink','tab:grey','tab:olive','tab:cyan','gold','lightgreen']

# half_len = 7
# norm = matplotlib.colors.Normalize(vmin=0,vmax=2*half_len)
# pal = plt.cm.rainbow_r
# color_map = matplotlib.cm.ScalarMappable(norm=norm, cmap=pal)
# colors = [color_map.to_rgba(i) for i in np.arange(2*half_len)]


outdir = 'f4b_sameCT_cor_AICAP_sep_by_CP_SE'
os.makedirs(outdir,exist_ok=True)

# ==== read AICAP data
idrfile = 'data/13059_2021_2456_MOESM2_ESM.xlsx'
idrdf = pd.read_excel(idrfile,sheet_name='condition 2(1,6-HD-2)',index_col=0)

# ==== TFMS CP
tfms = pd.read_csv('data/data_TFMS_CP_motif_Num_Len.csv',index_col=0)
tfms.index = [i.split('_')[0] for i in tfms.index]

# ==== TFBS info
tfbs_dir = '../f1_TF_cluster_potential/f3_CP_cor_SE/f1_CP_heatmap/_csv/'
per_data_dir = '../f1_TF_cluster_potential/f2_TFMS_TFBS_CP/'
or_key = ['CP','SE']
fisher_P_key = ['data_median_fisher','data_fisher_fisher']
tfbs_cp_types = ['CP_TFBS_vs_TFMS',
                 'CP_TFBS_with_motif_vs_TFMS',
                 'CP_TFBS_without_motif_vs_TFMS']


for tfbs_cp_type in tfbs_cp_types[:1]:
    # basename = '{}_{}_{}'.format(or_key,fisher_P_key,indir)
    tfbs_cp_s = pd.read_csv('{}/CP_data_median_fisher_{}_OR.csv'.format(tfbs_dir,tfbs_cp_type),index_col=0)
    tfbs_se_s = pd.read_csv('{}/SE_data_median_fisher_{}_OR.csv'.format(tfbs_dir,tfbs_cp_type),index_col=0)
    tfbs_per_data = pd.read_csv('{}/{}/per_data_{}.csv'.format(per_data_dir,tfbs_cp_type,tfbs_cp_type),index_col=0)
    # tfms = tfms.loc[tfbs_cp_s.index]
    
    for ct in tfbs_cp_s.columns[:]:
        y = tfbs_cp_s[ct].dropna()
        if len(y)<=3:
            continue
        x = tfbs_se_s[ct].loc[y.index]
        x = np.log2(x)
        z = idrdf.loc[y.index].AICAP
        ylabel = 'TFBS CP'
        xlabel = 'TFBS SE enrichment (log2 OR)'
        figname = outdir+os.sep+'{}_{}.pdf'.format(tfbs_cp_type,ct)
        r_value, p_value,df_out = scatter_plot(x,y,z,xlabel,ylabel,figname,ct,tfbs_per_data)
    
        df_out = df_out.sort_values(by='TFBS CP',ascending=False)
        df_out.to_csv(outdir+os.sep+'_{}_{}.csv'.format(tfbs_cp_type,ct))
        
    
    
    
    
    
    
    
    
    
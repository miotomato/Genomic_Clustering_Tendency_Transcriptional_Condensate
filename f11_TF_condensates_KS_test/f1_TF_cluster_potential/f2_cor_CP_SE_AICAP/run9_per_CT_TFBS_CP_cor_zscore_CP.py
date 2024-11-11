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



def scatter_plot(df,xlabel,ylabel,figname,title):

    x = df[xlabel]
    y = df[ylabel]
    # ==== scatter plot with density
    plt.figure(figsize=(2.2,2.2))
    plt.scatter(x,y,c='k',s=9)
    # label the TFs
    label_i=0
    top_y = df['avg rank'].sort_values(ascending=True)[:half_len].index
    # label_indexes = x.loc[top_y.intersection(x[x>0].index)].sort_values(ascending=False).index[:half_len]
    for label_index in top_y:
        plt.scatter(x[label_index],y[label_index],
                    c=colors1[label_i],
                    s=25,label=label_index)
        label_i+=1

    # label_i=0
    # top_y = y[x>0].sort_values(ascending=False)[:half_len].index
    # label_indexes = x.loc[top_y.intersection(x[x>0].index)].sort_values(ascending=False).index[:half_len]
    # for label_index in label_indexes:
    #     plt.scatter(x[label_index],y[label_index],
    #                 c=colors1[label_i],
    #                 s=25,label=label_index)
    #     label_i+=1
    
    # label_i=0
    # top_y = y.sort_values(ascending=False)[:half_len].index
    # label_indexes = x.loc[top_y.intersection(x[x<0].index)].sort_values(ascending=True).index[:half_len]
    # for label_index in label_indexes[::-1]:
    #     plt.scatter(x[label_index],y[label_index],
    #                 c=colors2[label_i],
    #                 s=25,label=label_index)
    #     label_i+=1
    
    plt.legend(bbox_to_anchor=[0,-.3],
           markerscale=1.2,fontsize=11,borderaxespad=0.2,labelspacing=.2,
           handletextpad=0.2,handlelength=1.,loc="upper left",frameon=False)
    # ==== linear regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    s,p = stats.spearmanr(x,y)    
    x_sort = np.sort(x)
    # plt.plot(x_sort,x_sort*slope+intercept,c = 'grey',ls='--',lw=.6)
    plt.text(.6,.06,r'$r = {:.2f}$'.format(r_value),fontsize=11,
             transform=plt.axes().transAxes)
    plt.axhline(y=0,color='grey',lw=.7,ls='--')
    plt.axvline(x=0,color='grey',lw=.7,ls='--')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    # plt.ylim(0)
    # plt.ylim([0,130])
    # if not title=='GM12878':
        # plt.yticks([])
    plt.ylabel('')
    xabs_max = abs(max(plt.axes().get_xlim(), key=abs))
    plt.axes().set_xlim(xmin=-xabs_max, xmax=xabs_max)    
    # plt.xlim([-2.55,2.55])
    # plt.axes().set_xticks([-2,0,2])
    # plt.axes().set_xticklabels(['-200'.rjust(9),0,'200'.ljust(6)])

    # if not title=='PANC-1':
    #     plt.axes().spines['right'].set_visible(False)
    plt.title(title,fontsize=13)
    plt.savefig('{}'.format(figname),bbox_inches='tight',pad_inches=0.02,transparent=True)
    plt.show()
    plt.close()
    


# ==== color panel
colors = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple',
          'tab:brown','tab:pink','tab:grey','tab:olive','tab:cyan','gold','lightgreen']

def return_colors(pal,half_len,a,b,c=1):
    norm = matplotlib.colors.Normalize(vmin=a*half_len,vmax=b*half_len)
    color_map = matplotlib.cm.ScalarMappable(norm=norm, cmap=pal)
    colors = [color_map.to_rgba(i) for i in np.arange(c*half_len)]
    return colors


half_len = 5
# norm = matplotlib.colors.Normalize(vmin=-.5,vmax=2.5*half_len)
# pal = plt.cm.rainbow_r
# color_map = matplotlib.cm.ScalarMappable(norm=norm, cmap=pal)
# colors = [color_map.to_rgba(i) for i in np.arange(2*half_len)]

# pal1 = plt.cm.BuPu_r
# pal2 = plt.cm.BuGn
# colors1 = return_colors(pal1,half_len,-.1,1.1)
# colors2 = return_colors(pal2,half_len,-.1,1.1)
pal1 = plt.cm.RdPu_r
pal2 = plt.cm.GnBu
colors1 = return_colors(pal1,half_len,-.2,1.0)
colors2 = return_colors(pal2,half_len,-.2,1.0)
colors = colors1 + colors2 
# colors = return_colors(plt.cm.summer,0,2,2)

outdir = 'f9_per_CT_TFBS_CP_cor_zscore_CP'

# ==== read AICAP data
# idrfile = 'data/13059_2021_2456_MOESM2_ESM.xlsx'
# idrdf = pd.read_excel(idrfile,sheet_name='condition 2(1,6-HD-2)',index_col=0)

# ==== TFMS CP
tfms = pd.read_csv('data/TFMS_CP_SE_enrich.csv',index_col=0)
tfms.index = [i.split('_')[0] for i in tfms.index]

# ==== TFBS info
tfbs_dir = 'f3_TFBS_CP_heatmap/_csv/'
# or_key = ['CP','SE']
# fisher_P_key = ['data_median','data_fisher']
tfbs_cp_types = ['CP_TFBS_all_vs_TFMS',
                 'CP_TFBS_overlap_motif_vs_TFMS',
                 'CP_TFBS_NOT_overlap_motif_vs_TFMS']

cp_types = ['TFMS_CP','TFBS_CP']

for tfbs_cp_type in tfbs_cp_types[:1]:
    # basename = '{}_{}_{}'.format(or_key,fisher_P_key,indir)
    # tfbs_cp_s = pd.read_csv('{}/data_fisher_{}_CP_RankSum_statistics.csv'.format(tfbs_dir,tfbs_cp_type),index_col=0)
    # tfbs_cp_s_z = pd.read_csv('{}/data_fisher_{}_CP_RankSum_statistics_zscored.csv'.format(tfbs_dir,tfbs_cp_type),index_col=0)
    tfbs_cp_s = pd.read_csv('{}/data_fisher_{}_CP_KStest_statistics.csv'.format(tfbs_dir,tfbs_cp_type),index_col=0)
    tfbs_cp_s_z = pd.read_csv('{}/data_fisher_{}_CP_KStest_statistics_zscored.csv'.format(tfbs_dir,tfbs_cp_type),index_col=0)
    tfbs_se_s = pd.read_csv('{}/data_fisher_{}_SE_statistics.csv'.format(tfbs_dir,tfbs_cp_type),index_col=0)
    tfms = tfms.loc[tfbs_cp_s.index]
    # print(tfbs_se_s.min().min())
    for cp_type in cp_types[1:]:
        os.makedirs(outdir+os.sep+cp_type,exist_ok=True)
        # ==== TFBS CP for each CT
        for ct in tfbs_cp_s.columns[:]:
            x = tfbs_cp_s_z[ct].dropna()
            if len(x)<3:
                continue
            y = tfbs_cp_s[ct].loc[x.index]
            # x = tfbs_se_s[ct].loc[y.index]
            # x = np.log2(x)
            xlabel = 'Z-scored CP'
            ylabel = 'TFBS CP'
            # if cp_type == 'TFMS_CP':
            #     y = tfms['log10-dis Wilcoxon-rank-sum-s'].loc[y.index]
            #     ylabel = 'TFMS CP'
            
            df = pd.concat([x,y],axis=1)
            df.columns = [xlabel,ylabel]
            df['Z-scored CP rank'] = df[xlabel].rank(ascending=False)
            df['TFBS CP rank'] = df[ylabel].rank(ascending=False)
            df['avg rank'] = (df['Z-scored CP rank'] + df['TFBS CP rank'])/2
            df = df.sort_values(by='avg rank')
            df.to_csv(outdir+os.sep+'{}/_{}_{}.csv'.format(cp_type,tfbs_cp_type,ct))
            
            figname = outdir+os.sep+'{}/{}_{}.pdf'.format(cp_type,tfbs_cp_type,ct)
            scatter_plot(df,xlabel,ylabel,figname,ct)
            
    
    
    
    
    
    
    
    
    
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



def scatter_plot(x,y,xlabel,ylabel,figname,title):

    # ==== scatter plot with density
    plt.figure(figsize=(1.8,2.5))
    plt.scatter(x,y,c='k',s=9)
    # label the TFs
    label_i=0
    label_indexes = y.sort_values(ascending=False)[:half_len].index
    label_indexes = label_indexes.append(x.sort_values(ascending=False)[:half_len].index)
    label_indexes = y.loc[set(label_indexes)].sort_values(ascending=False).index
    for label_index in label_indexes:
        plt.scatter(x[label_index],y[label_index],c=colors[label_i],
                    s=25,label=label_index)
        label_i+=1
    plt.legend(bbox_to_anchor=[.0,-.3],
           markerscale=1.2,fontsize=11,borderaxespad=0.2,labelspacing=.2,
           handletextpad=0.2,handlelength=1.,loc="upper left",frameon=False)
    # ==== linear regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)       
    x_sort = np.sort(x)
    plt.plot(x_sort,x_sort*slope+intercept,c = 'grey',ls='--',lw=.6)
    plt.text(.04,.88,'$R^2$ = {:.2f}'.format(r_value**2),fontsize=11,transform=plt.axes().transAxes)
    plt.axhline(y=0,color='grey',lw=.5,ls='--')
    plt.axvline(x=0,color='grey',lw=.5,ls='--')
    plt.ylim([0,135])
    plt.xlim([-230,230])
    plt.axes().set_xticks([-200,0,200])
    plt.axes().set_xticklabels(['-200'.rjust(9),0,'200'.ljust(6)])

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    if not title=='GM12878':
        plt.yticks([])
        plt.ylabel('')

    plt.savefig('{}'.format(figname),bbox_inches='tight',pad_inches=0.02,transparent=True)
    plt.show()
    plt.close()
    


def return_colors(pal,half_len,a,b):
    norm = matplotlib.colors.Normalize(vmin=a*half_len,vmax=b*half_len)
    color_map = matplotlib.cm.ScalarMappable(norm=norm, cmap=pal)
    colors = [color_map.to_rgba(i) for i in np.arange(half_len)]
    return colors

# ==== color panel
# colors = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple',
#           'tab:brown','tab:pink','tab:grey','tab:olive','tab:cyan','gold','lightgreen']

half_len = 6
# norm = matplotlib.colors.Normalize(vmin=0,vmax=2*half_len)
# pal = plt.cm.rainbow_r
# color_map = matplotlib.cm.ScalarMappable(norm=norm, cmap=pal)
# colors = [color_map.to_rgba(i) for i in np.arange(2*half_len)]

pal1 = plt.cm.RdPu_r
pal2 = plt.cm.GnBu
colors1 = return_colors(pal1,half_len,-.2,1.0)
colors2 = return_colors(pal2,half_len,-.2,1.0)
colors = colors1 + colors2 



# == check the TF motifs
outdir = 'f4_TFBS_CP_cor_TFMS_CP'
os.makedirs(outdir+os.sep+'_per_CT',exist_ok=True)

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
tfbs_cp_types = ['CP_TFBS_vs_TFMS',
                 'CP_TFBS_overlap_motif_vs_TFMS',
                 'CP_TFBS_NOT_overlap_motif_vs_TFMS']


for tfbs_cp_type in tfbs_cp_types[:1]:
    # basename = '{}_{}_{}'.format(or_key,fisher_P_key,indir)
    basename = 'data_fisher_{}_CP_RankSum'.format(tfbs_cp_type)
    tfbs_s = pd.read_csv('{}/{}_statistics.csv'.format(tfbs_dir,basename),index_col=0)
    tfbs_p = pd.read_csv('{}/{}_P.csv'.format(tfbs_dir,basename),index_col=0)
    # tfbs_se_s = pd.read_csv('{}/data_fisher_{}_SE_statistics.csv'.format(tfbs_dir,tfbs_cp_type),index_col=0)
    # df = pd.concat([tfms,tfbs_s],axis=1,join='inner')
    tfms = tfms.loc[tfbs_s.index]
    # continue
    # ==== TFBS CP across all Cell Types
    x = tfms['log10-dis Wilcoxon-rank-sum-s']
    y = tfbs_s.mean(axis=1)
    xlabel = 'TFMS CP'
    ylabel = 'avg TFBS CP'
    figname = outdir+os.sep+'{}.pdf'.format(tfbs_cp_type)
    scatter_plot(x,y,xlabel,ylabel,figname,tfbs_cp_type)
    
    # ==== TFs to compare AICAP
    # high_TFMS_CP_TFs = []
    # high_TFBS_CP_TFs = []
    
    # ==== TFBS CP for each CT
    for ct in tfbs_s.columns[:]:
        y = tfbs_s[ct].dropna()
        if len(y)<=3:
            continue
        x = tfms['log10-dis Wilcoxon-rank-sum-s'].loc[y.index]
        xlabel = 'TFMS CP'
        ylabel = 'TFBS CP'
        figname = outdir+os.sep+'_per_CT/{}_{}.pdf'.format(tfbs_cp_type,ct)
        scatter_plot(x,y,xlabel,ylabel,figname,ct)
        
        # high_TFMS_CP_TFs = np.append(high_TFMS_CP_TFs,x.sort_values(ascending=False)[:half_len].index)
        # high_TFBS_CP_TFs = np.append(high_TFBS_CP_TFs,y.sort_values(ascending=False)[:half_len].index)        
    
    # ==== compare AICAP of high_TFMS_CP_TFs vs high_TFBS_CP_TFs
    # plt.figure(figsize=(2.6,2.6))
    # high_TFMS_AICAP = np.log2(idrdf.loc[set(high_TFMS_CP_TFs)].dropna().AICAP)
    # high_TFBS_AICAP = np.log2(idrdf.loc[set(high_TFBS_CP_TFs)].dropna().AICAP)
    # box_vals = [high_TFMS_AICAP,high_TFBS_AICAP]
    # positions = [0,1]
    # box_colors = ['blue','red']
    # g = plt.boxplot(box_vals,positions=positions,widths = .6,patch_artist=True,\
    #             boxprops=dict(color='k',facecolor='w',fill=None,lw=1),\
    #             medianprops=dict(color='grey'),showfliers=False)    
    # for patch, color in zip(g['boxes'], colors):
    #     patch.set_facecolor(color)
    # scatter_X = []
    # for position_id in np.arange(len(positions)):
    #     scatter_x = np.random.normal(positions[position_id],0.07,len(box_vals[position_id]))
    #     plt.scatter(scatter_x,box_vals[position_id],color=box_colors[position_id],s=20,zorder=0,alpha=0.99)
    # for compr_pos in [[0,1,'t']]:
    #     mark_pvalue(compr_pos,positions,box_vals)
    # plt.axes().set_xticklabels(['high TFMS CP','high TFBS CP'],rotation=30, ha='right',fontsize=14,color='k')
    # plt.ylabel('log2 AICAP',fontsize=14)
    # plt.title(tfbs_cp_type)
    # plt.savefig('{}/AICAP_{}.pdf'.format(outdir,tfbs_cp_type),bbox_inches='tight',pad_inches=0.02,transparent=True)
    # plt.show()
    # plt.close()
    # ==== save log2 AICAP score
    # pd.concat([high_TFMS_AICAP,high_TFBS_AICAP]).to_csv('{}/AICAP_{}.csv'.format(outdir,tfbs_cp_type))
        
    
    
    
    
    
    
    
    
    
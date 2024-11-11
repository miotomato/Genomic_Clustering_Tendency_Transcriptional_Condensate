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
matplotlib.rcParams['font.size']=13
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

import warnings
warnings.filterwarnings("ignore")



def mark_pvalue(compr_pos,positions,box_vals):
    s,p = stats.ttest_ind(box_vals[compr_pos[0]],box_vals[compr_pos[1]],nan_policy='omit')
    y, h, col = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),90)*.99 ,1.05, 'k'
    y2 = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),3)*0.99
    x1,x2 = positions[compr_pos[0]],positions[compr_pos[1]]
    p_label='{:.1e}'.format(p)
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
    plt.figure(figsize=(2.2,2.5))
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


from matplotlib_venn import venn3,venn2

def venn3_plot(a,b,c,la,lb,lc,figname):
    plt.figure(figsize=(4,4))
    venn3([a,b,c],set_labels=(la,lb,lc))    
    plt.savefig(figname,bbox_inches='tight',pad_inches=0.1)
    plt.close()



def return_colors(pal,half_len,a,b):
    norm = matplotlib.colors.Normalize(vmin=a*half_len,vmax=b*half_len)
    color_map = matplotlib.cm.ScalarMappable(norm=norm, cmap=pal)
    colors = [color_map.to_rgba(i) for i in np.arange(half_len)]
    return colors




# ==== main
outdir = 'f8c_per_CT_cor_AICAP_cor_CP_sep_by_num_add_RPB'
os.makedirs(outdir+os.sep+'_csv',exist_ok=True)

# ==== read AICAP data
idrfile = '../../data/public/13059_2021_2456_MOESM2_ESM.xlsx'
idrdf_ori = pd.read_excel(idrfile,sheet_name='condition 2(1,6-HD-2)',index_col=0)


rbpfile = '../../data/public/NRG-20-113 Hentze Supplementary Table 1.xlsx'
rbpdf = pd.read_excel(rbpfile,sheet_name='Table',index_col=0)
# rbp = rbpdf[rbpdf['3470_human_RBPs']=='YES'].index
rbp = rbpdf.index
idrdf_w_rbp = idrdf_ori.loc[idrdf_ori.index.intersection(rbp)]
idrdf_wo_rbp = idrdf_ori.loc[idrdf_ori.index.difference(rbp)]

# ==== TFMS CP
tfms = pd.read_csv('data/TFMS_CP_SE_enrich.csv',index_col=0)
tfms.index = [i.split('_')[0] for i in tfms.index]
tfms.rename(index={'RBPJ':'NOTCH1'},inplace=True)


# ==== TFBS info
tfbs_dir = 'f3_TFBS_CP_heatmap/_csv/'
tfbs_cp_types = ['CP_TFBS_all_vs_TFMS',
                 'CP_TFBS_overlap_motif_vs_TFMS',
                 'CP_TFBS_NOT_overlap_motif_vs_TFMS']


half_len = 5
pal1 = plt.cm.RdPu_r
pal2 = plt.cm.GnBu
# pal1 = plt.cm.cool_r
# pal2 = plt.cm.summer_r
colors1 = return_colors(pal1,half_len,-.2,1.2)
colors2 = return_colors(pal2,half_len,-.5,1.1)
colors = colors1 + colors2 
# cmap = plt.cm.PiYG
# vmax=11
# vmin=-1*vmax; 
# norm = matplotlib.colors.Normalize(vmin=vmin,vmax=vmax)
# color_map = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)

# ylim_min = -.5
# ylim_max = 130
ylim_min = -.05
ylim_max = .35

for tfbs_cp_type in tfbs_cp_types[:1]:
    # basename = '{}_{}_{}'.format(or_key,fisher_P_key,indir)
    tfbs_cp_s = pd.read_csv('{}/data_fisher_{}_CP_KStest_statistics.csv'.format(tfbs_dir,tfbs_cp_type),index_col=0)
    # tfbs_cp_s = pd.read_csv('{}/data_fisher_{}_CP_RankSum_statistics.csv'.format(tfbs_dir,tfbs_cp_type),index_col=0)
    tfbs_se_s = pd.read_csv('{}/data_fisher_{}_SE_statistics.csv'.format(tfbs_dir,tfbs_cp_type),index_col=0)
    # tfms = tfms.loc[tfbs_cp_s.index]


a,b,c = set(tfbs_cp_s.index),set(idrdf_ori.index),set(rbpdf.index)
la,lb,lc = 'TF','IDR','RBP'
figname = outdir+os.sep+'venn3.pdf'
venn3_plot(a,b,c,la,lb,lc,figname)



def idrdf_by_rbp(idrdf,basename):
   
    x_all,y_all,z_all = [],[],[]
    for ct in tfbs_cp_s.columns[:]:
        y = tfbs_cp_s[ct].dropna()
        # if len(y)<3:
        #     continue
        x = idrdf.loc[y.index.intersection(idrdf.index)].AICAP
        if len(x.dropna())<1:
            continue
        # ==== save the TFBS CP SE AICAP info for each cell type
        df = pd.concat([y,tfbs_se_s[ct].loc[y.index],np.log2(x)],axis=1)
        df.columns = ['TFBS CP','log2 OR','log2 AICAP',]
        df.sort_values(by='TFBS CP',ascending=False).round(2).to_csv(outdir+os.sep+'_csv/{}.csv'.format(ct))
        
        # ==== and plot TFBS CP vs. AICAP for each cell type
        x = np.log2(x.dropna())
        y = y.loc[x.index]
        xlabel = 'log2 AICAP'
        ylabel = 'TFBS CP'
        figname = outdir+os.sep+'{}_{}.pdf'.format(tfbs_cp_type,ct)
        # scatter_plot(x,y,xlabel,ylabel,figname,ct)
        x_all = x.append(x_all)
        y_all = y.append(y_all)
        # z_all = z.append(z_all)   
        
    # ==== data merged from all cell types
    df = pd.concat([x_all,y_all],axis=1)
    df.columns = ['log2 AICAP','TFBS CP']
    df = df.sort_values(by='TFBS CP',ascending=True)
    df.to_csv(outdir+os.sep+'all_{}.csv'.format(tfbs_cp_type))
    
    
    # ==== merged plot of AICAP vs. TFBS CP
    # plt.figure(figsize=(10,2))
    # height_ratios = [2.2,.1,]
    # width_ratios = [1,1]
    # gs = gridspec.GridSpec(1,2,width_ratios = width_ratios,
    #                        # height_ratios =height_ratios, 
    #                        # hspace=1,
    #                        wspace=1)
    # ==== scatter plot
    plt.figure(figsize=(4,2))
    ax = plt.axes()
    # df = df.dropna()
    x,y = df['log2 AICAP'],df['TFBS CP']
    ax.scatter(x,y,c='k',s=9)
    # slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)   
    # x_sort = np.sort(x)
    # ax.plot(x_sort,x_sort*slope+intercept,c = 'grey',ls='--',lw=.6)
    s,p = stats.spearmanr(x,y)    
    ax.text(.7,.15,r'$r = {:.2f}$'.format(s),fontsize=11,transform=ax.transAxes)
    s,p = stats.pearsonr(x,y)    
    ax.text(.7,.05,r'$\rho = {:.2f}$'.format(s),fontsize=11,transform=ax.transAxes)
    # ax.axhline(y=0,color='k',lw=1.2,ls='--')
    ax.axvline(x=0,color='grey',lw=.8,ls='--')
    ax.set_xlabel('log2 AICAP')
    ax.set_ylabel('TFBS CP')
    ax.set_xlim([-10,10])
    # ax.set_xticklabels(['-10'.rjust(10),0,'10'.ljust(7)])
    # ax.set_ylim([-5,130])
    # ax.set_yticks([0,25,50,75,100,125])
    ax.set_ylim([ylim_min,ylim_max])
    # ax.set_yticks([0,25,50,75,100,125])
    ax.set_title(basename)
    figname = outdir+os.sep+'all_{}_{}_scatter.pdf'.format(tfbs_cp_type,basename) 
    plt.savefig('{}'.format(figname),bbox_inches='tight',pad_inches=0.02,transparent=True)
    plt.show()
    plt.close()
    
  

    # ==== box plot
    plt.figure(figsize=(4.5,2))
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
    colors = ['k']*len(box_step)
    positions = np.arange(len(box_step)-1)
    box_vals = []
    for ii in positions:
        df_ii = df.iloc[box_step[ii]:box_step[ii+1]]
        # df[(df['TFBS CP']>=box_step[ii]) & (df['TFBS CP']<box_step[ii+1])]
        box_val = df_ii['log2 AICAP'].dropna().values
        # box_val = np.log2(box_val)
        box_vals.append(box_val)
        s,p = scipy.stats.ttest_1samp(box_val,0);print(s,p)
        p_val = '*' if (p<.05 and s<0) else ''
        # ax.text(np.percentile(box_val,50), box_step[ii]+0.03,  p_val, fontsize=25,c='red',ha='center',va='center')
        ax.text(8, ii-.3,  p_val, fontsize=30,c='red',ha='center',va='center')

    box_color = 'gray'
    g = ax.boxplot(box_vals,
                    positions = positions,
                    # positions=np.array(box_step[:-1])+0.05,
                    widths = .3,patch_artist=True,\
                boxprops=dict(color=box_color,facecolor='w',fill=None,lw=1),\
                capprops=dict(color=box_color),\
                whiskerprops=dict(color=box_color),\
                medianprops=dict(color='red',lw=2),showfliers=False,vert=False)  
    g = ax.violinplot(box_vals,
                   positions = positions,showmedians=False,
                   # positions=np.array(box_step[:-1])+0.05,
                   widths = .6,vert=False,showextrema=False,
                   )
                # boxprops=dict(color='k',facecolor='w',fill=None,lw=1),\
                # medianprops=dict(color='grey'),showfliers=False,vert=False)  
    ax.axvline(x=0,color='grey',lw=.8,ls='--')
    ax.set_xlim([-10,10])
    # ax.set_xticklabels(['-10'.rjust(10),0,'10'.ljust(7)])
    # ax.set_ylim([ylim_min,ylim_max])
    # ax.tick_params(axis='y',length=3 ,right=True,left=False,labelright=True,labelleft=False)
    ax.set_yticklabels(yticklabels)
    ax.set_xlabel('log2 AICAP')
    ax.set_title(basename)
    figname = outdir+os.sep+'all_{}_{}_box.pdf'.format(tfbs_cp_type,basename) 
    
    plt.savefig('{}'.format(figname),bbox_inches='tight',pad_inches=0.02,transparent=True)
    plt.show()
    plt.close()
    
        

# ==== main    
idrdf_by_rbp(idrdf_w_rbp,'w_RBP')
idrdf_by_rbp(idrdf_wo_rbp,'wo_RBP')



# == compr cp of RBP vs. non-RBP TFs
rbp_tf_df = tfbs_cp_s.loc[tfbs_cp_s.index.intersection(rbpdf.index)]
non_rbp_tfs_df = tfbs_cp_s.loc[tfbs_cp_s.index.difference(rbpdf.index)]
rbp_tf_cps = [i for j in rbp_tf_df.values for i in j if not np.isnan(i)]
non_rbp_tf_cps = [i for j in non_rbp_tfs_df.values for i in j if not np.isnan(i)]

s,p = scipy.stats.ttest_ind(rbp_tf_cps,non_rbp_tf_cps,nan_policy='omit')
print(s,p)


# rbp_tfms_cps = tfms.loc[rbp_tf_df.index]['log10-dis ks_2samp-s signed']
# non_rbp_tfms_cps = tfms.loc[non_rbp_tfs_df.index]['log10-dis ks_2samp-s signed']


# ==== box compar TFBS CPs
plt.figure(figsize=(2.2,2.2))
ax = plt.axes()
# ax = plt.subplot(gs[0,1])
positions = [1,2]
box_vals = [non_rbp_tf_cps,rbp_tf_cps]

box_color = 'k'
g = ax.boxplot(box_vals,
                positions = positions,
                # positions=np.array(box_step[:-1])+0.05,
                widths = .3,patch_artist=True,\
            boxprops=dict(color=box_color,facecolor='w',fill=None,lw=1),\
            capprops=dict(color=box_color),\
            whiskerprops=dict(color=box_color),\
            medianprops=dict(color='red',lw=2),showfliers=False,vert=True)  

mark_pvalue([0,1,'t'],positions,box_vals)
ax.set_xticklabels(['non RBPs (n={})'.format(non_rbp_tfs_df.shape[0]),
                    'RBPs (n={})'.format(rbp_tf_df.shape[0])],rotation=45)
# ax.set_ylim([ylim_min,ylim_max])
# ax.tick_params(axis='y',length=3 ,right=True,left=False,labelright=True,labelleft=False)
# ax.set_yticklabels('')
ax.set_ylabel('TFBS CP')
# ax.set_title(basename)
figname = outdir+os.sep+'RBP_compr_TFBS_CP.pdf' 

plt.savefig('{}'.format(figname),bbox_inches='tight',pad_inches=0.02,transparent=True)
plt.show()
plt.close()



rbp_tfms_cps = tfms.loc[rbp_tf_df.index]['log10-dis ks_2samp-s signed']
non_rbp_tfms_cps = tfms.loc[non_rbp_tfs_df.index]['log10-dis ks_2samp-s signed']


# ==== box compar TFMS CPs
plt.figure(figsize=(2.2,2.2))
ax = plt.axes()
# ax = plt.subplot(gs[0,1])
positions = [1,2]
box_vals = [non_rbp_tfms_cps,rbp_tfms_cps]

box_color = 'k'
g = ax.boxplot(box_vals,
                positions = positions,
                # positions=np.array(box_step[:-1])+0.05,
                widths = .3,patch_artist=True,\
            boxprops=dict(color=box_color,facecolor='w',fill=None,lw=1),\
            capprops=dict(color=box_color),\
            whiskerprops=dict(color=box_color),\
            medianprops=dict(color='red',lw=2),showfliers=False,vert=True)  

mark_pvalue([0,1,'t'],positions,box_vals)
ax.set_xticklabels(['non RBPs (n={})'.format(non_rbp_tfs_df.shape[0]),
                    'RBPs (n={})'.format(rbp_tf_df.shape[0])],rotation=45)
# ax.set_ylim([ylim_min,ylim_max])
# ax.tick_params(axis='y',length=3 ,right=True,left=False,labelright=True,labelleft=False)
# ax.set_yticklabels('')
ax.set_ylabel('TFMS CP')
# ax.set_title(basename)
figname = outdir+os.sep+'RBP_compr_TFMS_CP.pdf' 

plt.savefig('{}'.format(figname),bbox_inches='tight',pad_inches=0.02,transparent=True)
plt.show()
plt.close()




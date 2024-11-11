import os,sys,argparse,glob,re
import numpy as np
import pandas as pd
# import find_overlap_keep_info_NOT_sep_strand_asimport
from collections import Counter
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=14
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams["mathtext.rm"] = "Arial"
import seaborn as sns
sns.set(font_scale=1.1)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks")
from scipy import stats
from scipy.interpolate import interpn
from scipy.stats import gaussian_kde





def scatter_plot(x,y,xlabel,ylabel,figname,title):

    # ==== scatter plot with density
    plt.figure(figsize=(2.6,2.6))
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
    plt.legend(bbox_to_anchor=[.98,1.02],
           markerscale=1.2,fontsize=9,borderaxespad=0.2,labelspacing=.2,
           handletextpad=0.2,handlelength=1.,loc="upper left",frameon=False)
    # ==== linear regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)       
    x_sort = np.sort(x)
    plt.plot(x_sort,x_sort*slope+intercept,c = 'grey',ls='--',lw=.6)
    plt.text(.62,.08,'$R^2$ = {:.2f}'.format(r_value**2),fontsize=12,transform=plt.axes().transAxes)
    plt.axhline(y=0,color='k',lw=1.2,ls='--')
    plt.axvline(x=0,color='k',lw=1.2,ls='--')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.savefig('{}'.format(figname),bbox_inches='tight',pad_inches=0.02,transparent=True)
    plt.show()
    plt.close()
    return r_value, p_value
    




# ==== main 

project_dir='/Volumes/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/'
# project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR'
# atac_overlap_SE_dir = '{}/f7_TF_condensates_test/f6_revised_TCGA_ATAC_cor_SE/f3_clinical_hicor_SE_overlapped_peaks/f1_clinical_at_each_peak/f1_ATAC_overlap_SE_caseID'.format(project_dir)
clinical_per_peak_dir = '{}/f7_TF_condensates_test/f6_revised_TCGA_ATAC_cor_SE/f3_clinical_hicor_SE_overlapped_peaks/f1_clinical_at_each_peak/f2_caseID_each_SE_vs_clinical/'.format(project_dir)
clinical_per_gene_dir = '{}/f7_TF_condensates_test/f6_revised_TCGA_ATAC_cor_SE/f5_clinical_per_gene/f1_clinical_per_gene/'.format(project_dir)
# atac_file = '{}/f7_TF_condensates_test/f6_revised_TCGA_ATAC_cor_SE/data/TCGA/tcga_atac.bed'.format(project_dir)
# atac_df = pd.read_csv(atac_file,sep='\t',header=None,index_col=3)
atac_overlap_SE_TFMS_dir = '{}/f8_TF_condensates_V2/f3_clinical_outcome/data/ATAC_overlap_SE_overlap_TFMS_avg_sig/'.format(project_dir)
rp_dir = '{}/f8_TF_condensates_V2/f3_clinical_outcome/data/ATAC_overlap_SE_overlap_TFMS_avg_sig_RP'.format(project_dir)

# ==== read matched names between TCGA ATAC and SE data
name_match = pd.read_excel('../../data/TCGA/TCGA-ATAC_SE_cancerType_match.xlsx',index_col=0)   
name_match = name_match.dropna()
# ==== read TFBS CP info
tfbs_dir = '../../f1_TF_cluster_potential/f2_cor_CP_SE_AICAP/f3_TFBS_CP_heatmap/_csv/'
tfbs_cp_type = 'CP_TFBS_vs_TFMS'
# cp_type = 'TFBS_CP'
tf_thre = 6
rp_thre = 1
p_thre = 0.05


# outdir = 'f2_clinical_gene_vs_peak_gene_treatLTcontrol'
outdir = 'f2_clinical_gene_vs_peak_gene'
os.makedirs(outdir,exist_ok=True)


# norm = matplotlib.colors.Normalize(vmin=-.5,vmax=2.5*half_len)
# pal = plt.cm.rainbow_r
# color_map = matplotlib.cm.ScalarMappable(norm=norm, cmap=pal)


for cancertype in ['BRCA','COAD'][:]:
    cancertype_SE = name_match.loc[cancertype].SE
    cancertype_SE_rename = name_match.loc[cancertype].SE_rename
    
    df = pd.read_csv('{}/{}_clinical_Gene_vs_Peaks.csv'.format(outdir,cancertype),index_col=0)
    df = df.dropna()
    df = df.sort_values(by='TFBS CP',ascending=False)

    norm = matplotlib.colors.Normalize(vmin=0,vmax=df.dropna().shape[0])
    pal = plt.cm.PiYG
    color_map = matplotlib.cm.ScalarMappable(norm=norm, cmap=pal)
    colors = [color_map.to_rgba(i) for i in np.arange(df.dropna().shape[0])]
    
    x = df['SE-overlapped-peaks #P<0.05']/df['SE-overlapped-peaks total']
    y = df['Target-genes #P<0.05']/df['Target-genes total']
    avg = .5*(x+y)
    # x = df['avg %P<0.05']
    z = df['TFBS CP']
                
    plt.figure(figsize=(2.6,2.6))
    # xy = np.vstack([x,y])
    # z = gaussian_kde(xy)(xy)
    # z[np.where(np.isnan(z))] = 0.0
    # idx = z.argsort()
    # x, y, z = x[idx], y[idx], z[idx]
    # g = plt.scatter(x,y,c=z,cmap = plt.cm.bwr,s=3,marker='o')
    x = x
    y = y
    
    plt.scatter(x,y,s=10,c='k')
    
    # plt.scatter(np.arange(len(x)),avg)
    
    label_i = 0
    for label_index in df.dropna().index:
        plt.scatter(x[label_index],y[label_index],c=colors[label_i],
                s=25,label=label_index)
        if label_i<10:
            plt.text(x[label_index],y[label_index],label_index,fontsize=7)
        label_i+=1
    plt.legend(bbox_to_anchor=[1.1,1.02],
            markerscale=1.2,fontsize=8,borderaxespad=0.2,labelspacing=.2,
            ncol = 2,
            handletextpad=0.2,handlelength=1.,loc="upper left",frameon=False)
        
        
    plt.title(cancertype)
    plt.xlabel('% SE overlapped peaks \n w/ rank-sum p<0.05')
    plt.ylabel('% Target genes \n w/ rank-sum p<0.05')
    figname = '{}/{}.pdf'.format(outdir,cancertype)
    plt.savefig('{}'.format(figname),bbox_inches='tight',pad_inches=0.02,transparent=True)
    plt.show()
    plt.close()

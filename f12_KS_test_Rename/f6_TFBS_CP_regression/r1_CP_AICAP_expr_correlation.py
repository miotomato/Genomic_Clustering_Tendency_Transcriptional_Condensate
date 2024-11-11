import os,sys,argparse,glob,re,bisect
import numpy as np
import pandas as pd
from collections import Counter
import operator
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
#matplotlib.rcParams['agg.path.chunksize'] = 10000
matplotlib.rcParams['font.size']=16
import seaborn as sns
sns.set(font_scale=1.1)
sns.set_style("whitegrid", {'axes.grid' : False})
import scipy
import scipy.optimize
sns.set_style("ticks")
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
from scipy import stats
from scipy.stats import gamma


hg38_chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',
             'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',
             'chr18','chr19','chr20','chr21','chr22','chrX','chrY']



def scatter_plot(df,xlabel,ylabel,outdir):
    
    plt.figure(figsize=(2.6,2.6))
    # xy = np.vstack([x,y])
    # z = gaussian_kde(xy)(xy)
    # z[np.where(np.isnan(z))] = 0.0
    # idx = z.argsort()
    # x, y, z = x[idx], y[idx], z[idx]
    # g = plt.scatter(x,y,c=z,cmap = plt.cm.GnBu_r,s=3,marker='o')
    x = df[xlabel]
    y = df[ylabel]
    plt.scatter(x,y,c='k',s=3)
    # ==== linear regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)       
    x_sort = np.sort(x)
    plt.plot(x_sort,x_sort*slope+intercept,c = 'grey',ls='--',lw=.9)
    plt.text(.55,1.15,'$r$ = {:.2f}'.format(r_value),fontsize=12,transform=plt.axes().transAxes)
    plt.text(.55,1.05,'$p$ = {:.1e}'.format(p_value),fontsize=12,transform=plt.axes().transAxes)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    # plt.axhline(y=0,color='k',lw=.5,ls='--')
    # plt.axvline(x=0,color='k',lw=.5,ls='--')
    figname = '{}/fig_{}_vs_{}.pdf'.format(outdir,xlabel,ylabel)
    plt.savefig('{}'.format(figname),bbox_inches='tight',pad_inches=0.02,transparent=True)
    plt.show()
    plt.close()




outdir = 'f1_CP_AICAP_expr_correlation'
os.makedirs(outdir,exist_ok=True)

project_dir = '/nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/f12_KS_test_Rename'
project_dir = '/Volumes/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/f12_KS_test_Rename'

expr_df = pd.read_csv('{}/f6_TFBS_CP_regression/data/CellType_with_AICAP_TF_GT3_geneExpr_TPM.csv'.format(project_dir),index_col=0)
tfms_cp = pd.read_csv('{}/f5_gamma_fit/f3_TFMS_CP_cor_SE_heatmap/TFMS_CP_SeEnrich_GammaFit.csv'.format(project_dir),index_col=0)
tfbs_cp = pd.read_csv('{}/f1_TF_cluster_potential/f2_cor_CP_SE_AICAP/f3_TFBS_CP_heatmap_with_motif_only/_csv/data_fisher_CP_TFBS_all_vs_TFMS_CP_KStest_statistics.csv'.format(project_dir),index_col=0)
tfbs_zcp = pd.read_csv('{}/f1_TF_cluster_potential/f2_cor_CP_SE_AICAP/f3_TFBS_CP_heatmap_with_motif_only/_csv/data_fisher_CP_TFBS_all_vs_TFMS_CP_KStest_statistics_zscored.csv'.format(project_dir),index_col=0)
tfbs_se = pd.read_csv('{}/f1_TF_cluster_potential/f2_cor_CP_SE_AICAP/f3_TFBS_CP_heatmap_with_motif_only/_csv/data_fisher_CP_TFBS_all_vs_TFMS_SE_statistics.csv'.format(project_dir),index_col=0)
idr_df = pd.read_excel('{}/data/public/13059_2021_2456_MOESM2_ESM.xlsx'.format(project_dir),sheet_name='condition 2(1,6-HD-2)',index_col=0)

tfms_cp.index = [i.split('_')[0] for i in tfms_cp.index]
tfms_cp.rename(index={'RBPJ':'NOTCH1'},inplace=True)


df_all = pd.DataFrame()
for ct in tfbs_cp.columns[:]:
    df_ct = pd.DataFrame()
    tfbs = tfbs_cp[ct].dropna()
    if ct not in expr_df.columns or len(tfbs)<1:
        continue
    idr = idr_df.loc[tfbs.index.intersection(idr_df.index)].AICAP
    idr = idr[idr<1]
    if len(idr.dropna())<1:
        continue
    print(ct)
    
    df_ct['TFMS_CP'] = tfms_cp['log10-dis ks_2samp-s signed'].loc[idr.index]
    df_ct['TFMS_alpha'] = tfms_cp['alpha'].loc[idr.index]  
    df_ct['TFBS_CP'] = tfbs.loc[idr.index]
    df_ct['TFBS_ZCP'] = tfbs_zcp[ct].loc[idr.index]
    a = df_ct['TFBS_CP'].rank(ascending=False)
    b = df_ct['TFBS_ZCP'].rank(ascending=False)
    df_ct['TFBS_CP_avgRank'] = (a+b)/2
    df_ct['TFBS_SE'] = tfbs_se[ct].loc[idr.index]
    df_ct['Expr_TPM'] = expr_df[ct].loc[idr.index]
    df_ct['Log2_TPM'] = np.log2(expr_df[ct].loc[idr.index]+1)
    df_ct['negLog2_AICAP'] = -np.log2(idr)
    df_ct.to_csv('{}/data_{}.csv'.format(outdir,ct))
    df_ct.index = ['{}_{}'.format(ct,ii) for ii in df_ct.index]
    df_all = pd.concat([df_all,df_ct])
    
df_all.to_csv('{}/data_all.csv'.format(outdir))    
    


xlabel = 'TFBS_CP'
ylabel = 'TFMS_CP'
scatter_plot(df_all,xlabel,ylabel,outdir)

# xlabel = 'TFBS_CP'
# ylabel = 'negLog2_AICAP'
# scatter_plot(df_all,xlabel,ylabel,outdir)

xlabel = 'TFMS_CP'
ylabel = 'negLog2_AICAP'
scatter_plot(df_all,xlabel,ylabel,outdir)

xlabel = 'TFBS_CP'
ylabel = 'Log2_TPM'
scatter_plot(df_all,xlabel,ylabel,outdir)

xlabel = 'TFMS_CP'
ylabel = 'Log2_TPM'
scatter_plot(df_all,xlabel,ylabel,outdir)

# xlabel = 'TFBS_SE'
# ylabel = 'Expr_log2_TPM_Plus1'
# scatter_plot(df_all,xlabel,ylabel,outdir)

# xlabel = 'TFBS_CP'
# ylabel = 'TFBS_SE'
# scatter_plot(df_all,xlabel,ylabel,outdir)



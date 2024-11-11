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
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid' : False})
import scipy
import scipy.optimize
sns.set_style("ticks")
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams["mathtext.rm"] = "Arial"
from scipy.interpolate import interpn
from scipy.stats import gaussian_kde




def return_final_ks_statistics(df):
    
    for ii in df.index:
        if df.loc[ii,'log10-dis ks_2samp-s less'] > df.loc[ii,'log10-dis ks_2samp-s greater']:
            df.loc[ii,'log10-dis ks_2samp-s signed'] = df.loc[ii,'log10-dis ks_2samp-s less']
            df.loc[ii,'log10-dis ks_2samp-p signed'] = df.loc[ii,'log10-dis ks_2samp-p less']
        else:
            df.loc[ii,'log10-dis ks_2samp-s signed'] = -1*df.loc[ii,'log10-dis ks_2samp-s greater']
            df.loc[ii,'log10-dis ks_2samp-p signed'] = df.loc[ii,'log10-dis ks_2samp-p greater']
    return df



def scatter_plot(value_dic,key_pairs,mark=False):

    x = value_dic[key_pairs[0]]
    y = value_dic[key_pairs[1]]
    xlabel = ' '.join(key_pairs[0].split('_'))
    ylabel = ' '.join(key_pairs[1].split('_'))
    figname = outdir+os.sep+'{}_vs_{}.pdf'.format(key_pairs[0],key_pairs[1])

    # ==== scatter plot with density
    plt.figure(figsize=(2.6,2.6))
    xy = np.vstack([x,y])
    z = gaussian_kde(xy)(xy)
    z[np.where(np.isnan(z))] = 0.0
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]
    g = plt.scatter(x,y,c=z,cmap = plt.cm.GnBu_r,s=3,marker='o')
    # plt.scatter(x,y,c='k',s=1)
    # ==== linear regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)       
    x_sort = np.sort(x)
    plt.plot(x_sort,x_sort*slope+intercept,c = 'grey',ls='--',lw=.6)
    plt.text(.6,1.01,'$R^2$ = {:.2f}'.format(r_value**2),fontsize=12,transform=plt.axes().transAxes)
    # ==== mark the outliers
    y_new = x*slope+intercept
    marked_indexes = (np.abs((y-y_new))).sort_values(ascending=False)[:6].index
    if mark == True:
        for marked_index in marked_indexes:
            plt.text(x[marked_index],y[marked_index],marked_index.split('_')[0],c='r',fontsize=6)
    plt.axhline(y=0,color='k',lw=1.2,ls='--')
    plt.axvline(x=0,color='k',lw=1.2,ls='--')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig('{}'.format(figname),bbox_inches='tight',pad_inches=0.02,transparent=True)
    plt.show()
    plt.close()
    



# == check the TF motifs
outdir = 'f1_TFMS_CP_cor_length_number'
os.makedirs(outdir,exist_ok=True)
os.makedirs('data',exist_ok=True)

df1 = pd.read_csv('../f1_TFMS_TFBS_CP/CP_TFMS_vs_random/data_CP_TFMS_vs_random.csv',index_col=0)
df1 = return_final_ks_statistics(df1)
df2 = pd.read_csv('../f1_TFMS_TFBS_CP/CP_TFMS_vs_random/data_TFMS_enrich_at_SE.csv',index_col=0)
df = pd.concat([df1,df2],axis=1)
df.to_csv('data/TFMS_CP_SE_enrich.csv')



value_dic = {'Motif_num':df['#TFMS'],
             'Motif_length':df['len-of-TFMS'],
             'T_test_statistics':df['log10-dis T-test-s'],
             'Rank_Sum_statistics':df['log10-dis Wilcoxon-rank-sum-s'],
             'KS_statistics':df['log10-dis ks_2samp-s signed'],
             'P_LT_005_num':df['#p<0.05'],
             'P_LT_005_percentage':100*df['#p<0.05']/df['#TFMS'],
             'P_LT_005_log2_OR_at_SE':np.log2(df['TFMS-p<0.05 enrich-at-SE-fisher-exact-s']),
             }


# == plot 
key_pairs = ['Motif_num','Motif_length']
scatter_plot(value_dic,key_pairs)


key_pairs = ['KS_statistics','Motif_num']
scatter_plot(value_dic,key_pairs)

key_pairs = ['KS_statistics','Motif_length']
scatter_plot(value_dic,key_pairs)


key_pairs = ['T_test_statistics','Motif_num']
scatter_plot(value_dic,key_pairs)

key_pairs = ['T_test_statistics','Motif_length']
scatter_plot(value_dic,key_pairs)

key_pairs = ['T_test_statistics','KS_statistics']
scatter_plot(value_dic,key_pairs,mark=True)


key_pairs = ['Rank_Sum_statistics','Motif_num']
scatter_plot(value_dic,key_pairs)

key_pairs = ['Rank_Sum_statistics','Motif_length']
scatter_plot(value_dic,key_pairs)

key_pairs = ['Rank_Sum_statistics','KS_statistics']
scatter_plot(value_dic,key_pairs,mark=True)

key_pairs = ['Rank_Sum_statistics','T_test_statistics']
scatter_plot(value_dic,key_pairs,mark=True)


key_pairs = ['P_LT_005_num','Motif_num']
scatter_plot(value_dic,key_pairs,mark=True)

key_pairs = ['P_LT_005_num','Motif_length']
scatter_plot(value_dic,key_pairs,mark=False)

key_pairs = ['P_LT_005_num','KS_statistics']
scatter_plot(value_dic,key_pairs,mark=True)

key_pairs = ['P_LT_005_num','T_test_statistics']
scatter_plot(value_dic,key_pairs,mark=True)

key_pairs = ['P_LT_005_num','Rank_Sum_statistics']
scatter_plot(value_dic,key_pairs,mark=True)



key_pairs = ['P_LT_005_percentage','Motif_num']
scatter_plot(value_dic,key_pairs,mark=False)

key_pairs = ['P_LT_005_percentage', 'Motif_length']
scatter_plot(value_dic,key_pairs,mark=False)

key_pairs = ['P_LT_005_percentage', 'KS_statistics']
scatter_plot(value_dic,key_pairs,mark=True)

key_pairs = ['P_LT_005_percentage', 'T_test_statistics']
scatter_plot(value_dic,key_pairs,mark=True)

key_pairs = ['P_LT_005_percentage', 'Rank_Sum_statistics']
scatter_plot(value_dic,key_pairs,mark=True)



key_pairs = ['P_LT_005_log2_OR_at_SE','Motif_num']
scatter_plot(value_dic,key_pairs,mark=False)

key_pairs = ['P_LT_005_log2_OR_at_SE', 'Motif_length']
scatter_plot(value_dic,key_pairs,mark=False)

key_pairs = ['P_LT_005_log2_OR_at_SE', 'KS_statistics']
scatter_plot(value_dic,key_pairs,mark=True)

key_pairs = ['P_LT_005_log2_OR_at_SE', 'T_test_statistics']
scatter_plot(value_dic,key_pairs,mark=True)

key_pairs = ['P_LT_005_log2_OR_at_SE', 'Rank_Sum_statistics']
scatter_plot(value_dic,key_pairs,mark=True)

key_pairs = ['P_LT_005_log2_OR_at_SE', 'P_LT_005_num']
scatter_plot(value_dic,key_pairs,mark=True)

key_pairs = ['P_LT_005_log2_OR_at_SE', 'P_LT_005_percentage']
scatter_plot(value_dic,key_pairs,mark=True)







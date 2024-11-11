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
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid' : False})
import scipy
import scipy.optimize
sns.set_style("ticks")
matplotlib.rcParams["font.sans-serif"] = ["Arial"]

from scipy.stats import gamma


hg38_chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',
             'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',
             'chr18','chr19','chr20','chr21','chr22','chrX','chrY']



def fit_exp_nonlinear(t, y):
    opt_parms, parm_cov = scipy.optimize.curve_fit(model_function, t, y, maxfev=200000)
#    print(opt_parms, parm_cov)
    A,m, K= opt_parms
    return A,m, K

# def fit_exp_linear(t,y,C=0):
#     y = y-C
#     y = np.log(y)
#     K,log_A = np.polyfit(t,y,1)
#     A = np.exp(log_A)
#     return A,K,C

def model_function(t,A,m,K):
    #return A*np.exp(K*t)
#    return A*((t-m)**(K))
    return (1/A)*((t-m)**(K))




# project_dir='/Volumes/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/'
# # project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/'
# infiles = glob.glob('{}/f12_KS_test_Rename/f1_TF_cluster_potential/f0_bedtools_closest/data_TFMS_jarspar//*'.format(project_dir))
# infiles = [i for i in infiles if not re.search('~',i)]
# infiles.sort()
# indir = '{}/f12_KS_test_Rename/f1_TF_cluster_potential/f0_bedtools_closest/data_TFMS_jarspar_by_pvalue/'.format(project_dir)


outdir = 'f7_TFMS_gamma_alpha_by_num_fig'
os.makedirs(outdir,exist_ok=True)
genome_len = 3298912062
num_thres = [2, 5, 10, 20, 50, 100]

df = pd.read_csv('f6_TFMS_gamma_alpha_by_num/TFMS_gamma_alpha_combined.csv',index_col=0)
df = df.sort_values(by='alpha',ascending=True)
df = df.dropna()
colors = ['tab:red','tab:blue','tab:orange','tab:green','tab:purple','tab:brown','tab:cyan']
columns = ['alpha','alpha top2k','alpha top5k','alpha top10k','alpha top20k','alpha top50k','alpha top100k']
xticklabels = ['all','top2k','top5k','top10k','top20k','top50k','top100k']
s=5
alpha=.9

# plot ranked by alpha
plt.figure(figsize=(6,3))
x = np.arange(len(df.index))
a = plt.scatter(x,df['alpha'],c=colors[0],s=s,alpha=alpha,zorder=6)    
b = plt.scatter(x,df['alpha top2k'],c=colors[1],s=s,alpha=alpha,zorder=0)
c = plt.scatter(x,df['alpha top5k'],c=colors[2],s=s,alpha=alpha,zorder=1)
d = plt.scatter(x,df['alpha top10k'],c=colors[3],s=s,alpha=alpha,zorder=2)
e = plt.scatter(x,df['alpha top20k'],c=colors[4],s=s,alpha=alpha,zorder=3)
f = plt.scatter(x,df['alpha top50k'],c=colors[5],s=s,alpha=alpha,zorder=4)
g = plt.scatter(x,df['alpha top100k'],c=colors[6],s=s,alpha=alpha,zorder=5)

plt.legend([a,b,c,d,e,f,g],['all','top2k','top5k','top10k','top20k','top50k','top100k'],
           bbox_to_anchor=[1.2,1],markerscale=2,fontsize=11,\
           loc='upper right',frameon=False,borderaxespad=0.2,labelspacing=.2,handletextpad=0.2,handlelength=1)

plt.axhline(y=1,color='grey',lw=1,ls='--')
plt.axes().set_xticklabels([])
plt.axes().tick_params(axis='x',direction='out', length=0, width=0, colors='black')    
# plt.xlabel(xlabel)
# plt.title()
plt.ylabel('Gamma $k$')
plt.xlabel('Rank of TF')
plt.savefig('{}/gamma_alpha_by_num.pdf'.format(outdir),bbox_inches='tight',pad_inches=0.02,transparent=True)
plt.show()
plt.close()



# plot ranked by per condition
plt.figure(figsize=(6,3))
for ii in np.arange(len(columns)):
    x = np.arange(len(df.index))
    y = df[columns[ii]].sort_values(ascending=True)
    plt.scatter(x,y,c=colors[ii],s=s,alpha=alpha,label=xticklabels[ii])    

plt.legend( bbox_to_anchor=[1.2,1],markerscale=2,fontsize=11,\
            loc='upper right',frameon=False,borderaxespad=0.2,labelspacing=.2,handletextpad=0.2,handlelength=1)

plt.axhline(y=1,color='grey',lw=1,ls='--')
plt.axes().set_xticklabels([])
plt.axes().tick_params(axis='x',direction='out', length=0, width=0, colors='black')    
# plt.xlabel(xlabel)
# plt.title()
plt.ylabel('Gamma $k$')
plt.xlabel('Rank of TF')
plt.savefig('{}/gamma_alpha_order_per_num.pdf'.format(outdir),bbox_inches='tight',pad_inches=0.02,transparent=True)
plt.show()
plt.close()











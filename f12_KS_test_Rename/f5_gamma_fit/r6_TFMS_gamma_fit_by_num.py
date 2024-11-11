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




project_dir='/Volumes/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/'
project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/'
infiles = glob.glob('{}/f12_KS_test_Rename/f1_TF_cluster_potential/f0_bedtools_closest/data_TFMS_jarspar//*'.format(project_dir))
infiles = [i for i in infiles if not re.search('~',i)]
infiles.sort()
indir = '{}/f12_KS_test_Rename/f1_TF_cluster_potential/f0_bedtools_closest/data_TFMS_jarspar_by_num/'.format(project_dir)


outdir = 'f6_TFMS_gamma_alpha_by_num'
os.makedirs(outdir,exist_ok=True)
genome_len = 3298912062
num_thres = [2, 5, 10, 20, 50, 100]

out_df = pd.DataFrame()
for infile in infiles[:]:
    outname = os.path.basename(infile).split('.tsv')[0]
    print(outname)
    for num_thre in num_thres:
        try:
            df = pd.read_csv('{}/{}_top{}k.tsv'.format(indir,outname,num_thre),sep='\t',header=None)
            df = df[df[0].isin(hg38_chroms)]
            # len_of_motif = len(df.loc[0,3])
            col = df.columns[-1]
            df['dis'] = np.abs(df[col])
            values = df['dis']
            # values_log = np.log10(values.clip(1))
            alpha,loc,scale = gamma.fit(values,floc=0,fscale=genome_len/len(values))
            out_df.loc[outname,'alpha top{}k'.format(num_thre)] = alpha
            # out_df.loc[outname,'loc'] = loc
            out_df.loc[outname,'scale top{}k'.format(num_thre)] = scale
            out_df.loc[outname,'#TFMS top{}k'.format(num_thre)] = len(values)
        except:
            print('error',num_thre,outname) # in case of 0-line file

out_df.to_csv('{}/TFMS_gamma_alpha.csv'.format(outdir,))

df = pd.read_csv('f3_TFMS_CP_cor_SE_heatmap/TFMS_CP_SeEnrich_GammaFit.csv',index_col=0)
out_df = pd.concat([df,out_df],axis=1)
out_df.to_csv('{}/TFMS_gamma_alpha_combined.csv'.format(outdir,))





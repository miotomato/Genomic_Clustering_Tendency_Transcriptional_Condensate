import os,sys,argparse
# import fileinput,time
import glob
import re,bisect
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=14
import seaborn as sns
sns.set(font_scale=1.1)
sns.set_style("whitegrid", {'axes.grid' : False})
import scipy
from scipy import stats
import scipy.optimize
sns.set_style("ticks")
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
import statsmodels.stats.multitest as ssm



hg38_chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',
             'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',
             'chr18','chr19','chr20','chr21','chr22','chrX','chrY']


def compr_sites_distribution(df_cp,infile,outdir): 
    
    # == read the distance from the input
    outname = os.path.basename(infile).split('.tsv')[0]
    df = pd.read_csv(infile,sep='\t',header=None)
    df = df[df[0].isin(hg38_chroms)]
    len_of_motif = len(df.loc[0,3])
    col = df.columns[-1]
    df['dis'] = np.abs(df[col])
    values = df['dis']
    values_log = np.log10(values.clip(1))
    
    # == random distribute same number of sites as control
    # genome_len = 3298912062
    # np.random.seed(0)
    # random_loci = np.random.sample(df.shape[0])*genome_len
    # random_loci = random_loci.astype(int)
    # random_loci.sort()
    # dis_to_down = np.append(random_loci[1:],genome_len) - random_loci
    # dis_to_up = random_loci - np.append(0,random_loci[:-1]) 
    # dis_min = np.minimum(dis_to_down,dis_to_up)
    # dis_min_log = np.log10(dis_min.clip(1))
    
    # ==== exponential distribution as control
    nums = df.shape[0]
    mean_val = 3298912062/(nums+1)
    sample_vals = np.random.exponential(mean_val,nums)
    sample_vals_log = np.log10(sample_vals)

    # == plot the distribution of true/control data
    plt.figure(figsize=(2.6,2.6))
    sns.distplot(sample_vals_log, kde=True,color='tab:grey',label='Random',hist_kws={"lw": .0} )
    sns.distplot(values_log, kde=True,color='tab:green',label='TFMS',hist_kws={"lw": .0})
    plt.title(outname)
    plt.legend(fontsize=12,borderaxespad=0.1,labelspacing=.1,
               handletextpad=0.1,handlelength=1,
               # loc="upper left",
               frameon=False)
    plt.ylabel('Density of sites')
    plt.xlabel('Distance to closest site (log$_{10}$ bp)')
    plt.xlim([0,6])
    plt.savefig('{}/_fig/{}.pdf'.format(outdir,outname),bbox_inches='tight',pad_inches=0.02,transparent=True)
    # plt.show()
    plt.close()

    # == save the odds ratio between true and random data per data
    df_out = pd.DataFrame()
    for flag_distance in np.arange(20,4000,20):
        true_count = (values<flag_distance).sum()
        rand_count = (sample_vals<flag_distance).sum()
        odds_ratio = true_count/rand_count
        df_out.loc[flag_distance,'#true sites'] = true_count.astype(int)
        df_out.loc[flag_distance,'%true sites'] = (100*true_count/df.shape[0]).round(2)
        df_out.loc[flag_distance,'#random sites'] = rand_count.astype(int)
        df_out.loc[flag_distance,'%random sites'] = (100*rand_count/df.shape[0]).round(2)
        df_out.loc[flag_distance,'odds ratio'] = odds_ratio.round(2) 
    df_out.index.name = 'bp'
    df_out.to_csv('{}/_csv/{}_Odds_Ratio.csv'.format(outdir,outname))


    # == save the p-value from exponential distribution for each TFMS
    # quantile = -1*mean_val*np.log(1-p)
    lam = 1/mean_val
    df['pvalue'] = 1-np.exp(-1*lam*df['dis'])
    df['adj.p'] = ssm.fdrcorrection(df['pvalue'])[1]
    df = df[[0,1,2,'dis','pvalue','adj.p']]
    df.to_csv('{}/_csv/{}_Exponential_Pvalue.csv'.format(outdir,outname),index=False)
  
    # == save the cluster potential across all data
    df_cp.loc[outname,'#TFMS'] = df.shape[0]
    df_cp.loc[outname,'len-of-TFMS'] = len_of_motif
    df_cp.loc[outname,'#p<0.1'] =  df[df['pvalue']<.1].shape[0]
    df_cp.loc[outname,'#p<0.05'] = df[df['pvalue']<.05].shape[0]
    df_cp.loc[outname,'#p<0.01'] = df[df['pvalue']<.01].shape[0]
    df_cp.loc[outname,'#adj.p<0.1'] =  df[df['adj.p']<.1].shape[0]
    df_cp.loc[outname,'#adj.p<0.05'] = df[df['adj.p']<.05].shape[0]
    df_cp.loc[outname,'#adj.p<0.01'] = df[df['adj.p']<.01].shape[0]
    df_cp.loc[outname,'cutoff of OR>2'] = df_out[df_out['odds ratio']>2].index.max()
    df_cp.loc[outname,'true median'] = np.median(values)
    df_cp.loc[outname,'random median'] = np.median(sample_vals)
    # == compr the true vs. control
    s,p = stats.ttest_ind(sample_vals_log,values_log,)
    df_cp.loc[outname,'log10-dis T-test-s'] = s
    df_cp.loc[outname,'log10-dis T-test-p'] = p
    s,p = scipy.stats.ranksums(sample_vals_log,values_log,)
    df_cp.loc[outname,'log10-dis Wilcoxon-rank-sum-s'] = s
    df_cp.loc[outname,'log10-dis Wilcoxon-rank-sum-p'] = p
    s,p = scipy.stats.ks_2samp(sample_vals_log,values_log)
    df_cp.loc[outname,'log10-dis ks_2samp-s'] = s
    df_cp.loc[outname,'log10-dis ks_2samp-p'] = p
    s,p = scipy.stats.ks_2samp(sample_vals_log,values_log,alternative='less')
    df_cp.loc[outname,'log10-dis ks_2samp-s less'] = s
    df_cp.loc[outname,'log10-dis ks_2samp-p less'] = p
    s,p = scipy.stats.ks_2samp(sample_vals_log,values_log,alternative='greater')
    df_cp.loc[outname,'log10-dis ks_2samp-s greater'] = s
    df_cp.loc[outname,'log10-dis ks_2samp-p greater'] = p
    return df_cp

    

# ==== main
outdir = 'CP_TFMS_vs_random'
os.makedirs(outdir+os.sep+'_fig',exist_ok=True)
os.makedirs(outdir+os.sep+'_csv',exist_ok=True)

# == check the TF motifs
# infiles = glob.glob('../../../f9_TF_condensates_V3/f1_TF_cluster_potential/f0_bedtools_closest/data_TFMS_jarspar//*')
infiles = glob.glob('../../f1_TF_cluster_potential/f0_bedtools_closest/data_TFMS_jarspar//*')
infiles = [i for i in infiles if not re.search('~',i)]
infiles.sort()

df_cp = pd.DataFrame()
for infile in infiles[:]:
    try:
        # print(infile)
        df_cp = compr_sites_distribution(df_cp,infile,outdir)
    except:
        print('error',infile) # in case of 0-line file

df_cp.index.name = 'TF'
df_cp = df_cp.sort_values(by='log10-dis Wilcoxon-rank-sum-s',ascending=False)
df_cp.to_csv('{}/data_CP_TFMS_vs_random.csv'.format(outdir))




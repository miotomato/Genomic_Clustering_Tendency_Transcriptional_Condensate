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


hg38_chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',
             'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',
             'chr18','chr19','chr20','chr21','chr22','chrX','chrY']


def compr_sites_distribution(df_cp,infile,outdir): 
    
    # == read the distance from the input
    outname = os.path.basename(infile).split('.tsv')[0]
    df = pd.read_csv(infile,sep='\t',header=None)
    df = df[df[0].isin(hg38_chroms)]
    col = df.columns[-1]
    values = df[col]
    values_log = np.log10(values.clip(1))
    
    # == random distribute same number of sites as control
    genome_len = 3298912062
    np.random.seed(0)
    random_loci = np.random.sample(df.shape[0])*genome_len
    random_loci = random_loci.astype(int)
    random_loci.sort()
    dis_to_down = np.append(random_loci[1:],genome_len) - random_loci
    dis_to_up = random_loci - np.append(0,random_loci[:-1]) 
    dis_min = np.minimum(dis_to_down,dis_to_up)
    dis_min_log = np.log10(dis_min.clip(1))
    
    # == plot the distribution of true/control data
    plt.figure(figsize=(2.6,2.6))
    sns.distplot(dis_min_log,kde=True,color='tab:grey',label='Random',hist_kws={"lw": .0} )
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
        rand_count = (dis_min<flag_distance).sum()
        odds_ratio = true_count/rand_count
        df_out.loc[flag_distance,'#true sites'] = true_count.astype(int)
        df_out.loc[flag_distance,'%true sites'] = (100*true_count/df.shape[0]).round(2)
        df_out.loc[flag_distance,'#random sites'] = rand_count.astype(int)
        df_out.loc[flag_distance,'%random sites'] = (100*rand_count/df.shape[0]).round(2)
        df_out.loc[flag_distance,'odds ratio'] = odds_ratio.round(2) 
    df_out.index.name = 'bp'
    df_out.to_csv('{}/_csv/{}.csv'.format(outdir,outname))

    # == save the cluster potential across all data
    df_cp.loc[outname,'#TFMS'] = df.shape[0]
    df_cp.loc[outname,'cutoff of OR>2'] = df_out[df_out['odds ratio']>2].index.max()
    df_cp.loc[outname,'true median'] = np.median(values)
    df_cp.loc[outname,'random median'] = np.median(dis_min)
    # == compr the true vs. control
    s,p = stats.ttest_ind(dis_min_log,values_log,)
    df_cp.loc[outname,'log10 distance t-stats'] = s.round(2)
    df_cp.loc[outname,'log10 distance pvalue'] = '{:.2e}'.format(p) 
    return df_cp

    

# ==== main
outdir = 'CP_TFMS_vs_random'
os.makedirs(outdir+os.sep+'_fig',exist_ok=True)
os.makedirs(outdir+os.sep+'_csv',exist_ok=True)

# == check the TF motifs
infiles = glob.glob('../f1_bedtools_closest/data_TFMS_jarspar//*')
infiles = [i for i in infiles if not re.search('~',i)]
infiles.sort()

df_cp = pd.DataFrame()
for infile in infiles[:]:
    try:
        df_cp = compr_sites_distribution(df_cp,infile,outdir)
    except:
        print(infile) # in case of 0-line file

df_cp.index.name = 'TF'
df_cp = df_cp.sort_values(by='log10 distance t-stats',ascending=False)
df_cp.to_csv('{}/data_CP_TFMS_vs_random.csv'.format(outdir))




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
import scipy.optimize
sns.set_style("ticks")
matplotlib.rcParams["font.sans-serif"] = ["Arial"]




def compr_sites_distribution(infile,outdir,outname): 
    # == read the distance from the input
    df = pd.read_csv(infile,sep='\t',header=None)
    col = df.columns[-1]
    values = df[col]
    values_log = np.log10(values.clip(1))
    
    # == random distribute same number of sites as control
    genome_len=3*10**9
    np.random.seed(0)
    random_loci = np.random.sample(df.shape[0])*genome_len
    random_loci = random_loci.astype(int)
    random_loci.sort()
    dis_to_down = np.append(random_loci[1:],genome_len) - random_loci
    dis_to_up = random_loci - np.append(0,random_loci[:-1]) 
    dis_min = np.minimum(dis_to_down,dis_to_up)
    # dis_min = dis_min[dis_min>0]
    dis_min_log = np.log10(dis_min.clip(1))
    
    # == plot the distribution of true/control data
    plt.figure(figsize=(2.6,2.6))
    sns.distplot(dis_min_log,kde=True,color='tab:green',label='Random',hist_kws={"lw": .0} )
    sns.distplot(values_log, kde=True,color='tab:purple',label='True',hist_kws={"lw": .0})
    plt.title(outname)
    plt.legend(fontsize=12,borderaxespad=0.1,labelspacing=.1,
               handletextpad=0.1,handlelength=1,
               # loc="upper left",
               frameon=False)
    plt.ylabel('Density of sites')
    plt.xlabel('Distance to closest site (log$_{10}$ bp)')
    plt.xlim([0,6])
    plt.savefig('{}/{}.pdf'.format(outdir,outname),bbox_inches='tight',pad_inches=0.02,transparent=True)
    plt.show()
    plt.close()
    

    # == save the odds ratio between true and random data
    df_out = pd.DataFrame()
    for flag_distance in np.arange(20,10000,20):
        true_count = (values<flag_distance).sum()
        rand_count = (dis_min<flag_distance).sum()
        odds_ratio = true_count/rand_count
        df_out.loc[flag_distance,'log10 bp'] = np.log10(flag_distance).round(2) 
        df_out.loc[flag_distance,'#true sites'] = true_count.astype(int)
        df_out.loc[flag_distance,'%true sites'] = (100*true_count/df.shape[0]).round(2)
        df_out.loc[flag_distance,'#random sites'] = rand_count.astype(int)
        df_out.loc[flag_distance,'%random sites'] = (100*rand_count/df.shape[0]).round(2)
        df_out.loc[flag_distance,'odds ratio'] = odds_ratio.round(2) 
    df_out.index.name = 'bp'
    df_out.to_csv('{}/_csv/{}.csv'.format(outdir,outname))





outdir = 'f2_closest_distribution'
os.makedirs(outdir+os.sep+'_csv',exist_ok=True)


# == check the TF motifs
infiles=glob.glob('f1_bedtools_closest/*')
infiles.sort()
for infile in infiles:
    outname = os.path.basename(infile).split('.tsv')[0]
    try:
        compr_sites_distribution(infile,outdir,outname)
    except:
        print(outname) # in case of 0-line file




# df[df[col]<50].to_csv('test.bed',sep='\t',header=None,index=False)



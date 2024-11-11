import os,sys,argparse
# import fileinput,time
import glob
import re,bisect
import pandas as pd
import numpy as np
from scipy import stats
# import matplotlib
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt
# matplotlib.rcParams['font.size']=14
# import seaborn as sns
# sns.set(font_scale=1.1)
# sns.set_style("whitegrid", {'axes.grid' : False})
# import scipy
# import scipy.optimize
# sns.set_style("ticks")
# matplotlib.rcParams["font.sans-serif"] = ["Arial"]



# == check the TF motifs
outdir = 'f3_closest_distribution_stats'
os.makedirs(outdir,exist_ok=True)
infiles=glob.glob('f1_bedtools_closest/*')
infiles.sort()

df_out = pd.DataFrame()
for infile in infiles:
    outname = os.path.basename(infile).split('.tsv')[0]
    # print(outname)
    
    # == read the odds ratio vs. control
    or_file = 'f2_closest_distribution/_csv/{}.csv'.format(outname)
    if not os.path.isfile(or_file):
        print(or_file)
        continue
    or_df = pd.read_csv(or_file,index_col=0)
    or2_threshold = or_df[or_df['odds ratio']>2].index.max()
    or_of_100 = or_df.loc[100,'odds ratio']
    or_of_200 = or_df.loc[200,'odds ratio']

    # == read true values
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
    dis_min_log = np.log10(dis_min.clip(1))
    
    # == save the compare data
    df_out.loc[outname,'cutoff of OR>2'] = or2_threshold
    df_out.loc[outname,'OR of 100bp'] = or_of_100
    df_out.loc[outname,'OR of 200bp'] = or_of_200
    df_out.loc[outname,'true median'] = np.median(values)
    df_out.loc[outname,'random median'] = np.median(dis_min)
    # == compr the true vs. control
    s,p = stats.ttest_ind(dis_min_log,values_log,)
    df_out.loc[outname,'log10 distance t-stats'] = s.round(2)
    df_out.loc[outname,'log10 distance pvalue'] = '{:.2e}'.format(p) 
    # == raw distance    
    s,p = stats.ttest_ind(dis_min,values,)
    df_out.loc[outname,'distance t-stats'] = s.round(2)
    df_out.loc[outname,'distance pvalue'] = '{:.2e}'.format(p)


df_out.index.name = 'TF'
df_out.to_csv('{}/distance_to_closest_site_vs_control_stats.csv'.format(outdir))




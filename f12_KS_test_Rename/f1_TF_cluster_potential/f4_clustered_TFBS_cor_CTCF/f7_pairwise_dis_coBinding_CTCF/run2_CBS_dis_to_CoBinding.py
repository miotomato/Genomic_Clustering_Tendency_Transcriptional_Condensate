import os,sys,argparse
# import fileinput,time
import glob
import re,bisect
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=13
import seaborn as sns
sns.set(font_scale=1.0)
sns.set_style("whitegrid", {'axes.grid' : False})
import scipy
from scipy import stats
import scipy.optimize
sns.set_style("ticks")
matplotlib.rcParams["font.sans-serif"] = ["Arial"]


hg38_chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',
             'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',
             'chr18','chr19','chr20','chr21','chr22','chrX','chrY']


def get_lines_by_chrom(infile):
    df = pd.read_csv(infile,sep='\t',header=None,low_memory=False)
    df = df[df[0].isin(hg38_chroms)]
    return df.shape[0]


def bedtools_sample(infile,number,outdir,outname):    
    sample_file = '{}/tmp/{}.sampled.bed'.format(outdir,outname)
    sort_file = '{}/tmp/{}.sampled.sort.bed'.format(outdir,outname)
    # ==== bedtools sample, sort and closest
    commandLine = 'bedtools sample -n {} -i {} > {}\n'.format(number,infile,sample_file)
    os.system(commandLine);print(commandLine)   
    commandLine = 'bedtools sort -i {} > {}\n'.format(sample_file,sort_file)
    os.system(commandLine);print(commandLine)  
    # commandLine = 'rm {} \n'.format(sample_file)
    # os.system(commandLine);print(commandLine)  
    return sort_file


def return_distance_distribution(afile,bfile,outdir,outname,flag):
    closest_file = '{}/tmp/{}_{}.closest.bed'.format(outdir,outname,flag)
    commandLine = 'bedtools closest -a {} -b {} -D ref -fd -io -t first > {}\n'.format(afile,bfile,closest_file)
    os.system(commandLine);print(commandLine)  

    df = pd.read_csv(closest_file,sep='\t',header=None,low_memory=False)
    df = df[df[0].isin(hg38_chroms)]
    col = df.columns[-1]
    values = np.abs(df[col])
    values_log = np.log10(values.clip(1))
    return values,values_log

    
def compr_pairwise_distance(df,peak_file,treat_file,control_file,outdir,outname,labels,colors): 
    # == read the distance from the input fiels
    values_treat,values_log_treat = return_distance_distribution(peak_file,treat_file,outdir,outname,'treat')
    total_sites = get_lines_by_chrom(treat_file)
    # ==== sample from control file and do bedtools closest
    sort_file = bedtools_sample(control_file,total_sites,outdir,outname)
    values_control,values_log_control = return_distance_distribution(peak_file,sort_file,outdir,outname,'control')
    
    if int(outname.split('sample')[1]) ==0:
        # == plot the distribution of true/control data
        plt.figure(figsize=(2.6,2.6))
        sns.distplot(values_log_control, kde=True,color=colors[1],label=labels[1],hist_kws={"lw": .0})
        sns.distplot(values_log_treat,kde=True,color=colors[0],label=labels[0],hist_kws={"lw": .0} )
        plt.title(outname)
        plt.legend(fontsize=10,borderaxespad=0.1,labelspacing=.1,
                    handletextpad=0.1,handlelength=1,
                    loc="upper left",
                    frameon=False)
        plt.ylabel('Density of sites')
        plt.xlabel('Distance to closest Peak (log$_{10}$ bp)')
        # plt.xlim([0,9])
        plt.savefig('{}/{}.pdf'.format(outdir,outname),bbox_inches='tight',pad_inches=0.02,transparent=True)
        # plt.show()
        plt.close()
        
    # == save the cluster potential across all data
    df.loc[outname,'#{}'.format(labels[0])] = total_sites
    df.loc[outname,'#{}'.format(labels[1])] = total_sites
    df.loc[outname,'{} dis-to-neighbor median'.format(labels[0])] = np.median(values_treat)
    df.loc[outname,'{} dis-to-neighbor median'.format(labels[1])] = np.median(values_control)
    # df.loc[outname,'{} percentile-1'.format(labels[1])] = np.percentile(values_control,1).round(2)
    # df.loc[outname,'{} percentile-5'.format(labels[1])] = np.percentile(values_control,5).round(2)
    # == compr the true vs. control
    s,p = stats.ttest_ind(values_log_control,values_log_treat)
    df.loc[outname,'log10-dis T-test-s'] = s
    df.loc[outname,'log10-dis T-test-p'] = p 
    # remove temp closest file
    # commandLine = 'rm {} \n'.format(closest_file)
    # os.system(commandLine);print(commandLine)  
    return df



# ==== main
outdir = 'f2_CBS_dis_to_coBinding'
os.makedirs(outdir+os.sep+'tmp',exist_ok=True)

cancertypes = ['COAD','BRCA']
df_out_summary = pd.DataFrame()

for cancertype in cancertypes:
    # get a/b files for bedtools closest    
    peak_file = 'data/{}_gained.sort.bed'.format(cancertype)
    treat_file = 'data/percentile_T_{}_coBinding.sort.bed'.format(cancertype)
    control_file = 'data/hg38_unionDHS_fc4_50merge.sort.bed'
    # ==== outcome info
    labels = ['Co-binding'.format(cancertype),'Sampled UDHS']   
    colors = ['tab:red','tab:grey'] 
    df_out = pd.DataFrame()
    for ii in np.arange(10): 
        outname_ii = '{}_sample{}'.format(cancertype,ii)   
        df_out = compr_pairwise_distance(df_out,peak_file,treat_file,control_file,outdir,outname_ii,labels,colors)
    # df_out.index.name = 'TF'
    # df_cp = df_cp.sort_values(by='log10 distance t-stats',ascending=False)
    df_out.to_csv('{}/{}.csv'.format(outdir,cancertype))
    df_out_summary = pd.concat([df_out_summary,df_out.median()],axis=1)

df_out_summary.columns = cancertypes
df_out_summary.to_csv('{}/summary.csv'.format(outdir))
    



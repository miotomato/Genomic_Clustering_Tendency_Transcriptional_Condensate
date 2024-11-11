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
sns.set(font_scale=1.0)
sns.set_style("whitegrid", {'axes.grid' : False})
import scipy
from scipy import stats
import scipy.optimize
sns.set_style("ticks")
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
import subprocess


hg38_chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',
             'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',
             'chr18','chr19','chr20','chr21','chr22','chrX','chrY']

def return_distance_distribution(infile):
    df = pd.read_csv(infile,sep='\t',header=None,low_memory=False)
    df = df[df[0].isin(hg38_chroms)]
    col = df.columns[-1]
    values = np.abs(df[col])
    values_log = np.log10(values.clip(1))
    return values,values_log


def bedtools_sample_closest(infile,number,outdir,outname):    
    sample_file = '{}/{}.sampled.bed'.format(outdir,outname)
    sort_file = '{}/{}.sampled.sort.bed'.format(outdir,outname)
    closest_file = '{}/{}.sampled.closest.bed'.format(outdir,outname)
    # ==== bedtools sample, sort and closest
    commandLine = 'bedtools sample -n {} -i {} > {}\n'.format(number,infile,sample_file)
    os.system(commandLine);print(commandLine)   
    commandLine = 'bedtools sort -i {} > {}\n'.format(sample_file,sort_file)
    os.system(commandLine);print(commandLine)  
    commandLine = 'bedtools closest -a {} -b {} -D ref -fd -io -t first > {}\n'.format(sort_file,sort_file,closest_file)
    os.system(commandLine);print(commandLine)  
    commandLine = 'rm {} \n'.format(sample_file)
    os.system(commandLine);print(commandLine)  
    return sort_file,closest_file
    
    
def compr_cluster_potential(df,treat_file,control_file,outdir,outname,labels,colors): 
    # == read the distance from the input fiels
    values_treat,values_log_treat = return_distance_distribution(treat_file)
    total_sites = len(values_treat) 
    # ==== sample from control file and do bedtools closest
    sort_file,closest_file = bedtools_sample_closest(control_file,total_sites,outdir,outname)
    values_control,values_log_control = return_distance_distribution(closest_file)
    
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
        plt.xlabel('Distance to closest site (log$_{10}$ bp)')
        # plt.xlim([0,9])
        plt.savefig('{}/_fig/{}.pdf'.format(outdir,outname),bbox_inches='tight',pad_inches=0.02,transparent=True)
        # plt.show()
        plt.close()
        
        # == save the odds ratio between true and random data per data
        df_out = pd.DataFrame()
        for flag_distance in np.arange(20,10000,20):
            true_count = (values_treat<flag_distance).sum()
            ctrl_count = (values_control<flag_distance).sum()
            odds_ratio = true_count/max(ctrl_count,1)
            df_out.loc[flag_distance,'#true sites'] = true_count.astype(int)
            df_out.loc[flag_distance,'%true sites'] = (100*true_count/total_sites).round(2)
            df_out.loc[flag_distance,'#control sites'] = ctrl_count.astype(int)
            df_out.loc[flag_distance,'%control sites'] = (100*ctrl_count/total_sites).round(2)
            df_out.loc[flag_distance,'odds ratio'] = round(odds_ratio,2)
        df_out.index.name = 'bp'
        df_out.to_csv('{}/_csv_sample/{}.csv'.format(outdir,outname))

    # == save the cluster potential across all data
    df.loc[outname,'#{}'.format(labels[0])] = total_sites
    df.loc[outname,'#{}'.format(labels[1])] = total_sites
    # df.loc[outname,'cutoff of OR>2'] = df_out[df_out['odds ratio']>2].index.max()
    df.loc[outname,'{} dis-to-neighbor median'.format(labels[0])] = np.median(values_treat)
    df.loc[outname,'{} dis-to-neighbor median'.format(labels[1])] = np.median(values_control)
    df.loc[outname,'{} percentile-1'.format(labels[1])] = np.percentile(values_control,1).round(2)
    df.loc[outname,'{} percentile-5'.format(labels[1])] = np.percentile(values_control,5).round(2)
    # == compr the true vs. control
    s,p = stats.ttest_ind(values_log_control,values_log_treat)
    df.loc[outname,'log10-dis T-test-s'] = s
    df.loc[outname,'log10-dis T-test-p'] = p
    s,p = stats.ranksums(values_log_control,values_log_treat)
    df.loc[outname,'log10-dis Wilcoxon-rank-sum-s'] = s
    df.loc[outname,'log10-dis Wilcoxon-rank-sum-p'] = p
    s,p = stats.ks_2samp(values_log_control,values_log_treat)
    df.loc[outname,'log10-dis ks_2samp-s'] = s
    df.loc[outname,'log10-dis ks_2samp-p'] = p
    # remove temp closest file
    commandLine = 'rm {} \n'.format(closest_file)
    os.system(commandLine);print(commandLine)  
    return df,sort_file

    

def get_lines(infile):
    with open(infile,'rb') as f:
        lines = 0
        buf_size = 1024*1024
        buf = f.raw.read(buf_size)
        while buf:
            lines += buf.count(b'\n')
            buf = f.raw.read(buf_size)
    return lines

def get_lines_by_chrom(infile):
    df = pd.read_csv(infile,sep='\t',header=None,low_memory=False)
    df = df[df[0].isin(hg38_chroms)]
    return df.shape[0]


def run_bedtools_intersect(afile,bfile):
    cl = 'bedtools intersect -a {} -b {} -u -wa|wc -l'.format(afile,bfile)
    print(cl)
    wa_count = subprocess.check_output(cl,shell=True).decode(sys.stdout.encoding).strip()
    return int(wa_count)

    
def enrichment_odds_ratio(df,treat_file,control_file,base_file,outname,labels):
    # == whether the motif overlap with SE
    treat_total = get_lines_by_chrom(treat_file)
    treat_overlap = run_bedtools_intersect(treat_file,base_file)
    # == whether the cluster motif overlap with SE
    control_total = get_lines_by_chrom(control_file)
    control_overlap = run_bedtools_intersect(control_file,base_file)
    # == fisher exact test
    a = treat_overlap
    b = treat_total-treat_overlap
    c = control_overlap
    d = control_total - control_overlap
    s,p = stats.fisher_exact([[a, b], [c, d]])
    
    df.loc[outname,'{} overlap SE'.format(labels[0])] = treat_overlap
    df.loc[outname,'{} overlap SE'.format(labels[1])] = control_overlap
    # df.loc[outname,'{} total'.format(labels[0])] = treat_total
    # df.loc[outname,'{} total'.format(labels[1])] = control_total
    df.loc[outname,'fisher_exact_s'] = round(s,2)
    df.loc[outname,'fisher_exact_p'] = '{:.2e}'.format(p) 
    return df






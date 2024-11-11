import os,sys,argparse
# import fileinput,time
import glob
import re,bisect
import pandas as pd
import numpy as np
import matplotlib
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt
# matplotlib.rcParams['font.size']=14
# import seaborn as sns
# sns.set(font_scale=1.1)
# sns.set_style("whitegrid", {'axes.grid' : False})
# import scipy
# from scipy import stats
# import scipy.optimize
# sns.set_style("ticks")
# matplotlib.rcParams["font.sans-serif"] = ["Arial"]
# import statsmodels.stats.multitest as ssm
import subprocess
from scipy import stats



hg38_chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',
             'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',
             'chr18','chr19','chr20','chr21','chr22','chrX','chrY']



def run_bedtools_intersect(df,bfile,indir,prename):
    
    afile = '{}/{}.bed'.format(indir,prename)
    df.to_csv(afile,index=False,sep='\t',header=None)
    cl = 'bedtools intersect -a {} -b {} -u -wa|wc -l'.format(afile,bfile)
    # print(cl)
    wa_count = subprocess.check_output(cl,shell=True).decode(sys.stdout.encoding).strip()
    os.system('rm {}'.format(afile))
    return int(wa_count)
    

# ==== main
indir = 'CP_TFMS_vs_random'
outdir = 'CP_TFMS_vs_random_SE_enrichment'
os.makedirs(outdir,exist_ok=True)

infiles = glob.glob('{}/_csv/*Exponential*'.format(indir))
infiles = [i for i in infiles if not re.search('~',i)]
infiles.sort()

se_files = glob.glob('../../data/SE_hg38/*.bed')
p_thres = [0.05,0.01]

for se_file in se_files[:]:
    se_prename = os.path.basename(se_file).split('.bed')[0]
    df_out = pd.DataFrame()
    for infile in infiles[:]:
        prename = os.path.basename(infile).split('_Exponential_Pvalue')[0]
        # ==== total TFMS and those overlapped with SE
        df = pd.read_csv(infile,)
        total = df.shape[0]
        total_overlapped = run_bedtools_intersect(df,se_file,indir,prename)
        df_out.loc[prename,'total'] = total
        df_out.loc[prename,'SE overlapped'] = total_overlapped
    
        # ==== TFMS with p cutoff and those overlapped with SE
        for p_thre in p_thres:
            df_p = df[df['pvalue']<p_thre]
            sub = df_p.shape[0]
            sub_overlapped = run_bedtools_intersect(df_p,se_file,indir,prename)
            # print(total,total_overlapped,sub,sub_overlapped)
            s,p = stats.fisher_exact([[sub_overlapped, sub - sub_overlapped], 
                                      [total_overlapped - sub_overlapped, total - sub - total_overlapped + sub_overlapped]])
                
            df_out.loc[prename,'TFMS-p<{} total'.format(p_thre)] = sub
            df_out.loc[prename,'TFMS-p<{} SE overlapped'.format(p_thre)] = sub_overlapped
            df_out.loc[prename,'TFMS-p<{} enrich-at-SE-fisher-exact-s'.format(p_thre)] = round(s,2)
            df_out.loc[prename,'TFMS-p<{} enrich-at-SE-fisher-exact-p'.format(p_thre)] = '{:.2e}'.format(p) 
        
    
    df_out.to_csv('{}/data_TFMS_enrich_at_SE_{}.csv'.format(outdir,se_prename))
    
    


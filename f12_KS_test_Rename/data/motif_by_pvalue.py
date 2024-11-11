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



infiles = glob.glob('motif_fimo_jaspar/*')
infiles = [i for i in infiles if not re.search('~',i)]
infiles.sort()


outdir = 'motif_fimo_jaspar_by_pvalue'
os.makedirs(outdir,exist_ok=True)
# genome_len = 3298912062
p_thres = [5,6,7]


for infile in infiles[:]:
    outname = os.path.basename(infile).split('.bed')[0]
    print(outname)
    try:   
    # if 1:
        df = pd.read_csv(infile,sep='\t',header=None)
        df = df[df[0].isin(hg38_chroms)]
        for p_thre in p_thres:
            df_tmp = df[df[4]<10**(-1*p_thre)]
            df_tmp.to_csv('{}/{}_p{}.bed'.format(outdir,outname,p_thre),sep='\t',index=False,header=None)
            
    except:
        print('error',outname) # in case of 0-line file



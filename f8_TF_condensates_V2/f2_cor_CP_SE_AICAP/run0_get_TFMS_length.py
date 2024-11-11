import os,sys,argparse
# import fileinput,time
import glob
import re,bisect
import pandas as pd
import numpy as np
# import matplotlib
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


hg38_chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',
             'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',
             'chr18','chr19','chr20','chr21','chr22','chrX','chrY']


   

# ==== main
outdir = 'data'
os.makedirs(outdir,exist_ok=True)

# == check the TF motifs
# file_dir="/Volumes/zanglab/shared/Motif/sites/hg38_fimo_jarspar/results/"
file_dir="/nv/vol190/zanglab/shared/Motif/sites/hg38_fimo_jarspar/results/"
infiles = glob.glob('{}/*'.format(file_dir))
infiles = [i for i in infiles if not re.search('~',i)]
infiles.sort()

df_cp = pd.DataFrame()
for infile in infiles[:]:
    try:
        df = pd.read_csv(infile,sep='\t',header=None)
        outname = os.path.basename(infile).split('.bed')[0]
        df_cp.loc[outname,'motif_Length'] = len(df[3][0])
    except:
        print(infile) # in case of 0-line file

df_cp.index.name = 'TF'
# df_cp = df_cp.sort_values(by='log10 distance t-stats',ascending=False)
df_cp.to_csv('{}/data_motif_Length.csv'.format(outdir))




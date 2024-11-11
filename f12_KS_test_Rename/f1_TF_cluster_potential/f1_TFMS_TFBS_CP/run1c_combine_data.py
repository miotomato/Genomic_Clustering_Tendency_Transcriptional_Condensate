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




# ==== main

outdir = 'CP_TFMS_vs_random'
infiles = glob.glob('{}/_csv/*enrich_at_SE*'.format(outdir))
infiles = [i for i in infiles if not re.search('~',i)]
infiles.sort()


df_out = pd.DataFrame()
for infile in infiles[:]:
    basename = os.path.basename(infile).split('_enrich_at_SE')[0]
    df = pd.read_csv(infile,index_col=0)  
    df_out = pd.concat([df_out,df.median().rename(basename)],axis=1)

df_out = np.transpose(df_out)
df_out = df_out[['motif total', 'motif SE overlapped','random total', 'random SE overlapped', 
                  'enrich-at-SE-fisher-exact-s',
        'enrich-at-SE-fisher-exact-p', ]]
df_out.to_csv('{}/data_TFMS_enrich_at_SE.csv'.format(outdir))



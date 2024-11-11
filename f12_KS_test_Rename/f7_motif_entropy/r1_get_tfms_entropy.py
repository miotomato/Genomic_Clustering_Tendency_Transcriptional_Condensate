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
from pyjaspar import jaspardb



# project_dir='/Volumes/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/'

infiles = glob.glob('../data/motif_fimo_jaspar/*')
infiles = [i for i in infiles if not re.search('~',i)]
infiles.sort()

jdb_obj = jaspardb()
print(jdb_obj.release)

outdir = 'f1_TFMS_entropy'
os.makedirs(outdir,exist_ok=True)
# genome_len = 3298912062

out_df = pd.DataFrame()
for infile in infiles[:]:
    outname = os.path.basename(infile).split('.bed')[0]
    motif_id = outname.split('_')[1]
    name = outname.split('_')[0]
    try:
        # jaspar_file = 'data/JASPAR2022_CORE_redundant_pfms_jaspar/{}.jaspar'.format(motif_id)
        motif = jdb_obj.fetch_motif_by_id(motif_id)
        pfm = pd.DataFrame.from_dict(motif.counts)
        pwm = pfm.divide(pfm.sum(axis=1),axis=0)
        entropy = 2+(np.log2(pwm)*pwm).sum(axis=1)
        out_df.loc[outname,'entropy_sum'] = entropy.sum()
        out_df.loc[outname,'entropy_avg'] = entropy.mean()
        out_df.loc[outname,'motif_len'] = len(entropy)
    except:
        print(outname)
        # motifs = jdb_obj.fetch_motifs_by_name(name)
        # pfm = pd.DataFrame.from_dict(motifs[0].counts)
        # pwm = pfm.divide(pfm.sum(axis=1),axis=0)
        # entropy = 2+(np.log2(pwm)*pwm).sum(axis=1)
        # out_df.loc[outname,'entropy_sum'] = entropy.sum()
        # out_df.loc[outname,'entropy_avg'] = entropy.mean()
        # out_df.loc[outname,'motif_len'] = len(entropy)

out_df.to_csv('{}/TFMS_entropy.csv'.format(outdir,))






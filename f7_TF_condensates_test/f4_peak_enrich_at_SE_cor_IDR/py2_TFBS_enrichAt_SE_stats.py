import os,sys,argparse
# import fileinput,time
import glob
import re,bisect
import pandas as pd
import numpy as np
import subprocess
from scipy import stats


    


# ==== main
    
outdir = 'f2_TFBS_enrichAt_SE_stats'
os.makedirs(outdir,exist_ok=True)

df = pd.read_csv('f1_slurm_log/cistrome_id_log.csv',index_col=0)   
permutation_types = ['shuffleGenome','sampleUDHS','sampleMergePeak']

for celltype in df.columns[:]:
    out_df = pd.DataFrame()
    factor_df = df[celltype].dropna()
    for factor in factor_df.index[:]:
        cistrome_id = factor_df[factor].astype(int)
        permutation_file = 'f1_TFBS_enrich_at_SE/{}_{}_{}.csv'.format(factor,celltype,cistrome_id)
        permutation_df = pd.read_csv(permutation_file,index_col=0)
        for permutation_type in permutation_types:
            permutation_data = permutation_df['SE_overlapped_{}'.format(permutation_type)]
            true_data = permutation_data.peak_file_overlapped
            sample_data = permutation_data.values[2:]
            relative_rank = (sample_data>=true_data).sum()
            zscore = (true_data-np.mean(sample_data))/np.std(sample_data)
            pvalue = stats.norm.sf(np.abs(zscore))
            out_df.loc[factor,'{}_relative_rank'.format(permutation_type)] = relative_rank
            out_df.loc[factor,'{}_zscore'.format(permutation_type)] = zscore
            out_df.loc[factor,'{}_pvalue'.format(permutation_type)] = pvalue

    out_df.to_csv('{}/{}_stats.csv'.format(outdir,celltype))


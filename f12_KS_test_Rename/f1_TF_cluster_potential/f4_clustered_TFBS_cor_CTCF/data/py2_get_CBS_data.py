import os,sys,argparse
import glob
import re,bisect
import pandas as pd
import numpy as np
# from GenomeData import *




indir = '/Volumes/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/fz_union_combination/f2_union_binding_data/'
outdir = './'

# get the combined df and mid positions
combined_df = pd.read_csv('{}/union_binding.csv'.format(indir))
combined_df.iloc[:,:4].to_csv('{}/union_CBS.bed'.format(outdir),index=False,sep='\t',header=None)


# high quality bindings (occupancy score ≥3)
occupancy3_filtered = combined_df[combined_df['occupancy_score']>=3]
# occupancy3_filtered.to_csv('{}/union_CBS_occupancy_GT3.csv'.format(outdir),index=False)
occupancy3_filtered.iloc[:,:4].to_csv('{}/union_CBS_occupancy_GT3.bed'.format(outdir),index=False,sep='\t',header=None)


# constitutive bindings

occupancy_thre = 0.80
kept_df = combined_df[combined_df['occupancy_score']>=771*occupancy_thre]
# kept_df.to_csv('{}/union_CBS_constitutive_{}.csv'.format(outdir,occupancy_thre),index=False)
kept_df.iloc[:,:4].to_csv('{}/union_CBS_thre_{}.bed'.format(outdir,occupancy_thre),index=False,sep='\t',header=None)


for occupancy_thre in [0.1,0.2,0.3]:
    kept_df = combined_df[combined_df['occupancy_score']<=771*occupancy_thre]
    # kept_df.to_csv('{}/union_CBS_specificity_{}.csv'.format(outdir,occupancy_thre),index=False)
    kept_df.iloc[:,:4].to_csv('{}/union_CBS_thre_{}.bed'.format(outdir,occupancy_thre),index=False,sep='\t',header=None)

    kept_df = occupancy3_filtered[occupancy3_filtered['occupancy_score']<=771*occupancy_thre]
    kept_df.iloc[:,:4].to_csv('{}/union_CBS_occupancy_GT3_thre_{}.bed'.format(outdir,occupancy_thre),index=False,sep='\t',header=None)
    

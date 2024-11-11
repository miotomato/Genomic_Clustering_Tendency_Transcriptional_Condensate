import os,sys,argparse
# import fileinput,time
import glob
import re,bisect
import pandas as pd
import numpy as np




outdir = 'f5_motif_within_300bp'
os.makedirs(outdir,exist_ok=True)


# == check the TF motifs
infiles=glob.glob('f1_bedtools_closest/*')
infiles.sort()
for infile in infiles:
    outname = os.path.basename(infile).split('.tsv')[0]
    try:
        df = pd.read_csv(infile,sep='\t',header=None)
        col = df.columns[-1]
        df = df[df[col]<=300].iloc[:,:3]
        df.to_csv('{}/{}.bed'.format(outdir,outname),sep='\t',header=None,index=False)
    except:
        print(outname) # in case of 0-line file







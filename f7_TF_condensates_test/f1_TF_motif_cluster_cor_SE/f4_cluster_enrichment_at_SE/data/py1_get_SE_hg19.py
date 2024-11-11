import os,sys,argparse
# import fileinput,time
import glob
import re,bisect
import pandas as pd
import numpy as np




outdir = 'SE_hg19'
os.makedirs(outdir,exist_ok=True)


# == check the TF motifs
infiles=glob.glob('1-s2.0-S0092867413012270-mmc7/*')
infiles.sort()

for infile in infiles[:]:
    outname = os.path.basename(infile).split('.csv')[0]
    outname = re.split('^UCSD_|BI_',outname)[-1]
    df = pd.read_csv(infile,skiprows = 1, header=None,index_col=0)
    df.columns = ['chromosome', 'start', 'end','associated gene', 'enhancer rank', 
                  'is enhancer','H3K27ac density (rpm/bp)','read density']
    # == select SE
    df = df[df['is enhancer']==1]
    df.index.name='SE_id'
    df.to_csv('{}/{}.csv'.format(outdir,outname))
    bed_col = ['chromosome', 'start', 'end']
    df[bed_col].to_csv('{}/{}.bed'.format(outdir,outname),header=None,index=False,sep='\t')
    print(outname) # in case of 0-line file







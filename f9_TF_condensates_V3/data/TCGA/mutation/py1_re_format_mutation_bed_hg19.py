import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect

#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#matplotlib.rcParams['font.size']=16
#matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
#matplotlib.rcParams["font.family"] = "sans-serif"
#import seaborn as sns
#sns.set(font_scale=2)
#sns.set_style("whitegrid", {'axes.grid' : False})
#sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k'})
#sns.despine(offset=0, trim=True)

from collections import Counter






indir = 'download'
outdir = 'icgc_mutation_bed_hg19'
os.makedirs(outdir,exist_ok=True)

gene_df = pd.read_csv('/Volumes/zanglab/zw5j/env/HOMER/data/accession/human2gene.tsv',sep='\t',header=None)
gene_df = gene_df[[4,6]].dropna()
gene_df.columns = ['ensembl','symbol']
gene_df.index = gene_df.ensembl
# gene_df = gene_df.drop_duplicates()
gene_df = gene_df[~gene_df.index.duplicated(keep='first')]


for cancertype in ['BRCA', 'CESC', 'COAD', 'LIHC', 'PRAD'][:]:
    # for each cancertype, keep the chr pos ref/alt mutation_rate for each mutation
    # mutation_rate_file='icgc_mutation_rate/{}_mutation_rate.txt'.format(cancertype)
    # mutation_rate_df = pd.read_csv(mutation_rate_file,index_col=0,sep='\t')
    # mutation_rate_df['mutation_rate']=mutation_rate_df['#mutation_donor']/mutation_rate_df['#all_donor']#;print(mutation_rate_df);exit()
    infile = '{}/simple_somatic_mutation.open.{}-US.tsv'.format(indir,cancertype)
    outfile = '{}/{}_mutation.bed'.format(outdir,cancertype)
    inf = open(infile)
    outf = open(outfile,'w')
    #outf.write('{}\t{}\t{}\t{}\t{}\n'.format('chr','pos','ref','alt','genotype'))
    line = inf.readline()
    line = inf.readline()# start from the second line
    collected_mutation_ids = set()
    
    ii=0
    while line :
        sline = line.split('\t')
        icgc_mutation_id = sline[0]#;print(icgc_mutation_id)
        if icgc_mutation_id not in collected_mutation_ids:
            chr = sline[8]
            pos = int(sline[9])
            submitted_sample_id = sline[6]
            TCGA_id = '-'.join(submitted_sample_id.split('-')[:3])
            gene_affected = sline[28]
            try:
                gene_symbol = gene_df.loc[gene_affected].symbol
            except:
                gene_symbol = 'None'
            # print(gene_symbol)
            ref = sline[14]
            mutation_from = sline[15]
            alt = sline[16]
            mutation_type = '{}>{}'.format(ref,alt)
            #print(chr,pos,ref,alt,genotype);exit()
            # outf.write('chr{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chr,pos,pos+1,TCGA_id,mutation_type,'+',gene_affected,gene_symbol,))
            outf.write('chr{}\t{}\t{}\t{}\t{}\t{}\n'.format(chr,pos,pos+1,TCGA_id,gene_symbol,'+'))
            collected_mutation_ids.add(icgc_mutation_id)
            # ii+=1
            
        line = inf.readline()
        
    outf.close()
    inf.close()
    print(cancertype)
    #exit()


  

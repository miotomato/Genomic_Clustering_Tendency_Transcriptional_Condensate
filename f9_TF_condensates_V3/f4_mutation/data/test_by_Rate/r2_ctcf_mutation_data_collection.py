import sys,argparse
import os,glob
import numpy as np
import pandas as pd
from scipy import stats
import re,bisect

chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',\
             'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',\
             'chr18','chr19','chr20','chr21','chr22','chrX','chrY']




outdf = pd.DataFrame()


outfile = 'ctcf_constitutive_mutatioinCount.bed'
df = pd.read_csv(outfile,sep='\t')
# df.columns = ['chr','start','end','id','count']
total_len = (df['end']-df['start']).sum()
total_mutation_count = df['mutationCount'].sum()
total_mutation_rate = df['mutationRate'].sum()*6285
# df['mutation rate per bp'] = df['count']/(df['end']-df['start'])
outdf.loc['CTCF','total_len'] = total_len
outdf.loc['CTCF','total mutationCount'] = total_mutation_count
outdf.loc['CTCF','%mutationCount per bp'] = total_mutation_count/total_len
outdf.loc['CTCF','total mutationRate'] = total_mutation_rate
outdf.loc['CTCF','%mutationRate per bp'] = total_mutation_rate/total_len
                
outdf.to_csv('CTCF_murationRate.csv')










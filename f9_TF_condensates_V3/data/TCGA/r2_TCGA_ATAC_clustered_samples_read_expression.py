import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *
from scipy import stats
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#matplotlib.rcParams['font.size']=16
#matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
#matplotlib.rcParams["font.family"] = "sans-serif"
#import seaborn as sns
#sns.set(font_scale=2)
#sns.set_style("whitegrid", {'axes.grid' : False})

import re,bisect
#plus = re.compile('\+')
#minus = re.compile('\-')



# == clustered samples for survival analysis
filtered_df = pd.read_excel('TCGA-ATAC_clustered_samples.xlsx',index_col=1)   


# == ATAC identifier
identifier = pd.read_csv('TCGA_identifier_mapping.txt',sep='\t',index_col=3)
identifier = identifier[['Case_ID']].drop_duplicates(keep='first')
identifier['Case_ID_short'] = [i[:16] for i in identifier['Case_ID']]
kept_caseID = identifier.loc[filtered_df.index].Case_ID_short


# == expression data
expr_file='GDC-PANCAN.htseq_fpkm-uq.tsv'
expr_df = pd.read_csv(expr_file,sep='\t',index_col=0)
kept_expr_col = expr_df.columns.intersection(kept_caseID)
uuid_col = [identifier[identifier.Case_ID_short==col].index[0] for col in kept_expr_col]
expr_df = expr_df[kept_expr_col]
expr_df.columns = uuid_col


# === change to gene symbol
geneID_df = pd.read_csv('BRCA_1_DEseq2Results_geneSymbol.csv',index_col=1)
filtered_Ensemble = geneID_df.geneSymbol.dropna()
expr_df = expr_df.loc[expr_df.index.intersection(filtered_Ensemble.index)]
expr_df.index = filtered_Ensemble.loc[expr_df.index].values
expr_df.to_csv('TCGA-ATAC_clustered_samples.htseq_fpkm-uq.tsv',sep='\t')






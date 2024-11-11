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
# filtered_df = pd.read_excel('TCGA-ATAC_clustered_samples.xlsx',index_col=1)   


# == ATAC identifier
# identifier = pd.read_csv('TCGA_identifier_mapping.txt',sep='\t',index_col=3)
# identifier = identifier[['Case_ID']].drop_duplicates(keep='first')
# identifier['Case_ID_short'] = [i[:16] for i in identifier['Case_ID']]
# kept_caseID = identifier.loc[filtered_df.index].Case_ID_short


# == expression data
expr_file='GDC-PANCAN.htseq_fpkm-uq.tsv'
# expr_file='GDC-PANCAN.htseq_fpkm-uq.head100.tsv'
expr_df = pd.read_csv(expr_file,sep='\t',index_col=0)

# == Tumor types range from 01 - 09, normal types from 10 - 19 and control samples from 20 - 29. 
kept_expr_col = [i for i in expr_df.columns if i.split('-')[-1].startswith('0')]
expr_df = expr_df[kept_expr_col]
expr_df.columns = [i[:-4] for i in expr_df.columns]

# === change to gene symbol
geneID_df = pd.read_csv('BRCA_1_DEseq2Results_geneSymbol.csv',index_col=1)
filtered_Ensemble = geneID_df.geneSymbol.dropna()
expr_df = expr_df.loc[expr_df.index.intersection(filtered_Ensemble.index)]
expr_df.index = filtered_Ensemble.loc[expr_df.index].values

# ==== remove duplicated columns and index
duplicated_col = expr_df.columns[expr_df.columns.duplicated()]
expr_df = expr_df.drop(duplicated_col,axis=1)
expr_df = expr_df.loc[~expr_df.index.duplicated(keep='first')]


for cancertype in ['BRCA','COAD'][:]:    
    clinical_data = '2022-03-20'  
    infile = 'clinical.project-TCGA-{}.{}.csv'.format(cancertype,clinical_data)
    clinical_df = pd.read_csv(infile,index_col=1)
    # ==== samples with both expr and clinical info
    cancer_samples = clinical_df.index.intersection(expr_df.columns) 
    print(cancertype,len(cancer_samples))
    expr_df[cancer_samples].to_csv('GDC-PANCAN.htseq_fpkm-uq.{}.tsv'.format(cancertype),sep='\t')
    
    
# x=expr_df.iloc[:,33]
# y=expr_df.iloc[:,33]
# plt.scatter(x,y)


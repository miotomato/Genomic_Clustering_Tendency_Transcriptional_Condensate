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



# ==== main
indir='f1_get_RP_from_bedGraph'
outdir='f2_avg_RP_per_sample'
os.makedirs(outdir,exist_ok=True)

# atac_id_info = pd.read_csv('../../data/TCGA/TCGA_identifier_mapping.txt',sep='\t',index_col=3)
filtered_df = pd.read_excel('../../data/TCGA/TCGA-ATAC_clustered_samples.xlsx',index_col=0)   
data_types = ['tcga_atac','HCT-116.merge','MCF-7.merge']

for data_type in data_types[:]:
    for cancertype in ['BRCA','COAD'][:]:    
        # ==== averaged RP from all replicates per cancer per data type
        df_cancertype = pd.DataFrame()
        clinical_df = pd.read_csv('../../data/TCGA//clinical.project-TCGA-{}.2022-03-20.csv'.format(cancertype),index_col=1)
        filtered_ids = filtered_df[filtered_df.cohort==cancertype].index
    
        # ==== average RP from all replicates per patient
        for filtered_id in filtered_ids[:]:
            df_id = pd.DataFrame()
            infiles = glob.glob('{}/{}_*_{}_*_{}.csv'.format(indir,cancertype,filtered_id,data_type))
            for infile in infiles:
                df = pd.read_csv(infile,sep='\t',index_col = 0)
                df_id = pd.concat([df_id,df],axis=1)
            df_id = df_id.mean(axis=1).rename(filtered_id)
            print(data_type, cancertype, filtered_id) 
            # ==== average RP per patient 
            df_cancertype = pd.concat([df_cancertype,df_id],axis=1)
        
        df_cancertype.to_csv('{}/{}_avg_RP_on_{}.csv'.format(outdir,cancertype,data_type))    
        
        

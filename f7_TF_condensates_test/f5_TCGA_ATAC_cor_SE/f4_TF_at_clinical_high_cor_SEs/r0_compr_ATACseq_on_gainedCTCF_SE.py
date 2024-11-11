import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *

# import matplotlib
# # matplotlib.use('Agg')
# import matplotlib.pyplot as plt
# matplotlib.rcParams['font.size']=18
# import seaborn as sns
# sns.set(font_scale=1.2)

# sns.set_style("whitegrid", {'axes.grid' : False})
# sns.set_style('ticks')
# import re,bisect
# from lifelines.statistics import logrank_test
# from lifelines import KaplanMeierFitter
# matplotlib.rcParams["font.sans-serif"] = ["Arial"]




# ==== main 
outdir = 'f0_ATACseq_on_gained_CTCF_and_SE'
os.makedirs(outdir,exist_ok=True)
pardir = '/Volumes/zanglab/zw5j/'

cancer_cell_line_match = {'BRCA':'MCF-7','COAD':'HCT-116'}
outdf = pd.DataFrame()

for cancertype in ['BRCA','COAD'][:]:
    gained_ctcf_overlapped_peak_file = '{}/work2017/T_ALL_CTCF/updated_201906/f5_TCGA_ATAC/f4_clinical_survival_panCancer/f1_ATAC_sig_caseID_differential_score/{}_gained_caseID_ATAC_sig_diff_score.csv'.format(pardir,cancertype)
    se_overlapped_peak_file = '{}/since2019_projects/phase_separation_FEpiTR/f7_TF_condensates_test/f5_TCGA_ATAC_cor_SE/f2_TCGA_ATAC_at_SE_cor_survival/f1_ATAC_overlap_SE_caseID/{}_ATAC_overlap_{}_SE_caseID.bed'.format(pardir,cancertype,cancer_cell_line_match[cancertype])
    gained_overlapped_peaks = pd.read_csv(gained_ctcf_overlapped_peak_file,sep='\t',index_col = 3)
    se_overlapped_peaks = pd.read_csv(se_overlapped_peak_file,sep='\t',index_col = 3)
    
    outdf.loc[cancertype,'# gained CTCF overlapped ATACseq peaks'] = gained_overlapped_peaks.shape[0]
    outdf.loc[cancertype,'# SE overlapped ATACseq peaks'] = se_overlapped_peaks.shape[0]
    overlapped_index = gained_overlapped_peaks.index.intersection(se_overlapped_peaks.index)
    outdf.loc[cancertype,'# shared ATACseq peaks'] = len(overlapped_index)
    outdf.loc[cancertype,'ID shared ATACseq peaks'] = ', '.join(overlapped_index)
    
outdf.to_csv(outdir+os.sep+'summary.csv')
import os,sys,argparse,glob,re,bisect
import numpy as np
import pandas as pd
from collections import Counter
import operator
import matplotlib
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt
# #matplotlib.rcParams['agg.path.chunksize'] = 10000
# matplotlib.rcParams['font.size']=16
# import seaborn as sns
# sns.set(font_scale=1.2)
# sns.set_style("whitegrid", {'axes.grid' : False})
# import scipy
# import scipy.optimize
# sns.set_style("ticks")
# matplotlib.rcParams["font.sans-serif"] = ["Arial"]

from scipy.stats import gamma


project_dir='/Volumes/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/f12_KS_test_Rename'
# project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/'

writer = pd.ExcelWriter('data/CTFBS_clinical_annotation.xlsx')
infile = '{}/f1_TF_cluster_potential/f3_clustered_TFBS/f6_atac_overlap_coBinding_TFBS_figs/HCT-116_top_zscored_TFBSCP/_COAD_percentile_T_ExtendMerge_JUND-CEBPB-SRF_logRank.csv'.format(project_dir)
df = pd.read_csv(infile,index_col=0)
df = df.drop(['fdr'],axis=1)
df = df.sort_values(by='log rank p',ascending=True)
df.to_excel(writer,'COAD_coBinding_overlapped_Peak')

infile = '{}/f1_TF_cluster_potential/f3_clustered_TFBS/f6_atac_overlap_coBinding_TFBS_figs/HCT-116_top_zscored_TFBSCP/COAD_percentile_T_ExtendMerge_clinical.csv'.format(project_dir)
df = pd.read_csv(infile,index_col=0)
df.to_excel(writer,'COAD_data_summary')

writer.close()
                



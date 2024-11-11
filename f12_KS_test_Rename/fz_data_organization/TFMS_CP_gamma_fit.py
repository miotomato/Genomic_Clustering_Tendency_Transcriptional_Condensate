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

# ==== read TFMS CP
infile = '{}/f1_TF_cluster_potential/f2_cor_CP_SE_AICAP/data/TFMS_CP_SE_enrich.csv'.format(project_dir)
df1 = pd.read_csv('../f1_TF_cluster_potential/f2_cor_CP_SE_AICAP/data//TFMS_CP_SE_enrich.csv',index_col=0)
df1 = df1[['#TFMS','len-of-TFMS','log10-dis ks_2samp-s signed','motif SE overlapped','enrich-at-SE-fisher-exact-s', 'enrich-at-SE-fisher-exact-p']]
df1.columns = ['#TFMS','length of TFMS','TFMS CP','#TFMS overlap SE','SE enrichment odds ratio', 'SE enrichment pvalue']

# ==== read Gamma k by p-value
infile = '{}/f5_gamma_fit/f4_TFMS_gamma_alpha_by_pvalue/TFMS_gamma_alpha_combined.csv'.format(project_dir)
df2 = pd.read_csv(infile,index_col=0)
# df2 = df2[['alpha','scale',
#            '#TFMS p5','alpha p5','scale p5',
#            '#TFMS p6','alpha p6','scale p6',
#            '#TFMS p7','alpha p7','scale p7',]]
# df2.columns = ['Gamma k','Gamma theta',
#            '#TFMS p<1e-5','Gamma k p<1e-5','Gamma theta p<1e-5',
#            '#TFMS p<1e-6','Gamma k p<1e-6','Gamma theta p<1e-6',
#            '#TFMS p<1e-7','Gamma k p<1e-7','Gamma theta p<1e-7',]

df2 = df2[['alpha','scale',
           'alpha p5','alpha p6','alpha p7']]
df2.columns = ['Gamma k','Gamma theta',
           'Gamma k p<1e-5','Gamma k p<1e-6','Gamma k p<1e-7',]
      
# ==== read Gamma k by number
infile = '{}/f5_gamma_fit/f6_TFMS_gamma_alpha_by_num/TFMS_gamma_alpha_combined.csv'.format(project_dir)
df3 = pd.read_csv(infile,index_col=0)
# df3 = df3[['#TFMS top2k','alpha top2k','scale top2k',
#            '#TFMS top5k','alpha top5k','scale top5k',
#            '#TFMS top10k','alpha top10k','scale top10k',
#            '#TFMS top20k','alpha top20k','scale top20k',
#            '#TFMS top50k','alpha top50k','scale top50k',
#            '#TFMS top100k','alpha top100k','scale top100k',]]
 
# df3.columns = ['#TFMS top2k','Gamma k top2k','Gamma theta top2k',
#                '#TFMS top5k','Gamma k top5k','Gamma theta top5k',
#                '#TFMS top10k','Gamma k top10k','Gamma theta top10k',
#                '#TFMS top20k','Gamma k top20k','Gamma theta top20k',
#                '#TFMS top50k','Gamma k top50k','Gamma theta top50k',
#                '#TFMS top100k','Gamma k top100k','Gamma theta top100k',] 

df3 = df3[['alpha top2k','alpha top5k','alpha top10k',
           'alpha top20k','alpha top50k','alpha top100k']]
 
df3.columns = ['Gamma k top2k','Gamma k top5k','Gamma k top10k',
               'Gamma k top20k','Gamma k top50k','Gamma k top100k'] 

# ==== combine all data
df = pd.concat([df1,df2,df3],axis=1)
df = df.sort_values(by='TFMS CP',ascending=False)

writer = pd.ExcelWriter('data/TFMS_CP_Gamma_fit.xlsx')
df.to_excel(writer)
writer.close()
                



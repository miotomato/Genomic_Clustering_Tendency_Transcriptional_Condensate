import os,sys,argparse
# import fileinput,time
import glob
import re,bisect
import pandas as pd
import numpy as np
import matplotlib
from scipy import stats

    

# ==== main
indirs = ['CP_TFBS_top10k_vs_TFMS','CP_TFBS_overlap_motif_vs_TFMS','CP_TFBS_NOT_overlap_motif_vs_TFMS']

for indir in indirs[:]:
    infiles = glob.glob('{}/_csv_cp/*csv'.format(indir))
    infiles.sort()
    
    # ==== combine statistics/p-values per data
    df_out = pd.DataFrame()
    for infile in infiles[:]:
        basename = os.path.basename(infile).split('.csv')[0]
        df = pd.read_csv(infile,index_col=0)
        df.replace([np.inf, -np.inf], np.nan, inplace=True)
        df = df.dropna()
        # ==== save the median values
        df_out[basename] = df.median()
        # for pcol in ['log10-distance CP pvalue','fisher_exact_p']:
        #     freedom = 2*df.shape[0]
        #     sum_score = (-2*np.log10(df[pcol])).sum()
        #     fisher_p = stats.chi2.sf(sum_score,freedom)
        #     df_out.loc['{} Fisher-data-combined'.format(pcol),basename] = fisher_p
    df_out = np.transpose(df_out)  
    df_out.to_csv('{}/per_data_{}.csv'.format(indir,indir))
    
    
    # ==== combined per TF
    df_out2 = pd.DataFrame(columns = df_out.columns)
    tfs = set([i.split('_')[1] for i in df_out.index])
    for tf in sorted(tfs):
        data_index = [i for i in df_out.index if re.search('_{}_'.format(tf),i)]
        df = df_out.loc[data_index]; # print(tf,df.shape)
        df_out2.loc[tf] = df.median()
        for pcol in ['log10-dis T-test-p',
                      'log10-dis Wilcoxon-rank-sum-p',
                      'fisher_exact_p',
                      'log10-dis ks_2samp-p',
                      ]:
            freedom = 2*df.shape[0]
            sum_score = (-2*np.log10(df[pcol])).sum()
            fisher_p = stats.chi2.sf(sum_score,freedom)
            df_out2.loc[tf,'{} Fisher-TF-combined'.format(pcol)] = fisher_p
        
    df_out2.to_csv('{}/per_TF_{}.csv'.format(indir,indir))
        

    # ==== combined per TF per celltype
    df_out3 = pd.DataFrame(columns = df_out.columns)
    tfs = set([i.split('_')[1] for i in df_out.index])
    celltypes = set([i.split('_')[0] for i in df_out.index])
    for tf in sorted(tfs):
        for celltype in celltypes:
            row_name = '{}_{}'.format(celltype,tf)
            data_index = [i for i in df_out.index if re.search('{}_{}_'.format(celltype,tf,),i)]
            # print(celltype,tf,data_index)
            df = df_out.loc[data_index]; # print(row_name,df.shape)
            if df.shape[0] ==0:
                continue
            df_out3.loc[row_name] = df.median()
            for pcol in ['log10-dis T-test-p',
                          'log10-dis Wilcoxon-rank-sum-p',
                          'fisher_exact_p',
                          'log10-dis ks_2samp-p',
                          ]:
                freedom = 2*df.shape[0]
                sum_score = (-2*np.log10(df['{}'.format(pcol)])).sum()
                fisher_p = stats.chi2.sf(sum_score,freedom)
                df_out3.loc[row_name,'{} Fisher-TF-Celltype-combined'.format(pcol)] = fisher_p
        
    df_out3.to_csv('{}/per_TF_per_Celltype_{}.csv'.format(indir,indir))

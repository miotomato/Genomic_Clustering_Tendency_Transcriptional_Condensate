import os,sys,argparse
# import fileinput,time
import glob
import re,bisect
import pandas as pd
import numpy as np
import matplotlib
from scipy import stats

    
tfbs_prenames = ['TFBS_all','TFBS_overlap_motif','TFBS_NOT_overlap_motif']
percentile_thres = ['percentile-5','percentile-1']

for percentile_thre in percentile_thres[:1]:
    outdir = 'f1_clustered_TFBS_{}'.format(percentile_thre)
    os.makedirs(outdir,exist_ok=True)
    
    for tfbs_prename in tfbs_prenames[:1]:
        cp_prename = 'CP_{}_vs_TFMS'.format(tfbs_prename)
        df = pd.read_csv('../f1_TFMS_TFBS_CP/{}/per_data_{}.csv'.format(cp_prename,cp_prename),index_col=0)
        df = df.dropna(how='all')
        
        # df = df[(df['#TFBS']>2000)&(df['#TFBS']<100000)]
        df = df[(df['#TFBS_all']>2000)]
        
        # ==== get the clustered TFBS
        for index in df.index[:]:
            data_id = index.split('_')[2]
            tf = index.split('_')[1]
            ct = index.split('_')[0]
            if ct not in ['MCF-7','HCT-116','HeLa','LNCaP','U87','HepG2']:
                continue
            # if tf not in ['ERG','E2F1','MYC','SRF','JUND','CEBPB','FOSL1','MAX']:
            #     continue
            # print(data_id)
            # ==== read the percentile threshold
            dis_thre = df.loc[index,'TFMS {}'.format(percentile_thre)].astype(int)
            # data_file = '../f0_bedtools_closest/data_{}/{}_sort_peaks.narrowPeak.tsv'.format(tfbs_prename,data_id)
            tfbs_prename_old = 'TFBS' if tfbs_prename == 'TFBS_all' else tfbs_prename
            data_file = '../../../f8_TF_condensates_V2/f1_TF_cluster_potential/f1_bedtools_closest//data_{}/{}_sort_peaks.narrowPeak.tsv'.format(tfbs_prename_old,data_id)
            if not os.path.isfile(data_file):
                print(index)
                # continue
            
            data_df = pd.read_csv(data_file,sep='\t',header=None)
            dis_col = data_df.columns[-1]
            data_df_T = data_df[data_df[dis_col]<=dis_thre]
            data_df_C = data_df[data_df[dis_col]>dis_thre]
            
            # ==== save the bed files of clustered TFBS
            os.makedirs(outdir+os.sep+ct,exist_ok=True)
            outname = '{}_{}_{}'.format(ct,tf,data_id)
            se_file = '../../data/SE_hg38/{}.bed'.format(ct)
            bed_T_file = '{}/{}/{}_T.bed'.format(outdir,ct,outname)
            bed_T_SE = '{}/{}/{}_T_on_SE.bed'.format(outdir,ct,outname)
            bed_C_file = '{}/{}/{}_C.bed'.format(outdir,ct,outname)
            bed_C_SE = '{}/{}/{}_C_on_SE.bed'.format(outdir,ct,outname)
            merge_file = '{}/{}/{}_T_ExtendMerge.bed'.format(outdir,ct,outname)
            merge_SE = '{}/{}/{}_T_ExtendMerge_on_SE.bed'.format(outdir,ct,outname)
            data_df_T.to_csv(bed_T_file,sep='\t',header=None,index=False)
            data_df_C.to_csv(bed_C_file,sep='\t',header=None,index=False)
            # commandLine = '''cat {}|awk '{{OFS="\\t";print$1,$2,$3+{}}}'>{}'''.format(bed_file,dis_thre.astype(int),bed_file_extend)
            # print(bed_T_file);print(bed_C_file)
            commandLine = 'bedtools merge -d {} -i {} > {}'.format(dis_thre,bed_T_file,merge_file)
            print(commandLine);os.system(commandLine)
            commandLine = 'bedtools intersect -a {} -b {} -wa -u > {}'.format(bed_T_file,se_file,bed_T_SE)
            print(commandLine);os.system(commandLine)
            commandLine = 'bedtools intersect -a {} -b {} -wa -u > {}'.format(bed_C_file,se_file,bed_C_SE)
            print(commandLine);os.system(commandLine)
            commandLine = 'bedtools intersect -a {} -b {} -wa -u > {}\n'.format(merge_file,se_file,merge_SE)
            print(commandLine);os.system(commandLine)
            df.loc[index,'#TFBS_new'] = data_df.shape[0]
    
        df = df.dropna()
        df.to_csv('{}/{}_new.csv'.format(outdir,cp_prename))
    








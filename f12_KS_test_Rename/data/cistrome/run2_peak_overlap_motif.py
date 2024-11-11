import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect



hg38_chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',
             'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',
             'chr18','chr19','chr20','chr21','chr22','chrX','chrY']


outdir = 'peaks_all'
outdir1 = 'peaks_overlap_motif'
outdir2 = 'peaks_NOT_overlap_motif'
os.makedirs(outdir,exist_ok=True)
os.makedirs(outdir1,exist_ok=True)
os.makedirs(outdir2,exist_ok=True)


# ==== keep high quality data
peak_dir = '/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor'
# peak_dir = '/Volumes/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor'
df_qc = pd.read_csv('cistrome2019_selected_QC.csv',index_col=0)           

# ==== read peak data with matched motif and SE 
motif_dir="/nv/vol190/zanglab/shared/Motif/sites/hg38_fimo_jarspar/results/"
# motif_dir="/Volumes/zanglab/shared/Motif/sites/hg38_fimo_jarspar/results/"
df_count = pd.read_csv('cistrome2019_data_Count_with_motif.csv',index_col=0,keep_default_na=False)

top_thre = None
# ==== get peaks with/without overlap with motif
for factor in df_count.index[:]:
    motif_files = glob.glob('{}/{}_*.bed'.format(motif_dir,factor))
    motif_files = [i for i in motif_files if not re.search('~',i)]
    # print(factor,motif_files)
    try:
        assert len(motif_files)==1
        for celltype in df_count.columns[:]:
            if df_count.loc[factor,celltype]==0:
                continue
            df_tmp = df_qc[(df_qc.Factor==factor)&(df_qc.Cell_line==celltype)]
            for dcID in df_tmp.index[:]:
                motif_file = motif_files[0]
                # select the top 10 peaks for each dataset
                peak_file = '{}/{}_sort_peaks.narrowPeak.bed'.format(peak_dir,dcID)
                # print(motif_file)
                # print(peak_file)
                peak_df = pd.read_csv(peak_file,sep='\t',header=None)
                peak_df = peak_df[peak_df[0].isin(hg38_chroms)][:top_thre]
                output_bed = '{}_{}_{}.bed'.format(celltype,factor,dcID)
                peak_df.to_csv('{}/{}'.format(outdir,output_bed),sep='\t',header=None,index=False)
                commandLine = 'bedtools intersect -a {}/{} -b {} -wa -u > {}/{}'.format(outdir,output_bed,motif_file,outdir1,output_bed)
                os.system(commandLine);print(commandLine)
                commandLine = 'bedtools intersect -a {}/{} -b {} -wa -v > {}/{}\n'.format(outdir,output_bed,motif_file,outdir2,output_bed)
                os.system(commandLine);print(commandLine)
    except:
        print('===ERROR===',factor,motif_files)        
                


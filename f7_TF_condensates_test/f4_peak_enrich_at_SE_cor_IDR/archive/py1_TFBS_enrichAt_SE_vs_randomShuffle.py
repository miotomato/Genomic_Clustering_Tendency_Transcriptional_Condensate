import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect
# import matplotlib
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt
# matplotlib.rcParams['font.size']=16
# import seaborn as sns
# sns.set(font_scale=1.2)
# sns.set_style("whitegrid", {'axes.grid' : False})
# sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k'})
# matplotlib.rcParams["font.sans-serif"] = ["Arial"]
# matplotlib.rcParams['mathtext.fontset'] = 'custom'
# matplotlib.rcParams["mathtext.rm"] = "Arial"
import subprocess

def get_lines(infile):
    with open(infile,'rb') as f:
        lines = 0
        buf_size = 1024*1024
        buf = f.raw.read(buf_size)
        while buf:
            lines += buf.count(b'\n')
            buf = f.raw.read(buf_size)
    return lines


def run_bedtools_intersect_A_on_B(afile,bfile):
    # module load gcc/9.2.0
    # module load bedtools/2.29.2
    cl = 'bedtools intersect -a {} -b {} -u -wa|wc -l'.format(afile,bfile)
    # print(cl)
    wa_count = subprocess.check_output(cl,shell=True).decode(sys.stdout.encoding).strip()
    return int(wa_count)


def return_cistrome_peak_file(cistrome_id):
    cistrome_dir = '/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new'
    # cistrome_dir = '/Volumes/zanglab/zw5j/data/cistrome/cistrome_db_2019_new'
    for batch in ['human_hm','human_factor']:
        peak_file='{}/{}/{}_sort_peaks.narrowPeak.bed'.format(cistrome_dir,batch,cistrome_id)
        if os.path.isfile(peak_file):
            return peak_file


outdir = 'f1_TFBS_enrich_at_SE_vs_shuffleGenome'
os.makedirs(outdir,exist_ok=True)

cistrome_id = 45963
data_df = pd.read_csv('../f3_CellType_specific_peak_cor_SE/f1_cellType_Factor_data_Count/cistrome2019_data_details.csv',index_col=0)   
factor = data_df.loc[cistrome_id,'Factor']
celltype = data_df.loc[cistrome_id,'Cell_line']
peak_file = return_cistrome_peak_file(cistrome_id)
celltype_SE_file = '../f1_TF_motif_cluster_cor_SE/f4_cluster_enrichment_at_SE/data/SE_hg38/{}.bed'.format(celltype)
genome = '/nv/vol190/zanglab/zw5j/data/Genome/ucsc/hg38/hg38_clean.chrom.sizes'
udhs_file = '/nv/vol190/zanglab/zw5j/data/unionDHS/hg38_unionDHS_fc4_50merge.bed'
  
# peak file overlapping SE
df = pd.DataFrame()
all_num = get_lines(peak_file)
overlapped_num = run_bedtools_intersect_A_on_B(peak_file,celltype_SE_file)
df.loc['peak_file_all','SE_overlapped'] = all_num
df.loc['peak_file_overlapped','SE_overlapped'] = overlapped_num
# random select background regions
for ii in np.arange(20):
    shuffle_file = '{}/shuffle.bed'.format(outdir)
    # commandline = 'bedtools shuffle -i {} -g {} -incl {} > {}'.format(peak_file,genome,udhs_file,shuffle_file)
    commandline = 'bedtools shuffle -i {} -g {} > {}'.format(peak_file,genome,shuffle_file)
    os.system(commandline)
    # print(commandline)
    overlapped_num = run_bedtools_intersect_A_on_B(shuffle_file,celltype_SE_file)
    df.loc['shuffle_{}'.format(ii),'SE_overlapped'] = overlapped_num

df.to_csv(outdir+os.sep+'{}_{}_{}.csv'.format(factor,celltype,cistrome_id))   
os.system('rm {}'.format(shuffle_file))




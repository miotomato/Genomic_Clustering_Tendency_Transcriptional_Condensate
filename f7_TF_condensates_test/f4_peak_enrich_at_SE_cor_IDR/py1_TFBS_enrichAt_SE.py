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
    cistrome_dir = '/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019'
    for batch in ['human_hm','human_factor']:
        df = pd.read_csv('{}/{}.txt'.format(cistrome_dir,batch),sep='\t',index_col=0)
        if cistrome_id in df.index:
            peak_file = '{}/{}/{}'.format(cistrome_dir,batch,df.loc[cistrome_id,'File'])
            return peak_file



def main(cistrome_id,factor,celltype):
    
    outdir = 'f1_TFBS_enrich_at_SE'
    os.makedirs(outdir,exist_ok=True)
    
    # read data
    peak_file = return_cistrome_peak_file(cistrome_id) 
    celltype_SE_file = '../../data/SE_hg38/{}.bed'.format(celltype)
    genome = '/nv/vol190/zanglab/zw5j/data/Genome/ucsc/hg38/hg38_clean.chrom.sizes'
    udhs_file = '/nv/vol190/zanglab/zw5j/data/unionDHS/hg38_unionDHS_fc4_50merge.bed'
    merge_file = '../../data/cistrome_data/{}.merge.bed'.format(factor)
    tmp_file = '{}/tmp_{}_{}_{}.bed'.format(outdir,factor,celltype,cistrome_id)
    
    # peak file overlapping SE
    columns = ['SE_overlapped_shuffleGenome',
               'SE_overlapped_sampleUDHS',
               'SE_overlapped_sampleMergePeak']
    df = pd.DataFrame(columns = columns)
    # overlapp between original data and SE
    all_num = get_lines(peak_file)
    overlapped_num = run_bedtools_intersect_A_on_B(peak_file,celltype_SE_file)
    df.loc['peak_file_all'] = all_num
    df.loc['peak_file_overlapped'] = overlapped_num
    
    # random select background regions
    for ii in np.arange(1000):
        # == random select background regions from genome
        column = 'SE_overlapped_shuffleGenome'
        commandline = 'bedtools shuffle -i {} -g {} > {}'.format(peak_file,genome,tmp_file)
        os.system(commandline)
        if ii==0:
            print(commandline)
        overlapped_num = run_bedtools_intersect_A_on_B(tmp_file,celltype_SE_file)
        df.loc['random_{}'.format(ii),column] = overlapped_num
    
        # == random select background regions from UDHS
        column = 'SE_overlapped_sampleUDHS'
        commandline = 'bedtools sample -i {} -n {} > {}'.format(udhs_file,all_num,tmp_file)
        os.system(commandline)
        if ii==0:
            print(commandline)
        overlapped_num = run_bedtools_intersect_A_on_B(tmp_file,celltype_SE_file)
        df.loc['random_{}'.format(ii),column] = overlapped_num

        # == random select background regions from merged peaks
        column = 'SE_overlapped_sampleMergePeak'
        commandline = 'bedtools sample -i {} -n {} > {}'.format(merge_file,all_num,tmp_file)
        os.system(commandline)
        if ii==0:
            print(commandline)
        overlapped_num = run_bedtools_intersect_A_on_B(tmp_file,celltype_SE_file)
        df.loc['random_{}'.format(ii),column] = overlapped_num
        # print(ii)
    
    
    df.to_csv(outdir+os.sep+'{}_{}_{}.csv'.format(factor,celltype,cistrome_id))   
    os.system('rm {}'.format(tmp_file))





if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', action = 'store', type = int, dest = 'cistrome_id')
    parser.add_argument('-f', action = 'store', type = str, dest = 'factor')
    parser.add_argument('-c', action = 'store', type = str, dest = 'celltype')
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main(args.cistrome_id,args.factor,args.celltype)

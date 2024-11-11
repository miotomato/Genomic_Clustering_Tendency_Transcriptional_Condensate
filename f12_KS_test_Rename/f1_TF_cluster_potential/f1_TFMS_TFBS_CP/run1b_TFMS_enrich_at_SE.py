import os,sys,argparse
# import fileinput,time
import glob
import re,bisect
import pandas as pd
import numpy as np
import matplotlib
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt
# matplotlib.rcParams['font.size']=14
# import seaborn as sns
# sns.set(font_scale=1.1)
# sns.set_style("whitegrid", {'axes.grid' : False})
# import scipy
# from scipy import stats
# import scipy.optimize
# sns.set_style("ticks")
# matplotlib.rcParams["font.sans-serif"] = ["Arial"]
# import statsmodels.stats.multitest as ssm
import subprocess
from scipy import stats



hg38_chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',
             'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',
             'chr18','chr19','chr20','chr21','chr22','chrX','chrY']




def run_bedtools_intersect(df,bfile,outdir,prename):
    
    motif_file = '{}/{}.bed'.format(outdir,prename)
    df.to_csv(motif_file,index=False,sep='\t',header=None)
    cl = 'bedtools intersect -a {} -b {} -u -wa|wc -l'.format(motif_file,bfile)
    # print(cl)
    wa_count = subprocess.check_output(cl,shell=True).decode(sys.stdout.encoding).strip()
    # os.system('rm {}'.format(afile))
    return int(wa_count),motif_file
    
    

def run_bedtools_shuffle_intersect(ii,motif_file,se_file,outdir,prename):

    motif_shuffle = '{}/{}_shuffle_{}.bed'.format(outdir,prename,ii)
    # motif_se_overlapped = '{}/{}_se_overlappe_{}.bed'.format(outdir,prename,ii)
    genome_file = '/nv/vol190/zanglab/zw5j/data/Genome/ucsc/hg38/hg38_clean.chrom.sizes'
    
    cl = 'bedtools shuffle -i {} -g {} > {}'.format(motif_file,genome_file,motif_shuffle)
    print(cl);os.system(cl)
    cl = 'bedtools intersect -a {} -b {} -u -wa|wc -l'.format(motif_shuffle,se_file)
    print(cl);os.system(cl)
    wa_count = subprocess.check_output(cl,shell=True).decode(sys.stdout.encoding).strip()
    
#     os.system('rm {}'.format(motif_file))
    os.system('rm {}'.format(motif_shuffle))
#     os.system('rm {}'.format(motif_se_overlapped))
    
    return int(wa_count)
    




# ==== main
# outdir = 'CP_TFMS_vs_random'
# 
# se_file = '../../data/SE_hg38/all_hg38_SE.bed'
# infiles = glob.glob('{}/_csv/*Exponential*'.format(outdir))
# infiles = [i for i in infiles if not re.search('~',i)]
# infiles.sort()

# p_thres = [0.05,0.01]
# df_out = pd.DataFrame()


# for infile in infiles[:3]:

def main(infile):

    outdir = 'CP_TFMS_vs_random'
    se_file = '../../data/SE_hg38/all_hg38_SE.bed'
    
    df_out_tmp = pd.DataFrame()
    
    prename = os.path.basename(infile).split('_Exponential_Pvalue')[0]
    # ==== total TFMS and those overlapped with SE
    df = pd.read_csv(infile,)
    total = df.shape[0]
    total_overlapped,motif_file = run_bedtools_intersect(df,se_file,outdir,prename)
    
    # ==== TFMS with p cutoff and those overlapped with SE
    for ii in np.arange(100):
        random_overlapped = run_bedtools_shuffle_intersect(ii,motif_file,se_file,outdir,prename)
#         df_p = df[df['pvalue']<p_thre]
#         sub = df_p.shape[0]
#         sub_overlapped = run_bedtools_intersect(df_p,se_file,outdir,prename)
        # print(total,total_overlapped,sub,sub_overlapped)
        s,p = stats.fisher_exact([[total_overlapped, total - total_overlapped], 
                                  [random_overlapped, total - random_overlapped]])
            
        df_out_tmp.loc[ii,'random total'] = total
        df_out_tmp.loc[ii,'random SE overlapped'] = random_overlapped
        df_out_tmp.loc[ii,'enrich-at-SE-fisher-exact-s'] = round(s,2)
        df_out_tmp.loc[ii,'enrich-at-SE-fisher-exact-p'] = '{:.2e}'.format(p) 

    os.system('rm {}'.format(motif_file))
    
    df_out_tmp['motif total'] = total
    df_out_tmp['motif SE overlapped'] = total_overlapped
    df_out_tmp.to_csv('{}/_csv/{}_enrich_at_SE.csv'.format(outdir,prename))    

# df_out.loc[prename] = df_out_tmp.median()
# df_out.to_csv('{}/data_TFMS_enrich_at_SE.csv'.format(outdir))
# 



if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'infile', metavar = '<file>')
    
    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main(args.infile)



import os,sys,argparse
# import fileinput,time
import glob
import re,bisect
import pandas as pd
import numpy as np
import subprocess
from scipy import stats


def get_lines(infile):
    with open(infile,'rb') as f:
        lines = 0
        buf_size = 1024*1024
        buf = f.raw.read(buf_size)
        while buf:
            lines += buf.count(b'\n')
            buf = f.raw.read(buf_size)
    return lines


def run_bedtools_intersect(afile,bfile):
    cl = 'bedtools intersect -a {} -b {} -u -wa|wc -l'.format(afile,bfile)
    # print(cl)
    wa_count = subprocess.check_output(cl,shell=True).decode(sys.stdout.encoding).strip()
    return int(wa_count)
    



def main(se_basename):
    
    outdir = 'f1_cluster_enrichment_at_SE'
    os.makedirs(outdir,exist_ok=True)
    
    pdir = '/nv/vol190/zanglab/'
    project_dir = '{}/zw5j/since2019_projects/phase_separation_FEpiTR/'.format(pdir)
    motif_dir = '{}/shared/Motif/sites/hg38_fimo_jarspar/results'.format(pdir)
    cluster_dir = '{}/f7_TF_condensates_test/f1_TF_motif_cluster_cor_SE/f2_TF_fimo_jarspar_cluster/f5_motif_within_300bp'.format(project_dir)
    se_dir = '{}/f7_TF_condensates_test/f1_TF_motif_cluster_cor_SE/f4_cluster_enrichment_at_SE/data/SE_hg38'.format(project_dir)
    
    # se_files = glob.glob('{}/*bed'.format(se_dir))
    # se_files.sort()
    # se_basename = os.path.basename(se_file).split('.bed')[0]
    se_file = '{}/{}.bed'.format(se_dir,se_basename)
    cluster_files = glob.glob('{}/*bed'.format(cluster_dir))
    cluster_files.sort()
    

    df = pd.DataFrame()
    for cluster_file in cluster_files[:]:
        motif_basename = os.path.basename(cluster_file)
        motif = motif_basename.split('_')[0]
        motif_file = '{}/{}'.format(motif_dir,motif_basename)
        print(se_basename,motif_basename,motif)
        
        # == whether the motif overlap with SE
        motif_total = get_lines(motif_file)
        motif_overlap = run_bedtools_intersect(motif_file,se_file)
    
        # == whether the cluster motif overlap with SE
        cluster_total = get_lines(cluster_file)
        cluster_overlap = run_bedtools_intersect(cluster_file,se_file)
        
        # == fisher exact test
        a = cluster_overlap
        b = cluster_total-cluster_overlap
        c = motif_overlap - cluster_overlap
        d = motif_total - cluster_total - motif_overlap + cluster_overlap
        print(a,b,c,d)
        s,p = stats.fisher_exact([[a, b],[c, d]])

        df.loc[motif,'motif_total'] = motif_total
        df.loc[motif,'motif_overlap_SE'] = motif_overlap
        df.loc[motif,'cluster_total'] = cluster_total
        df.loc[motif,'cluster_overlap_SE'] = cluster_overlap
        df.loc[motif,'fisher_exact_s'] = s
        df.loc[motif,'fisher_exact_p'] = p
        
    df.index.name='TF'
    df.to_csv('{}/{}.csv'.format(outdir,se_basename))
        
                          
                          

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    #parser.add_argument('-d', '--data', action = 'store', type = str,dest = 'data', help = 'input file of', metavar = '<str>')
    #parser.add_argument('-n', '--normalization', action = 'store', type = str,dest = 'normalization', help = 'input file of', metavar = '<str>')
    #parser.add_argument('-c', '--chrom', action = 'store', type = str,dest = 'chrom', help = 'input file of', metavar = '<str>')
    #parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of', metavar = '<file>')
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input dir of ', metavar = '<file>')
    #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    #parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main(args.infile)



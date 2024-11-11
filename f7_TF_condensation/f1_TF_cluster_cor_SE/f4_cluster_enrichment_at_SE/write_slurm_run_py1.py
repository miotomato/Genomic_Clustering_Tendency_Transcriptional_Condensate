import sys,argparse
import os,glob
import numpy as np
import pandas as pd
from scipy import stats
import re,bisect

chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',\
             'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',\
             'chr18','chr19','chr20','chr21','chr22','chrX','chrY']


def write_slurm(slurm_dir,basename):

    slurmfile = slurm_dir+os.sep+'run_{}.slurm'.format(basename)
    with open(slurmfile,'w') as slurmout:
        slurmout.write('''#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=6000
#SBATCH -t 12:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
''')

        #######################
        # REMEMBER TO CHANGE THE MEMORY !!!
        #########################
            
        slurmout.write('#SBATCH -o {}/slurm_{}.out\n\n'.format(slurm_dir,basename))
        slurmout.write('time python py1_cluster_cor_SE_fisher_exact.py -i {}'.format(basename))



def main():
    
    slurm_dir = 'run_py1_slurms'
    os.makedirs(slurm_dir,exist_ok=True)
    
    pdir = '/nv/vol190/zanglab/'
    project_dir = '{}/zw5j/since2019_projects/phase_separation_FEpiTR/'.format(pdir)
    se_dir = '{}/f7_TF_condensation/f1_TF_cluster_cor_SE/f4_cluster_enrichment_at_SE/data/SE_hg38'.format(project_dir)
    se_files = glob.glob('{}/*bed'.format(se_dir))
    se_files.sort()
    
    for se_file in se_files[:]:
        basename = os.path.basename(se_file).split('.bed')[0]
        write_slurm(slurm_dir,basename)



    
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    #parser.add_argument('-d', '--data', action = 'store', type = str,dest = 'data', help = 'input file of', metavar = '<str>')
    #parser.add_argument('-n', '--normalization', action = 'store', type = str,dest = 'normalization', help = 'input file of', metavar = '<str>')
    #parser.add_argument('-c', '--chrom', action = 'store', type = str,dest = 'chrom', help = 'input file of', metavar = '<str>')
    #parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of', metavar = '<file>')
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    #parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main()

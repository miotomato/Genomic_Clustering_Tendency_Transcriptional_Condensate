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


def write_slurm(slurm_dir,basename,commandline):

    slurmfile = slurm_dir+os.sep+'run_{}.slurm'.format(basename)
    with open(slurmfile,'w') as slurmout:
        slurmout.write('''#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 12:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
''')

        #######################
        # REMEMBER TO CHANGE THE MEMORY !!!
        #########################
            
        slurmout.write('#SBATCH -o {}/slurm_{}.out\n\n'.format(slurm_dir,basename))
#         slurmout.write('module load gcc/7.1.0\n')
#         slurmout.write('module load bedtools/2.26.0\n\n')
        slurmout.write(commandline)

    



# ==== main
slurm_dir = 'slurm_run1b'
os.makedirs(slurm_dir,exist_ok=True)

outdir = 'CP_TFMS_vs_random'
infiles = glob.glob('{}/_csv/*Exponential*'.format(outdir))
infiles = [i for i in infiles if not re.search('~',i)]
infiles.sort()


for infile in infiles[:]:
    basename = os.path.basename(infile).split('_Exponential_Pvalue')[0]
    commanLine = 'time python run1b_TFMS_enrich_at_SE.py -i {}\n\n'.format(infile)
    write_slurm(slurm_dir,basename,commanLine)
    



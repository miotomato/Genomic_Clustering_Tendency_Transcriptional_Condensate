import os,sys,argparse
# import fileinput,time
import glob
import re,bisect
import pandas as pd
import numpy as np
import matplotlib


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
        slurmout.write('module load gcc/7.1.0\n')
        slurmout.write('module load bedtools/2.26.0\n\n')
        slurmout.write(commandline)

    

# ==== main
slurm_dir = 'slurm_run2'
os.makedirs(slurm_dir,exist_ok=True)

# ==== read cistrome TFBS and TFMS info
df_qc = pd.read_csv('../../data/cistrome2019_selected_QC.csv',index_col=0)
df_count = pd.read_csv('../../data/cistrome2019_data_Count_with_SE_motif.csv',index_col=0)
# peak_dir = '../f1_bedtools_closest/data_TFBS'
# motif_dir = '../f1_bedtools_closest/data_TFMS_jarspar/'


# ==== compare TFBS and TFMS per TF per celltype
# df_cp = pd.DataFrame()
for factor in df_count.index[:]:
    for celltype in df_count.columns[:]:
        if df_count.loc[factor,celltype]==0:
            continue
        df_tmp = df_qc[(df_qc.Factor==factor)&(df_qc.Cell_line==celltype)]
        for dcID in df_tmp.index:
            basename = '{}_{}_{}'.format(dcID,factor,celltype)
            commanLine = 'time python run2_TFBS_vs_TFMS_sample.py -i {}\n\n'.format(dcID)
            write_slurm(slurm_dir,basename,commanLine)
            

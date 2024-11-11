import sys,argparse
import os,glob
import numpy as np
import pandas as pd
from scipy import stats
import re,bisect

chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',\
             'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',\
             'chr18','chr19','chr20','chr21','chr22','chrX','chrY']


def write_slurm(slurm_dir,basename,commandline):

    slurmfile = slurm_dir+os.sep+'run_{}.slurm'.format(basename)
    with open(slurmfile,'w') as slurmout:
        slurmout.write('''#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
''')

        #######################
        # REMEMBER TO CHANGE THE MEMORY !!!
        #########################
            
        slurmout.write('#SBATCH -o {}/slurm_{}.out\n\n'.format(slurm_dir,basename))
        slurmout.write(commandline)





    
slurm_dir = 'processed_cell_Regions_slurms'
outdir = 'processed_cell_Regions_per_sample_per_cellType'
os.makedirs(slurm_dir,exist_ok=True)
os.makedirs(outdir,exist_ok=True)

# == collect all count files
sample_countFiles = glob.glob('GSE201336_RAW/*.tsv.gz')
indir = 'processed_cell_ID_per_sample_per_cellType'
catalogTypes = ['immune','stromal']

for catalog in catalogTypes[:]:
    info_file = 'final_cells_included_and_cell_types/{}_celltypes_atac.tsv'.format(catalog)
    info_df = pd.read_csv(info_file,sep='\t')
    samples = list(set(info_df.Sample))
    for sample in sorted(samples)[:]:
        # ==== collect countFile from all replicates
        sample_uniq = sample.split('-R')[0]
        countFiles = [i for i in sample_countFiles if re.search('_{}_'.format(sample_uniq),i)]
        # print(sample,len(countFiles))
        # assert len(countFile)>=1
        df_sample = info_df[info_df.Sample==sample]
        cellTypes = list(set(df_sample.CellType))
        for cellType in cellTypes[:]:
            cellType_rename = '_'.join(re.split(r'\s|\/',cellType))
            os.makedirs('{}/{}/{}'.format(outdir,catalog,cellType_rename),exist_ok=True)
            id_file = '{}/{}/{}/{}.txt'.format(indir,catalog,cellType_rename,sample)
            # print(id_file)
            outfile = '{}/{}/{}/{}.bed'.format(outdir,catalog,cellType_rename,sample)
            basename = '{}_{}_{}'.format(catalog,sample,cellType_rename)
            commandLine = 'zcat {} | grep -f {} > \\\n{}\n'.format(' \\\n'.join(countFiles),id_file,outfile)
            # print(commandLine)
            write_slurm(slurm_dir,basename,commandLine)
            
            










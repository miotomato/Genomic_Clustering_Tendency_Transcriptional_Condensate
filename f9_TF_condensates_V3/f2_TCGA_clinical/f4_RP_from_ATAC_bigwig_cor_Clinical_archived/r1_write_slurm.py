import sys,argparse
import os,glob
import numpy as np
import pandas as pd
from scipy import stats
import re,bisect


    



# ==== main 
    
slurm_dir = 'f1_slurms'
outdir='f1_get_RP_from_bedGraph'
os.makedirs(slurm_dir,exist_ok=True)
os.makedirs(outdir,exist_ok=True)

tcga_bedfile='../../data/TCGA/tcga_atac.bed'
coad_merge_bedfile='../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merge/HCT-116.merge.bed'
brca_merge_bedfile='../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merge/MCF-7.merge.bed'

bdgfiles = glob.glob('../../data/TCGA/bigwig/*.bdg')
for bdgfile in bdgfiles:
    bdg_prename=os.path.basename(bdgfile).split('.bdg')[0]
    print(bdg_prename)
    
    # ==== write the .slurm file
    slurmfile = slurm_dir+os.sep+'run_{}.slurm'.format(bdg_prename)
    with open(slurmfile,'w') as slurmout:
        slurmout.write('''#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 12:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
''')

    #######################
    # REMEMBER TO CHANGE THE MEMORY !!!
    #########################
        
        slurmout.write('#SBATCH -o {}/slurm_{}.out\n\n'.format(slurm_dir,bdg_prename))
        afile = tcga_bedfile
        afile_prename = os.path.basename(afile).split('.bed')[0]
        slurmout.write('time python get_RP_from_bedGraph.py \\\n-a {} \\\n-b {} \\\n-s hg38 -o {}/{}_on_{}.csv\n\n'.format(afile,bdgfile,outdir,bdg_prename,afile_prename))
    
        afile = coad_merge_bedfile
        afile_prename = os.path.basename(afile).split('.bed')[0]
        slurmout.write('time python get_RP_from_bedGraph.py \\\n-a {} \\\n-b {} \\\n-s hg38 -o {}/{}_on_{}.csv\n\n'.format(afile,bdgfile,outdir,bdg_prename,afile_prename))

        afile = brca_merge_bedfile
        afile_prename = os.path.basename(afile).split('.bed')[0]
        slurmout.write('time python get_RP_from_bedGraph.py \\\n-a {} \\\n-b {} \\\n-s hg38 -o {}/{}_on_{}.csv\n\n'.format(afile,bdgfile,outdir,bdg_prename,afile_prename))


    
    

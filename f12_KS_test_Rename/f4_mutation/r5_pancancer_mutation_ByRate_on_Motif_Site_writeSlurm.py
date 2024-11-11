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
#SBATCH --mem=40000
#SBATCH -t 4:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
''')

        #######################
        # REMEMBER TO CHANGE THE MEMORY !!!
        #########################
            
        slurmout.write('#SBATCH -o {}/slurm_{}.out\n\n'.format(slurm_dir,basename))
        slurmout.write('\n'.join(commandline))





    
slurm_dir = 'f5_pancancer_mutation_ByRate_on_Motif_Site_slurm'
outdir = 'f5_pancancer_mutation_ByRate_on_Motif_Site'
os.makedirs(slurm_dir,exist_ok=True)
os.makedirs(outdir,exist_ok=True)

            
            
cts = ['MCF-7','HCT-116','HeLa','LNCaP','U87','HepG2']
datadir = '../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged/'
tfbs_dir = '../f1_TF_cluster_potential/f2_cor_CP_SE_AICAP/f3_TFBS_CP_heatmap//_csv/'
tfbs_cp_s = pd.read_csv('{}/data_fisher_CP_TFBS_all_vs_TFMS_CP_KStest_statistics.csv'.format(tfbs_dir),index_col=0)

mutation_file='../../f9_TF_condensates_V3/f4_mutation/data//pancancer_mutation.bed'
# df = pd.DataFrame()
for ct in cts[:]:
    os.makedirs(outdir+os.sep+ct,exist_ok=True)
    factors = tfbs_cp_s[ct].dropna().index
    for factor in factors[:]:

        motif_dir="/nv/vol190/zanglab/shared/Motif/sites/hg38_fimo_jarspar/results/"
        motif_files = glob.glob('{}/{}_*'.format(motif_dir,factor)) 
        motif_files = [i for i in motif_files if not re.search('~',i)]
        assert len(motif_files)==1
        motif_file = motif_files[0]

        for flag in ['percentile_T','percentile_T_ExtendMerge','percentile_C']:
            mergefile = '{}/{}/{}_{}_{}.merge.bed'.format(datadir,ct,ct,factor,flag)
            if os.path.isfile(mergefile):
                commandLine = []
                motif_overlapped = '{}/{}/{}_{}_{}_Motif.overlapped.bed'.format(outdir,ct,ct,factor,flag)                
                motif_overlapped_merged = '{}/{}/{}_{}_{}_Motif.overlapped.merge.bed'.format(outdir,ct,ct,factor,flag)                
                commandLine.append('time bedtools intersect \\\n-a {} \\\n-b {} \\\n-wa -u > {}\n'.format(motif_file,mergefile,motif_overlapped))
                commandLine.append('time bedtools merge \\\n-i {} \\\n> {}\n'.format(motif_overlapped,motif_overlapped_merged))
                outfile = '{}/{}/{}_{}_{}_mutatioinCount.bed'.format(outdir,ct,ct,factor,flag)
                commandLine.append('time python find_overlap_keep_info_NOT_sep_strand_lastColMarked_keep_bfile_info.py \\\n-a {} \\\n-b {} \\\n-s hg38 -p {}\n'.format(motif_overlapped_merged,mutation_file,outfile))
                # print(ct,factor,flag,mergefile)
                basename = '{}_{}_{}'.format(ct,factor,flag)
                write_slurm(slurm_dir,basename,commandLine)
                
            # mergefile_SE = '{}/{}/{}_{}_{}.merge.SE_overlapped.bed'.format(outdir,ct,ct,factor,flag)










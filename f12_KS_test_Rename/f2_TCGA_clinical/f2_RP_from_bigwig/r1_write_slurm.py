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

mergefile_indir = '../../f1_TF_cluster_potential/f3_clustered_TFBS/f2_bedfiles_merged/'
tcga_bedfile='../../../f9_TF_condensates_V3/data//TCGA/tcga_atac.bed'
name_match = pd.read_excel('../../../f9_TF_condensates_V3/data/TCGA/TCGA-ATAC_SE_cancerType_match.xlsx',index_col=0)   
name_match = name_match.dropna()

# cts = ['MCF-7','HCT-116','HeLa','LNCaP','U87','HepG2']
tfbs_dir = '../..//f1_TF_cluster_potential/f2_cor_CP_SE_AICAP/f3_TFBS_CP_heatmap//_csv/'
tfbs_cp_s = pd.read_csv('{}/data_fisher_CP_TFBS_all_vs_TFMS_CP_KStest_statistics.csv'.format(tfbs_dir),index_col=0)

halflifes = [10000,5000,20000]


for cancertype in name_match.index[:]:
    os.makedirs(outdir+os.sep+cancertype,exist_ok=True)
    ct = name_match.loc[cancertype,'SE']
    factors = tfbs_cp_s[ct].dropna().index

    bdgfiles = glob.glob('../../../f9_TF_condensates_V3/data/TCGA/bigwig/{}/*.bdg'.format(cancertype))
    for bdgfile in bdgfiles[:]:
        bdg_prename=os.path.basename(bdgfile).split('.bdg')[0]
    
        for factor in factors[:]:
            slurmfile = slurm_dir+os.sep+'run_{}_on_{}.slurm'.format(bdg_prename,factor)
            with open(slurmfile,'w') as slurmout:
                slurmout.write('''#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=20000
#SBATCH -t 12:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
''')
                slurmout.write('#SBATCH -o {}/slurm_{}_on_{}.out\n\n'.format(slurm_dir,bdg_prename,factor))
    
                for halflife in halflifes:
                    for flag in ['percentile_T','percentile_T_ExtendMerge','percentile_C']:
                        mergefile = '{}/{}/{}_{}_{}.merge.bed'.format(mergefile_indir,ct,ct,factor,flag)
                        mergefile_SE = '{}/{}/{}_{}_{}.merge.SE_overlapped.bed'.format(mergefile_indir,ct,ct,factor,flag)
                        
                        if os.path.isfile(mergefile):                            
                            afile = mergefile
                            afile_prename = os.path.basename(afile).split('.bed')[0]
                            slurmout.write('time python get_RP_from_bedGraph.py \\\n-a {} \\\n-b {} \\\n-s hg38 -u {} -o {}/{}/{}_halflife_{}_on_{}.csv\n\n'.format(afile,bdgfile,halflife,outdir,cancertype,bdg_prename,halflife,afile_prename))
        
                            afile = mergefile_SE
                            afile_prename = os.path.basename(afile).split('.bed')[0]
                            slurmout.write('time python get_RP_from_bedGraph.py \\\n-a {} \\\n-b {} \\\n-s hg38 -u {} -o {}/{}/{}_halflife_{}_on_{}.csv\n\n'.format(afile,bdgfile,halflife,outdir,cancertype,bdg_prename,halflife,afile_prename))
            

    
    

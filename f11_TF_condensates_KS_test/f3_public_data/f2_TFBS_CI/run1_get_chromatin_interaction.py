import sys,argparse
import os,glob,shutil
import numpy as np
import pandas as pd
import random,string
from scipy import stats
# from bart3d import local_interaction_from_hicpro, local_interaction_from_hic, local_interaction_from_cool
import local_interaction_from_hicpro



def normalized_mirror_concat(df):
    # ==== normalize the contact matrix by dividing averaged reads of bins with the same distance 
    df = df/df.mean()  # median =0 in most cases
    init_cols = df.columns
    # ==== get up-stream chromatin interactions
    for column_pos in np.arange(1,len(init_cols)):
        column = init_cols[column_pos]
        mir_column = '-{}'.format(column)
        new_values = np.append([0]*column_pos,df[column].values)[:len(df.index)]
        df.insert(0,mir_column,new_values)
    return df



def return_matrix_file(workdir,cellType):
    
    dir1 = '{}/since2019_projects/phase_separation_FEpiTR/f11_TF_condensates_KS_test/f3_public_data/data_processed/public_hic_wt/bglii'.format(workdir)
    dir2 = '{}/since2019_projects/phase_separation_FEpiTR/f11_TF_condensates_KS_test/f3_public_data/data_processed/public_hic_wt/dpnii'.format(workdir)
    dir3 = '{}/since2019_projects/phase_separation_FEpiTR/f11_TF_condensates_KS_test/f3_public_data/data_processed/public_hic_wt/hindiii'.format(workdir)
    dir4 = '{}/since2019_projects/phase_separation_FEpiTR/f11_TF_condensates_KS_test/f3_public_data/data_processed/public_hic_wt/mboi'.format(workdir)
    dir5 = '{}/since2019_projects/phase_separation_FEpiTR/f11_TF_condensates_KS_test/f3_public_data/data_processed/GSE104333_hct116_hic/hct116_rad21_auxin/'.format(workdir)
    dir6 = '{}/work2017/T_ALL_CTCF/11_HiC_analysis/f1_preprocess/HiC_Pro/hicPro_raw_fq/panos_hic/'.format(workdir)

    for matrix_dir in [dir1,dir2,dir3,dir4,dir5,dir6]: 
        matrix_file="{}/hic_results/matrix/{}/raw/5000/{}_5000.matrix".format(matrix_dir,cellType,cellType)      
        ord_file="{}/hic_results/matrix/{}/raw/5000/{}_5000_abs.bed".format(matrix_dir,cellType,cellType)   
        if os.path.isfile(matrix_file) and os.path.isfile(ord_file):
            return matrix_file,ord_file
    # return None,None




# ==== main
workdir = '/nv/vol190/zanglab/zw5j/'
# workdir = '/Volumes/zanglab/zw5j/'
outdir = "f1_CI"
os.makedirs(outdir,exist_ok=True)

# choose the right chrom and chrom_len
chromsize = '{}/env/bart3d/bart3d/utility/hg38_clean.chrom.sizes'.format(workdir)
chrom_len = pd.read_csv(chromsize,header=None,index_col=0,sep='\t')[1].to_dict()
chroms = list(chrom_len.keys())
species="hg38"
genomic_KB_distances = [20,50,100,200,500]
cellTypes = ['LNCaP','H1','HepG2','MCF7','Panc1','GM12878','HeLa','IMR90','K562','HCT116_no_auxin_rep1','Jurkat']


    
for genomic_KB_distance in genomic_KB_distances:
    genomic_distance = genomic_KB_distance*1000
    
    for cellType in cellTypes:
        # get hicpro files
        matrix_file,ord_file = return_matrix_file(workdir,cellType)
        # read hicpro results
        interaction_dfs,resolution = local_interaction_from_hicpro.get_local_interaction_from_matrix(ord_file,matrix_file,outdir,chrom_len,genomic_distance)

        # write out the chromatin interactions
        cellType_name = cellType.split('_')[0]
        outfile_name = '{}/{}_CI_{}KB.bed'.format(outdir,cellType_name,genomic_KB_distance)
        outfile = open(outfile_name,'w')
        
        for chr in chroms[:]:
            if chr in interaction_dfs.keys():
                interaction_df = interaction_dfs[chr].astype(float).round(2)
                # ==== get up-stream interactions
                interaction_df = normalized_mirror_concat(interaction_df)
        
                for ii in interaction_df.index:
                    start,end = ii,ii+resolution
                    values = '\t'.join(interaction_df.loc[ii].round(2).astype(str))
                    outfile.write('{}\t{}\t{}\t{}\n'.format(chr,start,end,values))
        
        outfile.close()
            




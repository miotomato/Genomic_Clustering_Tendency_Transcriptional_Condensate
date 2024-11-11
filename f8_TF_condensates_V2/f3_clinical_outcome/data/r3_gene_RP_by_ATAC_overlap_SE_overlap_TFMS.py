import os,sys,argparse,glob,re
import numpy as np
import pandas as pd
# import find_overlap_keep_info_NOT_sep_strand_asimport



# == main 
indir = 'ATAC_overlap_SE_overlap_TFMS_avg_sig'
outdir = 'ATAC_overlap_SE_overlap_TFMS_avg_sig_RP'
os.makedirs(outdir,exist_ok=True)


# project_dir='/Volumes/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/f7_TF_condensates_test'
project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/f7_TF_condensates_test'
motif_file_dir='/nv/vol190/zanglab/shared/Motif/sites/hg38_fimo_jarspar/results/'
# motif_file_dir='/Volumes/zanglab/shared/Motif/sites/hg38_fimo_jarspar/results/'
motif_infiles = glob.glob('{}/*'.format(motif_file_dir))
motif_infiles = [i for i in motif_infiles if not re.search('~',i)]
motif_infiles.sort()


# read matched names between TCGA ATAC and SE data
name_match = pd.read_excel('TCGA-ATAC_SE_cancerType_match.xlsx',index_col=0)   
name_match = name_match.dropna()
for cancertype in ['BRCA','COAD'][:]:
    cancertype_SE = name_match.loc[cancertype].SE
    cancertype_SE_rename = name_match.loc[cancertype].SE_rename

    for motif_file in motif_infiles[:]:
        motif_name = os.path.basename(motif_file).split('.bed')[0]
        bed_file = '{}/{}_ATAC_overlap_{}_SE_overlap_{}.bed'.format(indir,cancertype,cancertype_SE,motif_name)
        out_file = '{}/{}_ATAC_overlap_{}_SE_overlap_{}_RP.tsv'.format(outdir,cancertype,cancertype_SE,motif_name)

        if os.path.isfile(bed_file):
            commandLine = 'python2 get-regulatory-potential-on-genes_peak_level.py \
                -b {} -s hg38 -g hg38_unique_geneSymbol.ucsc -o {}\n'.format(bed_file,out_file,)
            os.system(commandLine);print(commandLine)










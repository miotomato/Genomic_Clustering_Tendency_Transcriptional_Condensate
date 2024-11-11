import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect

# import matplotlib
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt
# matplotlib.rcParams['font.size']=16
# import seaborn as sns
# sns.set(font_scale=1.2)
# sns.set_style("whitegrid", {'axes.grid' : False})
# sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k'})
# matplotlib.rcParams["font.sans-serif"] = ["Arial"]
# matplotlib.rcParams['mathtext.fontset'] = 'custom'
# matplotlib.rcParams["mathtext.rm"] = "Arial"


outdir = 'f1_continuum_changes_at_coBinding'
os.makedirs(outdir,exist_ok=True)

project_dir = '/nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/'
# project_dir = '/Volumes/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/'
public_data_dir = '{}/f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022/'.format(project_dir)
scdata_indir = 'processed_cell_Regions_MERGE_by_GrossPathology'

# ==== binding info in HCT116
tf_binding_dir = '{}/f9_TF_condensates_V3/f1_TF_cluster_potential/f3_clustered_TFBS/f4_cobinding_TFBS_venn/'.format(project_dir)
tf_binding_df = pd.read_csv('{}/HCT-116_cobinding.csv'.format(tf_binding_dir))

co_binding_col = ['HCT-116_SRF', 'HCT-116_JUND', 'HCT-116_CEBPB']
co_binding_df = tf_binding_df[np.sign(tf_binding_df[co_binding_col]).sum(axis=1)==3]
co_binding_bed = '{}/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed'.format(outdir)
co_binding_df[['chr', 'start', 'end']].to_csv(co_binding_bed,index=False,header=None,sep='\t')

se_file = '{}/f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed'.format(project_dir)

# ==== read sample GrossPathology
# sample_info = pd.read_excel('{}/41588_2022_1088_MOESM3_ESM.xlsx'.format(public_data_dir),sheet_name='TableS2',skiprows=[0],index_col=0)
GrossPathologys = ['Adenocarcinoma', 'Normal', 'Polyp', 'Unaffected']
catalogTypes = ['immune','stromal']

for catalog in catalogTypes[:]:
    # outdf = pd.DataFrame()
    info_file = '{}/final_cells_included_and_cell_types/{}_celltypes_atac.tsv'.format(public_data_dir,catalog)
    info_df = pd.read_csv(info_file,sep='\t')
    
    cellTypes = list(set(info_df.CellType))
    for cellType in cellTypes[:]:
        cellType_rename = '_'.join(re.split(r'\s|\/',cellType))
        os.makedirs('{}/{}/{}'.format(outdir,catalog,cellType_rename),exist_ok=True)
        
        for GrossPathology in GrossPathologys[:]:
            bedfile = '{}/{}/{}/{}/{}.bed'.format(public_data_dir,scdata_indir,catalog,cellType_rename,GrossPathology)
            if os.path.isfile(bedfile):
                print(catalog,cellType,GrossPathology)
                
                bfile = co_binding_bed
                outbed_prename = '{}/{}/{}/{}_co_binding'.format(outdir,catalog,cellType_rename,GrossPathology)
                commandLine = 'bedtools intersect -a {} \\\n-b {} \\\n-wa -u > {}_overlapped.bed\n'.format(bedfile,bfile,outbed_prename)
                print(commandLine);os.system(commandLine)
                commandLine = 'bedtools intersect -a {} \\\n-b {} \\\n-wa -v > {}_NOT_overlapped.bed\n'.format(bedfile,bfile,outbed_prename)
                print(commandLine);os.system(commandLine)
        
                bfile = se_file
                outbed_prename = '{}/{}/{}/{}_SE'.format(outdir,catalog,cellType_rename,GrossPathology)
                commandLine = 'bedtools intersect -a {} \\\n-b {} \\\n-wa -u > {}_overlapped.bed\n'.format(bedfile,bfile,outbed_prename)
                print(commandLine);os.system(commandLine)
                commandLine = 'bedtools intersect -a {} \\\n-b {} \\\n-wa -v > {}_NOT_overlapped.bed\n'.format(bedfile,bfile,outbed_prename)
                print(commandLine);os.system(commandLine)
        






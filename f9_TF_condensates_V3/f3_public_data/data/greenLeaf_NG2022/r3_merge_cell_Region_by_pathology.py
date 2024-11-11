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




indir = 'processed_cell_Regions_per_sample_per_cellType'
outdir = 'processed_cell_Regions_MERGE_by_GrossPathology'
# os.makedirs(outdir,exist_ok=True)

# ==== read sample GrossPathology
sample_info = pd.read_excel('41588_2022_1088_MOESM3_ESM.xlsx',sheet_name='TableS2',skiprows=[0],index_col=0)
GrossPathologys = ['Adenocarcinoma', 'Normal', 'Polyp', 'Unaffected']
catalogTypes = ['immune','stromal']

for catalog in catalogTypes[:]:
    # outdf = pd.DataFrame()
    info_file = 'final_cells_included_and_cell_types/{}_celltypes_atac.tsv'.format(catalog)
    info_df = pd.read_csv(info_file,sep='\t')
    # ==== process by each CT
    cellTypes = list(set(info_df.CellType))
    for cellType in cellTypes[:]:
        cellType_rename = '_'.join(re.split(r'\s|\/',cellType))
        os.makedirs('{}/{}/{}'.format(outdir,catalog,cellType_rename),exist_ok=True)
        df_CT = info_df[info_df.CellType==cellType]
        all_samples = list(set(df_CT.Sample))
        # ==== separate all samples by pathology type
        for GrossPathology in GrossPathologys[:]:
            sample_index = sample_info[sample_info.GrossPathology == GrossPathology].index
            samples = [i for i in all_samples if re.split(r'-D|-S',i)[0] in sample_index]
            bedfiles = ['{}/{}/{}/{}.bed'.format(indir,catalog,cellType_rename,sample) for sample in samples]
            bedfile_cat = '{}/{}/{}/{}.bed'.format(outdir,catalog,cellType_rename,GrossPathology)
            if len(bedfiles)>0:
                commandLine = 'cat {} > {}\n'.format(' \\\n'.join(bedfiles),bedfile_cat)
                print(commandLine)
                os.system(commandLine)
            

# for catalog in catalogTypes[:]:
#     for GrossPathology in GrossPathologys[:]:
#         commandLine = 'cat {}/{}/*/{}.bed > {}/{}/All_{}.bed'.format(outdir,catalog,GrossPathology,outdir,catalog,GrossPathology)
#         print(commandLine)
#         os.system(commandLine)





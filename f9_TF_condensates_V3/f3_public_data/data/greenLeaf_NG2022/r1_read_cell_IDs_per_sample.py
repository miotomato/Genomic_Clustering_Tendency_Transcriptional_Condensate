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


catalogTypes = ['immune','stromal','epithelial']
catalogTypes = ['immune','stromal']

outdir = 'processed_cell_ID_per_sample'
for catalog in catalogTypes[:]:
    os.makedirs(outdir+os.sep+catalog,exist_ok=True)
    outdf = pd.DataFrame()
    info_file = 'final_cells_included_and_cell_types/{}_celltypes_atac.tsv'.format(catalog)
    info_df = pd.read_csv(info_file,sep='\t')
    samples = list(set(info_df.Sample))
    for sample in samples[:]:
        df_sample = info_df[info_df.Sample==sample]
        cells = [i.split('#')[1] for i in df_sample.Cell]
        outdf.loc[sample,'#cells'] = len(cells)
        outfile = '{}/{}/{}.txt'.format(outdir,catalog,sample)
        with open(outfile,'w') as outf:
            outf.write('\n'.join(cells)+'\n')
    outdf.to_csv('{}/{}_cellCount.csv'.format(outdir,catalog))



outdir = 'processed_cell_ID_per_sample_per_cellType'
for catalog in catalogTypes[:]:
    outdf = pd.DataFrame()
    info_file = 'final_cells_included_and_cell_types/{}_celltypes_atac.tsv'.format(catalog)
    info_df = pd.read_csv(info_file,sep='\t')
    samples = list(set(info_df.Sample))
    for sample in samples[:]:
        df_sample = info_df[info_df.Sample==sample]
        cellTypes = list(set(df_sample.CellType))
        for cellType in cellTypes:
            cellType_rename = '_'.join(re.split(r'\s|\/',cellType))
            os.makedirs('{}/{}/{}'.format(outdir,catalog,cellType_rename),exist_ok=True)
            print(catalog,sample,cellType,cellType_rename)
            df_sample_CT = df_sample[df_sample.CellType==cellType]
            cells = [i.split('#')[1] for i in df_sample_CT.Cell]
            outdf.loc['{} {}'.format(sample,cellType),'#cells'] = len(cells)
            outfile = '{}/{}/{}/{}.txt'.format(outdir,catalog,cellType_rename,sample)
            with open(outfile,'w') as outf:
                outf.write('\n'.join(cells)+'\n')
    outdf.to_csv('{}/{}_cellCount.csv'.format(outdir,catalog))







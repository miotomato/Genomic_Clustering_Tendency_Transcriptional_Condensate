import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect

import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=16
import seaborn as sns
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k'})
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams["mathtext.rm"] = "Arial"
from scipy import stats
import warnings
warnings.filterwarnings("ignore")

def mark_pvalue(compr_pos,positions,box_vals):
    # s,p = stats.ttest_ind(box_vals[compr_pos[0]],box_vals[compr_pos[1]],nan_policy='omit')
    a,b = box_vals[compr_pos[0]],box_vals[compr_pos[1]]
    s,p = stats.ks_2samp(a,b)
    es = (np.mean(a)-np.mean(b))/(np.std(a+b)) 
    y, h, col = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),95)*.99 ,1.05, 'k'
    y2 = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),3)*0.99
    x1,x2 = positions[compr_pos[0]],positions[compr_pos[1]]
    # p_label='{} {:.1e}'.format(s,p);print(s,p)
    p_label='es:{:.2f}\nks:{}'.format(es,s.round(2))
    # if p_label[-2]=='0':
    #     p_label = p_label[:-2]+p_label[-1]
    # if p>=0.05:
    #     p_label = 'n.s.'
    if compr_pos[2] == 't':
        plt.plot([x1*1.03, x1*1.03, x2*0.97, x2*0.97], [y, y*h, y*h, y], lw=1, c=col)
        plt.text((x1+x2)*.5, y*h, p_label, ha='center', va='bottom', color=col,fontsize=12)
    else:
        plt.plot([x1*1.03, x1*1.03, x2*0.97, x2*0.97], [y2, y2*.91, y2*.91, y2], lw=1, c=col)
        plt.text((x1+x2)*.5, y2*.95, p_label, ha='center', va='top', color=col,fontsize=12)


def box_plot(box_vals,xticklabels,outdir,catalog,cellType):

    positions = np.arange(len(box_vals))
    colors = ['silver']*len(box_vals)
    
    plt.figure(figsize=(3,3))
    g = plt.boxplot(box_vals,positions=positions,widths = .6,patch_artist=True,\
                boxprops=dict(color='k',facecolor='w',fill=None,lw=1),\
                medianprops=dict(color='grey'),showfliers=False)    
                
    # for patch, color in zip(g['boxes'], colors):
    #     patch.set_facecolor(color)

    # scatter_X = []
    # for position_id in np.arange(len(positions)):
    #     scatter_x = np.random.normal(positions[position_id],0.07,len(box_vals[position_id]))
    #     plt.scatter(scatter_x,box_vals[position_id],color=colors[position_id],s=20,zorder=0,alpha=0.99,rasterized=True)

    # positions=[0,1,2]
    for compr_pos in [[0,1,'t'],[2,3,'t'],[4,5,'t'],[6,7,'t']][:int(len(box_vals)/2)]:
        mark_pvalue(compr_pos,positions,box_vals)


    #plt.legend((g[0],g[1]),['a','b'],borderaxespad=0.1,labelspacing=.1)
    #plt.axes().set_xlim(-.5,3)
    # plt.axes().set_xticklabels(xticklabels,fontsize=22)



    # plt.xlim([0,i+1])
    # plt.ylim([-.1,7])
    plt.title(cellType)
    plt.axes().set_xticks(positions)
    plt.axes().set_xticklabels(xticklabels,rotation=45, ha='right',fontsize=16,color='k')
    # plt.legend([g["boxes"][0],g["boxes"][1],g["boxes"][2]],['Control','ENCODE','GTEx'],fontsize=16,borderaxespad=0.1,labelspacing=.1,handletextpad=0.2,loc="upper left",frameon=False)
    # plt.axes().tick_params(axis='y',direction='out', length=4, width=1, colors='black')    
    # plt.ylabel('# junciton region\n intergenic genes',fontsize=18)
    figname = '{}/{}/{}.pdf'.format(outdir,catalog,cellType)
    plt.savefig(figname,bbox_inches='tight',pad_inches=0.1,dpi=600)
    plt.show()
    plt.close()


def add_count_info(bedfile,box_vals,xticklabels,dfout,cellType,flag,GrossPathology):

    df = pd.read_csv(bedfile,header=None,sep='\t')
    # vals = df[4].values
    vals = 1000*df[4]/(df[2]-df[1])
    box_vals.append(vals)
    xticklabels.append('{} {}'.format(GrossPathology,flag))
    
    dfout.loc['{} {}'.format(cellType,flag),'{} Total'.format(GrossPathology)] = len(vals)
    dfout.loc['{} {}'.format(cellType,flag),'{} >2'.format(GrossPathology)] = (vals>10).sum()
    dfout.loc['{} {}'.format(cellType,flag),'{} >2 %'.format(GrossPathology)] = (vals>10).sum()/len(vals)
    
    return box_vals,xticklabels,dfout





# ==== main 

indir = 'f1_continuum_changes_at_coBinding'
outdir = 'f2_continuum_changes_at_coBinding_figs'
os.makedirs(outdir,exist_ok=True)

project_dir = '/nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/'
project_dir = '/Volumes/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/'
public_data_dir = '{}/f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022/'.format(project_dir)
scdata_indir = 'processed_cell_Regions_MERGE_by_GrossPathology'


# ==== read sample GrossPathology
# sample_info = pd.read_excel('{}/41588_2022_1088_MOESM3_ESM.xlsx'.format(public_data_dir),sheet_name='TableS2',skiprows=[0],index_col=0)
GrossPathologys = ['Normal', 'Unaffected', 'Polyp', 'Adenocarcinoma', ]
catalogTypes = ['immune','stromal']

for catalog in catalogTypes[1:]:
    # outdf = pd.DataFrame()
    info_file = '{}/final_cells_included_and_cell_types/{}_celltypes_atac.tsv'.format(public_data_dir,catalog)
    info_df = pd.read_csv(info_file,sep='\t')
    os.makedirs('{}/{}'.format(outdir,catalog),exist_ok=True)
    
    dfout = pd.DataFrame()
    cellTypes = list(set(info_df.CellType))
    # cellTypes = ['CD4+']
    for cellType in cellTypes[:]:
        cellType_rename = '_'.join(re.split(r'\s|\/',cellType))
        
        box_vals = []
        xticklabels =[]
        for GrossPathology in GrossPathologys[:]:
            # bedfile = '{}/{}/{}/{}/{}.bed'.format(public_data_dir,scdata_indir,catalog,cellType_rename,GrossPathology)
            # if os.path.isfile(bedfile):
                # print(catalog,cellType,GrossPathology)
                # bfile = co_binding_bed
            bedfile_overlap = '{}/{}/{}/{}_co_binding_overlapped.bed'.format(indir,catalog,cellType_rename,GrossPathology)
            bedfile_NOT_overlap = '{}/{}/{}/{}_co_binding_NOT_overlapped.bed'.format(indir,catalog,cellType_rename,GrossPathology)
            if os.path.isfile(bedfile_overlap):
                
                box_vals,xticklabels,dfout = add_count_info(bedfile_overlap,box_vals,xticklabels,dfout,cellType,'overlapped',GrossPathology)
                box_vals,xticklabels,dfout = add_count_info(bedfile_NOT_overlap,box_vals,xticklabels,dfout,cellType,'NOT overlapped',GrossPathology)
                
                
                # commandLine = 'bedtools intersect -a {} \\\n-b {} \\\n-wa -u > {}\n'.format(bedfile,bfile,outbed)
                # print(commandLine)
                # os.system(commandLine)
                # dfout.format(cellType,GrossPathology) = (df[4].values>2).sum()
                
        box_plot(box_vals,xticklabels,outdir,catalog,cellType_rename)
        
    dfout.to_csv('{}/{}/data.csv'.format(outdir,catalog))
                
        
                
        
# a,b = box_vals[0], box_vals[1]  
# es = (np.mean(a)-np.mean(b))/(np.std(a+b)) 
# print(es)

# a,b = box_vals[2], box_vals[3] 
# es = (np.mean(a)-np.mean(b))/(np.std(a+b)) 
# print(es)

# a,b = box_vals[4], box_vals[5]  
# es = (np.mean(a)-np.mean(b))/(np.std(a+b)) 
# print(es)

# a,b = box_vals[6], box_vals[7]  
# es = (np.mean(a)-np.mean(b))/(np.std(a+b)) 
# print(es)

    
# s,p = stats.ks_2samp(box_vals[0], box_vals[1])        
# print(s,p)        
# s,p = stats.ranksums(box_vals[2], box_vals[3])        
# print(s,p)        
# s,p = stats.ranksums(box_vals[4], box_vals[5])        
# print(s,p)        
# s,p = stats.ranksums(box_vals[6], box_vals[7])        
# print(s,p)        


a = [1, 2, 1, 2, 1, 2]
b =  [3,33,3,444,2]
a = [1, 2, 1, 2, 1, 2]*88
b =  [3,33,3,444,2]*88
s,p = stats.ks_2samp(a,b)
print(s,p)
# es = (np.mean(a)-np.mean(b))/(np.std(a+b))
# print(es)

# a = [1, 2, 1, 2, 1, 2]*88
# b =  [3,33,3,444,2]*88
# es = (np.mean(a)-np.mean(b))/(np.std(a+b))
# print(es)




import sys,argparse
import os,glob,re
import numpy as np
import pandas as pd
#from GenomeData import *
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=13
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=1)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks")
from matplotlib_venn import venn3,venn2
from scipy import stats
import warnings
warnings.filterwarnings("ignore")


# indir = 'f3_coBinding_merge'
indir = '../../f1_TF_cluster_potential/f3_clustered_TFBS_V2_PeaksGT2k/f4_cobinding_TFBS_venn3'
outdir = 'f3_DCI_overlap_coBinding_TFBS'
os.makedirs(outdir,exist_ok=True)

# atac_file = '../../../f9_TF_condensates_V3/data/TCGA//tcga_atac.bed'
alpha = .7
selected_factors = {'MCF-7 top_TFBSCP':['ERG','ELK1','FOS'],
                    'MCF-7 top_zscored_TFBSCP':['ERG','JUND','ZNF143'],
                    'HCT-116 top_TFBSCP':['ELF1','JUND','SRF',],
                    'HCT-116 top_zscored_TFBSCP':['JUND','CEBPB','YY1']}

genomicDistances = [100000, 200000, 500000]
reps = ['rep1','rep2','all_reps']

for treat_flag in ['percentile_T','percentile_T_ExtendMerge'][:]:
    for ct in ['MCF-7','HCT-116'][:]:
        for factorType in ['top_TFBSCP','top_zscored_TFBSCP'][:]:
            subdir = '{}_{}'.format(ct,factorType)
            os.makedirs(outdir+os.sep+subdir,exist_ok=True)            
            factors = selected_factors['{} {}'.format(ct,factorType)]            
            prenames = factors + ['hg38_exons','hg38_introns','hg38_4k_promoter_geneID']            

            outdf = pd.read_csv('{}/{}/{}_{}_cobinding.csv'.format(indir,subdir,treat_flag,ct,))
       
            # ==== bar plot for genomic annotations
            dfs = [outdf,
                   outdf[(outdf[prenames[0]]!=0)],
                   outdf[(outdf[prenames[1]]!=0)],
                   outdf[(outdf[prenames[2]]!=0)],
                   outdf[(outdf[prenames[0]]!=0)&(outdf[prenames[1]]!=0)],
                   outdf[(outdf[prenames[0]]!=0)&(outdf[prenames[1]]!=0)&(outdf[prenames[2]]!=0)]]
            xticklabels = ['Union',
                           '{}'.format('-'.join([i for i in prenames[:1]])),
                           '{}'.format('-'.join([i for i in prenames[1:2]])),
                           '{}'.format('-'.join([i for i in prenames[2:3]])),
                           '{}'.format('-'.join([i for i in prenames[:2]])),
                           '{}'.format('-'.join([i for i in prenames[:3]])),
                           ]
            
            for ii in np.arange(len(xticklabels)):
                outbed = '{}/{}/{}_{}_{}.bed'.format(outdir,subdir,treat_flag,ct,xticklabels[ii])
                dfs[ii].to_csv(outbed,index=False,sep='\t',header=None)
                
                for genomicDis in genomicDistances: 
                    for rep in reps:
                        dci_file = 'f1_bart3d/RAD21_6hr_auxin_over_NT_{}_dis{}k_differential_score.bed'.format(rep,str(int(genomicDis/1000))) 
                        overlapped_bed = '{}/{}/DCI_RAD21_KO_{}_dis{}k_overlapped_{}_{}_{}.bed'.format(outdir,subdir,rep,str(int(genomicDis/1000)),treat_flag,ct,xticklabels[ii])
                        commandLine = 'bedtools intersect \\\n-a {} \\\n-b {} \\\n-wa -u > {}\n'.format(dci_file,outbed,overlapped_bed)
                        print(commandLine);os.system(commandLine)




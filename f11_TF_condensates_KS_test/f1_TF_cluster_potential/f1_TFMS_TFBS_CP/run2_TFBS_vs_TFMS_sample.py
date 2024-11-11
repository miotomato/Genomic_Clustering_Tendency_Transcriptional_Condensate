import os,sys,argparse
# import fileinput,time
import glob
import re,bisect
import pandas as pd
import numpy as np
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=14
import seaborn as sns
sns.set(font_scale=1.1)
sns.set_style("whitegrid", {'axes.grid' : False})
import scipy
from scipy import stats
import scipy.optimize
sns.set_style("ticks")
matplotlib.rcParams["font.sans-serif"] = ["Arial"]

from module import compr_cluster_potential,enrichment_odds_ratio


   

def main(dcID):
    
    tfbs_types = ['TFBS_all','TFBS_overlap_motif','TFBS_NOT_overlap_motif']
    
    for tfbs_type in tfbs_types:
        outdir = 'CP_{}_vs_TFMS'.format(tfbs_type)
        os.makedirs(outdir+os.sep+'_fig',exist_ok=True)
        os.makedirs(outdir+os.sep+'_csv_sample',exist_ok=True)
        os.makedirs(outdir+os.sep+'_csv_cp',exist_ok=True)
        
        # ==== cistrome info
        df_qc = pd.read_csv('../../data/cistrome/cistrome2019_selected_QC.csv',index_col=0)
        factor = df_qc.loc[dcID].Factor 
        celltype = df_qc.loc[dcID].Cell_line 
        # ==== peak file
        peak_dir = '../f0_bedtools_closest/data_{}'.format(tfbs_type)
        peak_file = '{}/{}_{}_{}.tsv'.format(peak_dir,celltype,factor,dcID)
        # ==== SE data    
        celltype_SE_dir = '../../data/SE_hg38/'
        celltype_SE_file = '{}/{}.bed'.format(celltype_SE_dir,celltype) 
        # ==== motif file
        motif_dir="/nv/vol190/zanglab/shared/Motif/sites/hg38_fimo_jarspar/results/"
        motif_files = glob.glob('{}/{}_*'.format(motif_dir,factor)) 
        motif_files = [i for i in motif_files if not re.search('~',i)]
        assert len(motif_files)==1
        motif_file = motif_files[0]
        # ==== outcome info
        outname = '{}_{}_{}'.format(celltype,factor,dcID)
        labels = [tfbs_type,'TFMS']   
        colors = ['tab:purple','tab:green'] 
        df_out = pd.DataFrame()
        for ii in np.arange(100): 
            outname_ii = '{}_sample{}'.format(outname,ii)   
            df_out,sort_file = compr_cluster_potential(df_out,peak_file,motif_file,outdir,outname_ii,labels,colors)
            df_out = enrichment_odds_ratio(df_out,peak_file,sort_file,celltype_SE_file,outname_ii,labels)
            commandLine = 'rm {} \n'.format(sort_file)
            os.system(commandLine);print(commandLine)  
        df_out.index.name = 'TF'
        # df_cp = df_cp.sort_values(by='log10 distance t-stats',ascending=False)
        df_out.to_csv('{}/_csv_cp/{}.csv'.format(outdir,outname))
        





if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--id', action = 'store', type = int,dest = 'id', help = 'cistrome ID', metavar = '<file>')
    
    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main(args.id)


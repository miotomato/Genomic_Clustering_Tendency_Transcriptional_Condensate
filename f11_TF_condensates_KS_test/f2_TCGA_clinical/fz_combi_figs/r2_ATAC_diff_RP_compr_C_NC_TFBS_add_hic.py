import os,sys,argparse
# import fileinput,time
import glob
import re,bisect
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=12
import seaborn as sns
sns.set(font_scale=1)
sns.set_style("whitegrid", {'axes.grid' : False})
import scipy
import scipy.optimize
sns.set_style("ticks")
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
from matplotlib import gridspec




project_dir = '/nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/f11_TF_condensates_KS_test'
project_dir = '/Volumes/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/f11_TF_condensates_KS_test'
diff_atac_dir = '{}/f2_TCGA_clinical/f1_diff_ATAC/f9b_diff_ATAC_overlap_TFBS_clustered_data_figs'.format(project_dir)
atac_rp_dir = '{}/f2_TCGA_clinical/f4_RP_from_bigwig/f4_avg_RP_per_sample_across_patients_figs'.format(project_dir)
hic_dir = '{}/f3_public_data/f2_TFBS_CI/f3_CI_figs'.format(project_dir)


tfbs_file = '{}/f1_TF_cluster_potential/f3_clustered_TFBS_V2_PeaksGT2k/f2_bedfiles_merged/data_merged_SE_overlapped.csv'.format(project_dir)
tfbs_df = pd.read_csv(tfbs_file,index_col = 0)

name_match = pd.read_excel('{}/../f9_TF_condensates_V3/data/TCGA/TCGA-ATAC_SE_cancerType_match.xlsx'.format(project_dir),index_col=0)
name_match = name_match.dropna()

outdir = 'f2_ATAC_HiC_compr_TFBS_heatmap'
os.makedirs(outdir,exist_ok=True)

removed_tfs = ['HSF1','T','NR2C2']
cancertypes = ['BRCA', 'CESC', 'COAD', 'LIHC', 'PRAD']
treat_flags = ['percentile_T','percentile_T_ExtendMerge']


rank_dir = '../../f1_TF_cluster_potential/f2_cor_CP_SE_AICAP_V2_PeaksGT2k/f9_per_CT_TFBS_CP_cor_zscore_CP/TFBS_CP/'
# tfbs_cp_s = pd.read_csv('{}/data_fisher_CP_TFBS_all_vs_TFMS_CP_KStest_statistics.csv'.format(tfbs_dir),index_col=0)

genomic_dis_kbs = [20,50,100,200,500]


for genomic_dis in genomic_dis_kbs:
    for treat_flag in treat_flags[:1]:
        for cancertype in cancertypes[:1]:
            ct = name_match.loc[cancertype,'SE']
            atac_diff = pd.read_csv('{}/{}_by_{}.csv'.format(diff_atac_dir,cancertype,treat_flag),index_col=0)
            atac_rp = pd.read_csv('{}/{}_halflife_10000_by_{}.merge.csv'.format(atac_rp_dir,cancertype,treat_flag),index_col=0)
            hic_df = pd.read_csv('{}/_{}_CI_{}KB.csv'.format(hic_dir,ct,genomic_dis),index_col=0)
            
            rank_df = pd.read_csv('{}/_CP_TFBS_all_vs_TFMS_{}.csv'.format(rank_dir,ct),index_col=0)
            # enrichment at SE
            out_df = pd.DataFrame()
            # kept_factors = [i for i in atac_rp.index if i not in removed_tfs]
            kept_factors = [i for i in rank_df.index if '{} {}'.format(ct,i) in tfbs_df.index]
            for kept_factor in kept_factors:
                tfbs_index = '{} {}'.format(ct,kept_factor)
                # tfbs_df_ct = tfbs_df.loc[tfbs_index]
                a = tfbs_df.loc[tfbs_index,'# {} on SE'.format(treat_flag)]
                b = tfbs_df.loc[tfbs_index,'# {}'.format(treat_flag)]
                c = tfbs_df.loc[tfbs_index,'# percentile_C on SE'] 
                d = tfbs_df.loc[tfbs_index,'# percentile_C']
                s,p = stats.fisher_exact([[a,b],[c,d]])
                out_df.loc[kept_factor,'SE enrich statistics'] = s
                out_df.loc[kept_factor,'SE enrich pvalue'] = p
            
            out_df['SE enrich pvalue'] = out_df['SE enrich pvalue'].replace(0,1e-299)
            out_df['differential atac statistics'] = atac_diff['statistics']
            out_df['differential atac pvalue'] = atac_diff['pvalue'].replace(0,1e-299)
            out_df['atac RP statistics'] = atac_rp['statistics']
            out_df['atac RP pvalue'] = atac_rp['pvalue'].replace(0,1e-299)
            out_df['hic CI statistics'] = hic_df['statistics']
            out_df['hic CI pvalue'] = hic_df['pvalue'].replace(0,1e-299)
            out_df.to_csv(outdir+os.sep+'data_{}_{}_{}KB.csv'.format(cancertype,treat_flag,genomic_dis))
    
    
            # == plot the heatmap
            plt.figure(figsize=(out_df.shape[0]/3.5,3))
            width_ratios = [1,.3,.5]
            height_ratios = [.35,.7,1,1,1,1]
            gs = gridspec.GridSpec(6,3,width_ratios = width_ratios,
                                   height_ratios = height_ratios, wspace=.15,hspace=.0)
            
            vmax = 5
            vmin = -5
            pal = plt.cm.bwr
            norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
            color_map = matplotlib.cm.ScalarMappable(norm=norm, cmap=pal)
            scatter_scale = 10
               
            ax = plt.subplot(gs[2,:])   
            for ii in np.arange(len(kept_factors)):
                statistics = out_df.loc[kept_factors[ii],'SE enrich statistics']
                pvalue = out_df.loc[kept_factors[ii],'SE enrich pvalue']
                color = color_map.to_rgba(statistics)
                size = np.sqrt(-1*np.log10(pvalue))#;print(size)
                ax.scatter(ii,0,s = scatter_scale*size, color = color)        
            ax.set_yticks([])
            ax.set_xticks([])
            ax.set_xlim([-.5,len(kept_factors)-.5])
            if cancertype=='BRCA':
                ax.set_ylabel('SE Enrichment \n C-TFBS/NC-TFBS',rotation=0,fontsize=13,ha='right',va='center')
            ax.set_title(cancertype,fontsize=14,loc='left')
            
    
            ax = plt.subplot(gs[3,:])   
            for ii in np.arange(len(kept_factors)):
                statistics = out_df.loc[kept_factors[ii],'atac RP statistics']
                pvalue = out_df.loc[kept_factors[ii],'atac RP pvalue']
                color = color_map.to_rgba(statistics)
                size = np.sqrt(-1*np.log10(pvalue))#;print(size)
                ax.scatter(ii,0,s = scatter_scale*size, color = color)        
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_xlim([-.5,len(kept_factors)-.5])
            if cancertype=='BRCA':
                ax.set_ylabel('ATAC PR \n C-TFBS/NC-TFBS',rotation=0,fontsize=13,ha='right',va='center')
    
    
            ax = plt.subplot(gs[4,:])   
            for ii in np.arange(len(kept_factors)):
                statistics = out_df.loc[kept_factors[ii],'differential atac statistics']
                pvalue = out_df.loc[kept_factors[ii],'differential atac pvalue']
                color = color_map.to_rgba(statistics)
                size = np.sqrt(-1*np.log10(pvalue))#;print(size)
                ax.scatter(ii,0,s = scatter_scale*size, color = color)        
            ax.set_yticks([])
            ax.set_xticks([])
            # ax.set_xticks(np.arange(len(kept_factors)))
            # ax.set_xticklabels(kept_factors,rotation=60,fontsize=13,ha='right')
            ax.set_xlim([-.5,len(kept_factors)-.5])
            if cancertype=='BRCA':
                ax.set_ylabel('$\Delta$ ATAC \n C-TFBS/NC-TFBS',rotation=0,fontsize=13,ha='right',va='center')
    
    
            ax = plt.subplot(gs[5,:])   
            for ii in np.arange(len(kept_factors)):
                statistics = out_df.loc[kept_factors[ii],'hic CI statistics']
                pvalue = out_df.loc[kept_factors[ii],'hic CI pvalue']
                color = color_map.to_rgba(statistics)
                size = np.sqrt(-1*np.log10(pvalue))#;print(size)
                ax.scatter(ii,0,s = scatter_scale*size, color = color)        
            ax.set_yticks([])
            ax.set_xticks(np.arange(len(kept_factors)))
            ax.set_xticklabels(kept_factors,rotation=60,fontsize=13,ha='right')
            ax.set_xlim([-.5,len(kept_factors)-.5])
            if cancertype=='BRCA':
                ax.set_ylabel('Hi-C CI \n C-TFBS/NC-TFBS',rotation=0,fontsize=13,ha='right',va='center')
    
                  
            if cancertype=='BRCA':
                # color bar
                ax = plt.subplot(gs[0,1]) 
                cbar = matplotlib.colorbar.ColorbarBase(ax,cmap=pal,norm=norm,orientation='horizontal')
                cbar.set_ticks([vmin,0,vmax])
                cbar.set_label('Statistis ',labelpad=-40, fontsize=13,va='center')
                # # pvalue legend
                ax = plt.subplot(gs[0,2])   
                ax.scatter(.2,0,s = scatter_scale*np.sqrt(-1*np.log10(1e-10)), color = 'k')
                ax.scatter(.5,0,s = scatter_scale*np.sqrt(-1*np.log10(1e-50)), color = 'k')
                ax.scatter(.8,0,s = scatter_scale*np.sqrt(-1*np.log10(1e-100)), color = 'k')
                ax.text(.2,-.2,'1e-10',ha='center',fontsize=11)
                ax.text(.5,-.2,'1e-50',ha='center',fontsize=11)
                ax.text(.8,-.2,'1e-100',ha='center',fontsize=11)
                # ax.set_ylim([-.5,.5])
                ax.set_xlim([-.0,1.])
                ax.set_title('P-value',ha='center',fontsize=13)
                ax.axis('off')
            
            plt.savefig('{}/fig_{}_{}_{}KB.pdf'.format(outdir,cancertype,treat_flag,genomic_dis),bbox_inches='tight',pad_inches=0.02,transparent=True)
            plt.show()
            plt.close()
            
            


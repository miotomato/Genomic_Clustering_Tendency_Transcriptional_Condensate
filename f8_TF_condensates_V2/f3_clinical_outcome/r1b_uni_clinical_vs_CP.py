import os,sys,argparse,glob,re
import numpy as np
import pandas as pd
# import find_overlap_keep_info_NOT_sep_strand_asimport
from collections import Counter
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=14
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams["mathtext.rm"] = "Arial"
import seaborn as sns
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks")
from scipy import stats



# ==== main 
# indir = 'data/ATAC_overlap_SE_overlap_TFMS'
outdir = 'f1_ATAC_overlap_SE_TFMS_clinical_uni'
os.makedirs(outdir+os.sep+'_csv',exist_ok=True)

# ==== read matched names between TCGA ATAC and SE data
name_match = pd.read_excel('data/TCGA-ATAC_SE_cancerType_match.xlsx',index_col=0)   
name_match = name_match.dropna()
# ==== read TFBS CP info
tfbs_dir = '../f1_TF_cluster_potential/f3_CP_cor_SE/f1_CP_heatmap/_csv/'
tfbs_cp_type = 'CP_TFBS_vs_TFMS'
cp_type = 'TFBS_CP'


for data_type in ['SE','CP'][:]:
    for cancertype in ['BRCA','COAD'][:]:
        cancertype_SE = name_match.loc[cancertype].SE
        cancertype_SE_rename = name_match.loc[cancertype].SE_rename
        clinical_file = '{}/{}_clinical_summary.csv'.format(outdir,cancertype)
        clinical_df = pd.read_csv(clinical_file,index_col=0)
        clinical_df.index = [i.split('_') [0] for i in clinical_df.index]
        # atac_overlap_SE_file='{}/{}_ATAC_overlap_{}_SE_caseID.bed'.format(atac_overlap_SE_dir,cancertype,cancertype_SE)
               
        # == then select TFs with high TFBS CP
        tfbs_cp_s = pd.read_csv('{}/CP_data_median_fisher_{}_OR.csv'.format(tfbs_dir,tfbs_cp_type),index_col=0)
        tfbs_se_s = pd.read_csv('{}/SE_data_median_fisher_{}_OR.csv'.format(tfbs_dir,tfbs_cp_type),index_col=0)
        y = tfbs_cp_s[[cancertype_SE]].dropna()#;print(cancertype,len(y))
        ylabel = 'TFBS CP'
        if data_type == 'SE':
            y = tfbs_se_s[[cancertype_SE]].dropna()#;print(cancertype,len(y))
            ylabel = 'TFBS enrichment at SE'
        
        # == save the data
        df = pd.concat([clinical_df,y],axis=1,join='inner')
        df = df.sort_values(by = [cancertype_SE],ascending=False)
        df.to_csv('{}/{}_{}_vs_clinical.csv'.format(outdir,cancertype,data_type))
        
        y = df[cancertype_SE]
        x = 100*df['#P<0.05']/df['total']
        xlabel = '% logrank P<0.05'
        
        vmax = df.shape[0]
        norm = matplotlib.colors.Normalize(vmin=0,vmax=vmax)
        pal = plt.cm.rainbow_r
        color_map = matplotlib.cm.ScalarMappable(norm=norm, cmap=pal)
        colors = [color_map.to_rgba(i) for i in np.arange(vmax)]

        # ==== scatter plot with density
        plt.figure(figsize=(2.6,2.6))
        plt.scatter(x,y,c='k',s=9)
        # label the TFs
        label_i=0
        for label_index in df.index:
            plt.scatter(x[label_index],y[label_index],c=colors[label_i],
                        s=25,label=label_index)
            label_i+=1
        plt.legend(bbox_to_anchor=[.99,1.03],ncol=3,
                markerscale=1.2,fontsize=9,borderaxespad=0.2,labelspacing=.2,
                handletextpad=0.2,handlelength=1.,loc="upper left",frameon=False)
        # ==== linear regression
        slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)       
        x_sort = np.sort(x)
        plt.plot(x_sort,x_sort*slope+intercept,c = 'grey',ls='--',lw=.6)
        plt.text(.02,.88,'$R^2$ = {:.2f}'.format(r_value**2),fontsize=12,transform=plt.axes().transAxes)
        # plt.axhline(y=0,color='k',lw=1.2,ls='--')
        # plt.axvline(x=0,color='k',lw=1.2,ls='--')
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title('{}'.format(cancertype))
        plt.savefig('{}/{}_{}_vs_clinical.pdf'.format(outdir,cancertype,data_type),bbox_inches='tight',pad_inches=0.02,transparent=True)
        plt.show()
        plt.close()
    



            

  


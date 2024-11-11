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


def mark_pvalue(compr_pos,positions,box_vals):
    s,p = stats.ttest_ind(box_vals[compr_pos[0]],box_vals[compr_pos[1]],nan_policy='omit')
    y, h, col = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),96)*.99 ,1.05, 'k'
    y2 = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),3)*0.99
    x1,x2 = positions[compr_pos[0]],positions[compr_pos[1]]
    p_label='{:.1e}'.format(p)
    if p_label[-2]=='0':
        p_label = p_label[:-2]+p_label[-1]
    if p>=0.05:
        p_label = 'n.s.'
    if compr_pos[2] == 't' and s<0 and p<0.05:
        # plt.plot([x1*1.03, x1*1.03, x2*0.97, x2*0.97], [y, y*h, y*h, y], lw=1, c=col)
        plt.text(x2, y*h, '*', ha='center', va='bottom', color='r',fontsize=22)
    # else:
    #     plt.plot([x1*1.03, x1*1.03, x2*0.97, x2*0.97], [y2, y2*.91, y2*.91, y2], lw=1, c=col)
    #     plt.text((x1+x2)*.5, y2*.95, p_label, ha='center', va='top', color=col,fontsize=12)


def box_plot(box_vals,xticklabels,figname,cancertype,data_type):

    positions = np.arange(len(box_vals))
    colors = ['grey']*len(box_vals)
    
    plt.figure(figsize=(4,3))
    g = plt.boxplot(box_vals,positions=positions,widths = .6,patch_artist=True,\
                boxprops=dict(color='k',facecolor='w',fill=None,lw=1),\
                medianprops=dict(color='grey'),showfliers=False)    
                
    for patch, color in zip(g['boxes'], colors):
        patch.set_facecolor(color)

    # scatter_X = []
    # for position_id in np.arange(len(positions)):
    #     scatter_x = np.random.normal(positions[position_id],0.07,len(box_vals[position_id]))
    #     plt.scatter(scatter_x,box_vals[position_id],color=colors[position_id],s=20,zorder=0,alpha=0.99,rasterized=True)

    for pos_test in positions[1:]:
        for compr_pos in [[0,pos_test,'t']]:
            mark_pvalue(compr_pos,positions,box_vals)

    plt.title('{} {}'.format(cancertype,data_type))
    plt.axes().set_xticklabels(xticklabels,rotation=45, ha='right',fontsize=14,color='k')
    # plt.legend([g["boxes"][0],g["boxes"][1],g["boxes"][2]],['Control','ENCODE','GTEx'],fontsize=16,borderaxespad=0.1,labelspacing=.1,handletextpad=0.2,loc="upper left",frameon=False)
    # plt.axes().tick_params(axis='y',direction='out', length=4, width=1, colors='black')    
    plt.ylabel('-log10(log-rank P) ',fontsize=14)
    plt.savefig(figname,bbox_inches='tight',pad_inches=0.1,dpi=600)
    plt.show()
    plt.close()


def add_logrank_info(out_df,df,motif_name):

    df_p = df[df['log rank p']<0.05]
    df_pt = df[(df['treat time']<df['ctrl time'])&(df['log rank p']<0.05)]
    out_df.loc[motif_name,'total'] = df.shape[0]
    out_df.loc[motif_name,'#P<0.05'] = df_p.shape[0]
    out_df.loc[motif_name,'%P<0.05'] = np.round(df_p.shape[0]/df.shape[0],2)
    out_df.loc[motif_name,'#TreatTime<CtrlTime&P<0.05'] = df_pt.shape[0]
    out_df.loc[motif_name,'%TreatTime<CtrlTime&P<0.05'] = np.round(df_pt.shape[0]/df.shape[0],2)
    return out_df




# ==== main 
indir = 'data/ATAC_overlap_SE_overlap_TFMS'
outdir = 'f2_ATAC_overlap_SE_TFMS_clinical_multi'

os.makedirs(outdir+os.sep+'_csv',exist_ok=True)

project_dir='/Volumes/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/f7_TF_condensates_test'
# project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/f7_TF_condensates_test'
atac_overlap_SE_dir = '{}/f6_revised_TCGA_ATAC_cor_SE/f3_clinical_hicor_SE_overlapped_peaks/f1_clinical_at_each_peak/f1_ATAC_overlap_SE_caseID'.format(project_dir)
clinical_dir = '{}/f6_revised_TCGA_ATAC_cor_SE/f3_clinical_hicor_SE_overlapped_peaks/f1_clinical_at_each_peak/f2_caseID_each_SE_vs_clinical'.format(project_dir)
atac_file = '{}/f6_revised_TCGA_ATAC_cor_SE/data/TCGA/tcga_atac.bed'.format(project_dir)
atac_df = pd.read_csv(atac_file,sep='\t',header=None,index_col=3)


# ==== read matched names between TCGA ATAC and SE data
name_match = pd.read_excel('data/TCGA-ATAC_SE_cancerType_match.xlsx',index_col=0)   
name_match = name_match.dropna()
# ==== read TFBS CP info
tfbs_dir = '../f1_TF_cluster_potential/f3_CP_cor_SE/f1_CP_heatmap/_csv/'
tfbs_cp_type = 'CP_TFBS_vs_TFMS'
cp_type = 'TFBS_CP'
tf_thre = 5

for data_type in ['SE','CP'][1:]:
    for cancertype in ['BRCA','COAD'][1:]:
        cancertype_SE = name_match.loc[cancertype].SE
        cancertype_SE_rename = name_match.loc[cancertype].SE_rename
        # atac_overlap_SE_file='{}/{}_ATAC_overlap_{}_SE_caseID.bed'.format(atac_overlap_SE_dir,cancertype,cancertype_SE)
        box_vals,xticklabels = [],[]
        clinical_file = '{}/{}_logrank_info.csv'.format(clinical_dir,cancertype)
        clinical_df = pd.read_csv(clinical_file,index_col=0)
        # atac_overlap_SE_file='{}/{}_ATAC_overlap_{}_SE_caseID.bed'.format(atac_overlap_SE_dir,cancertype,cancertype_SE)
        
        # == save box values
        # pvalues = clinical_df[clinical_df['treat time']<clinical_df['ctrl time']]['log rank p'].values
        
        pvalues = clinical_df['log rank p'].values
        box_vals.append(-1*np.log10(pvalues))
        xticklabels.append('All SEs')
        
        # == then select TFs with high TFBS CP
        tfbs_cp_s = pd.read_csv('{}/CP_data_median_fisher_{}_OR.csv'.format(tfbs_dir,tfbs_cp_type),index_col=0)
        tfbs_se_s = pd.read_csv('{}/SE_data_median_fisher_{}_OR.csv'.format(tfbs_dir,tfbs_cp_type),index_col=0)
        y = tfbs_cp_s[cancertype_SE].dropna()#;print(cancertype,len(y))
        selected_tfs = y.sort_values(ascending=False).index[:tf_thre]
        if data_type == 'SE':
            x = tfbs_se_s[cancertype_SE].dropna()#;print(cancertype,len(y))
            selected_tfs = x.sort_values(ascending=False).index[:tf_thre]
       
        # == save the log rank P info of each TF
        out_df = pd.DataFrame()
        out_df = add_logrank_info(out_df,clinical_df,'All SEs')
        region_count = []
        for motif_name in selected_tfs[:tf_thre]:
            overlapped_beds = glob.glob('{}/{}_ATAC_overlap_{}_SE_overlap_{}_*.bed'.format(indir,cancertype,cancertype_SE,motif_name))
            overlapped_beds = [i for i in overlapped_beds if not re.search('~',i)]
            assert len(overlapped_beds) ==1
            
            try:
                overlapped_df = pd.read_csv(overlapped_beds[0],sep='\t',index_col=3,header=None)
                region_count = np.append(region_count,overlapped_df.index)
                
                # == save box values
                df = clinical_df.loc[overlapped_df.index]
                # pvalues = df[df['treat time']<df['ctrl time']]['log rank p'].values
                pvalues = df['log rank p'].values
                box_vals.append(-1*np.log10(pvalues))
                xticklabels.append(motif_name)
                df = pd.concat([atac_df,df],axis=1,join='inner')
                df.to_csv('{}/_csv/top{}_{}_overlap_SE_{}_ranksum.csv'.format(outdir,data_type,cancertype,motif_name))
                out_df = add_logrank_info(out_df,df,motif_name)
            except:
                print(motif_name)
                
        for count_thre in np.arange(2,len(selected_tfs)):
            motif_name = 'CoBinding_{}'.format(int(count_thre))
            count = dict(Counter(region_count))
            count_key = [i for i in count.keys() if count[i]>=count_thre]
            
            # == save box values
            df = clinical_df.loc[count_key]
            # pvalues = df[df['treat time']<df['ctrl time']]['log rank p'].values
            pvalues = df['log rank p'].values
            box_vals.append(-1*np.log10(pvalues))
            xticklabels.append('CoBinding$≥{}$'.format(int(count_thre)))
            
            df = pd.concat([atac_df,df],axis=1,join='inner')
            df.to_csv('{}/_csv/top{}_{}_overlap_SE_{}_ranksum.csv'.format(outdir,data_type,cancertype,motif_name))
            out_df.loc[motif_name,'total'] = df.shape[0]
            try:
                out_df = add_logrank_info(out_df,df,motif_name)
            except:
                pass
        
        
        out_df.to_csv('{}/top{}_{}_clinical_summary.csv'.format(outdir,data_type,cancertype))
        figname = '{}/top{}_{}_clinical_summary.pdf'.format(outdir,data_type,cancertype)
        box_plot(box_vals,xticklabels,figname,cancertype,data_type)    
        
    



            

  


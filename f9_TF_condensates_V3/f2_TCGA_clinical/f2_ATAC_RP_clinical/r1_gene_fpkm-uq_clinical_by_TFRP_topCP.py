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
sns.set(font_scale=1.1)
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
    
    plt.figure(figsize=(3,2.5))
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
    plt.axes().set_xticklabels(xticklabels,rotation=45, ha='right',fontsize=13,color='k')
    # plt.legend([g["boxes"][0],g["boxes"][1],g["boxes"][2]],['Control','ENCODE','GTEx'],fontsize=16,borderaxespad=0.1,labelspacing=.1,handletextpad=0.2,loc="upper left",frameon=False)
    # plt.axes().tick_params(axis='y',direction='out', length=4, width=1, colors='black')    
    plt.ylabel('-log10(log-rank P) ',fontsize=13)
    plt.savefig(figname,bbox_inches='tight',pad_inches=0.1,dpi=600)
    plt.show()
    plt.close()



def stack_bar(box_vals,xticklabels,figname,cancertype,data_type):
    
    # print(annotation)#;exit()
    plt.figure(figsize=(7,3))
    # cancertypes = sorted(annotation.keys())
    # cancertypes = annotation.keys()
    colors = ['lightgreen','cornflowerblue','lightgrey','lightsalmon'][::-1]
    labels = ["Promoter","Exon","Intron","Distal"][::-1]
    positions = np.arange(len(box_vals))
    
    t_all = len(box_vals[0])
    p_all = (box_vals[0]<p_thre).sum() 
    
    for position in positions:
        vals = box_vals[position]
        total = len(vals)
        a = (vals<p_thre).sum() 
        
        s,p = stats.fisher_exact([[a,total-a],[p_all-a,t_all-p_all-total+a]])
        p_label='{:.1e}'.format(p)
        if p_label[-2]=='0':
            p_label = p_label[:-2]+p_label[-1]
        if p<0.05:
            # plt.text(position,100*a/total, p_label , ha='center', va='bottom', color='k',fontsize=10)
            plt.text(position,92*a/total, '*' , ha='center', va='bottom', color='tab:red',fontsize=20)
        
        g0 = plt.bar(position,100*a/total,bottom=0,width = .68, 
                     lw=0,color = 'tab:purple',alpha=.6,label = labels[0])
        

    plt.axes().set_xticks(positions) 
    plt.title('{} {}'.format(cancertype,data_type))
    plt.axes().set_xticklabels(xticklabels,rotation=45, ha='right',fontsize=13,color='k')
    plt.ylabel('% genes w/ logrank P<0.05 ',fontsize=13)
    plt.savefig(figname,bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
    plt.show()
    plt.close()




def add_logrank_info(out_df,df,motif_name,cp_fig=False):

    df_p = df[df['log rank p']<p_thre]
    if cp_fig:
        for gene in df_p.index[:10]:
            fig_file = '{}/fig/{}_{}_survival.pdf'.format(clinical_dir,cancertype,gene)
            os.system('cp {} {}/_fig/'.format(fig_file,outdir))
    df_pt = df[(df['treat time']<df['ctrl time'])&(df['log rank p']<p_thre)]
    out_df.loc[motif_name,'total'] = df.shape[0]
    out_df.loc[motif_name,'#P<{}'.format(p_thre)] = df_p.shape[0]
    out_df.loc[motif_name,'%P<{}'.format(p_thre)] = np.round(df_p.shape[0]/df.shape[0],2)
    out_df.loc[motif_name,'#TreatTime<CtrlTime&P<{}'.format(p_thre)] = df_pt.shape[0]
    out_df.loc[motif_name,'%TreatTime<CtrlTime&P<{}'.format(p_thre)] = np.round(df_pt.shape[0]/df.shape[0],2)
    return out_df




# ==== main 

project_dir='/Volumes/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR/'
# project_dir='/nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR'
# atac_overlap_SE_dir = '{}/f7_TF_condensates_test/f6_revised_TCGA_ATAC_cor_SE/f3_clinical_hicor_SE_overlapped_peaks/f1_clinical_at_each_peak/f1_ATAC_overlap_SE_caseID'.format(project_dir)
clinical_dir = '{}/f7_TF_condensates_test/f6_revised_TCGA_ATAC_cor_SE/f5_clinical_per_gene/f1_clinical_per_gene/'.format(project_dir)
# atac_file = '{}/f7_TF_condensates_test/f6_revised_TCGA_ATAC_cor_SE/data/TCGA/tcga_atac.bed'.format(project_dir)
# atac_df = pd.read_csv(atac_file,sep='\t',header=None,index_col=3)
rp_dir = '{}/f8_TF_condensates_V2/f3_clinical_outcome/data/ATAC_overlap_SE_overlap_TFMS_avg_sig_RP'.format(project_dir)


outdir = 'f1_gene_fpkm-uq_clinical_by_TFRP_topCP'
os.makedirs(outdir+os.sep+'_csv',exist_ok=True)
os.makedirs(outdir+os.sep+'_fig',exist_ok=True)


# ==== read matched names between TCGA ATAC and SE data
name_match = pd.read_excel('../../data/TCGA/TCGA-ATAC_SE_cancerType_match.xlsx',index_col=0)   
name_match = name_match.dropna()
# ==== read TFBS CP info
tfbs_dir = '../../f1_TF_cluster_potential/f2_cor_CP_SE_AICAP/f3_TFBS_CP_heatmap/_csv/'
tfbs_cp_type = 'CP_TFBS_vs_TFMS'
# cp_type = 'TFBS_CP'
tf_thre = None
rp_thre = 1
p_thre = 0.05
# motif_file_dir='/nv/vol190/zanglab/shared/Motif/sites/hg38_fimo_jarspar/results/'
# motif_file_dir='/Volumes/zanglab/shared/Motif/sites/hg38_fimo_jarspar/results/'
# motif_infiles = glob.glob('{}/*'.format(motif_file_dir))
# motif_infiles = [i for i in motif_infiles if not re.search('~',i)]
# motif_infiles.sort()


for data_type in ['SE','CP'][1:]:
    for cancertype in ['BRCA','COAD'][:]:
        cancertype_SE = name_match.loc[cancertype].SE
        cancertype_SE_rename = name_match.loc[cancertype].SE_rename
        # atac_overlap_SE_file='{}/{}_ATAC_overlap_{}_SE_caseID.bed'.format(atac_overlap_SE_dir,cancertype,cancertype_SE)
        
        # == logrank pvalue of all genes
        clinical_file = '{}/{}_logrank_info.csv'.format(clinical_dir,cancertype)
        clinical_df = pd.read_csv(clinical_file,index_col=0)
        # clinical_df = clinical_df[(clinical_df['treat time']<clinical_df['ctrl time'])]
        
        # == save the log rank P info of each gene
        out_df = pd.DataFrame()
        out_df = add_logrank_info(out_df,clinical_df,'All Genes')
        
        box_vals,xticklabels = [],[]
        pvalues = clinical_df['log rank p'].values
        box_vals.append(pvalues)
        # box_vals.append(-1*np.log10(pvalues))
        xticklabels.append('All Genes')
        
        # == then select TFs with high TFBS CP
        tfbs_cp_s = pd.read_csv('{}/data_fisher_{}_CP_RankSum_statistics.csv'.format(tfbs_dir,tfbs_cp_type),index_col=0)
        tfbs_se_s = pd.read_csv('{}/data_fisher_{}_SE_statistics.csv'.format(tfbs_dir,tfbs_cp_type),index_col=0)
        # tfbs_cp_s = pd.read_csv('{}/CP_data_median_fisher_{}_OR.csv'.format(tfbs_dir,tfbs_cp_type),index_col=0)
        # tfbs_se_s = pd.read_csv('{}/SE_data_median_fisher_{}_OR.csv'.format(tfbs_dir,tfbs_cp_type),index_col=0)
        y = tfbs_cp_s[cancertype_SE].dropna()#;print(cancertype,len(y))
        selected_tfs = y.sort_values(ascending=False).index[:tf_thre] ;print(selected_tfs)
        if data_type == 'SE':
            x = tfbs_se_s[cancertype_SE].dropna()#;print(cancertype,len(y))
            selected_tfs = x.sort_values(ascending=False).index[:tf_thre]
       
        
        region_count = []
        for motif_name in selected_tfs[:]:
            rp_files = glob.glob('{}/{}_ATAC_overlap_{}_SE_overlap_{}_*.tsv'.format(rp_dir,cancertype,cancertype_SE,motif_name))
            rp_files = [i for i in rp_files if not re.search('~',i)]
            assert len(rp_files) ==1

            # try:
            if 1:
                rp_df = pd.read_csv(rp_files[0],sep='\t',index_col=0,header=None)
                rp_df.columns=['RP']
                rp_df = rp_df.sort_values(by='RP',ascending=False)
                
                # ==== get genes with RP>1, read the gene clinical info
                # df = clinical_df.loc[rp_df.index[:int(rp_df[rp_df.RP>0].shape[0]*.4)]].dropna();print(cancertype,df.shape)
                df = clinical_df.loc[rp_df[rp_df.RP>rp_thre].index].dropna();print(cancertype,df.shape)
                df = df.sort_values(by='log rank p',ascending=True)
                df['RP'] = rp_df.loc[df.index].RP
                df.to_csv(outdir+os.sep+'_csv/{}_{}_targets.csv'.format(cancertype,motif_name))
                
                region_count = np.append(region_count,df.index)
                out_df = add_logrank_info(out_df,df,motif_name)
                pvalues = df['log rank p'].values
                box_vals.append(pvalues)
                # box_vals.append(-1*np.log10(pvalues))
                xticklabels.append(motif_name)

            # except:
                # print('error',motif_name)
                
        # == genes targeted by multiple TFs
        for count_thre in np.arange(2,5):
            motif_name = 'CoTarget_{}'.format(int(count_thre))
            count = dict(Counter(region_count))
            count_key = [i for i in count.keys() if count[i]>=count_thre]
            
        #     # == save box values
            df = clinical_df.loc[count_key]
            df.to_csv(outdir+os.sep+'_csv/{}_{}_targets.csv'.format(cancertype,motif_name))
        #     # pvalues = df[df['treat time']<df['ctrl time']]['log rank p'].values
            pvalues = df['log rank p'].values
            box_vals.append(pvalues)
            # box_vals.append(-1*np.log10(pvalues))
            xticklabels.append('CoTarget$≥{}$'.format(int(count_thre)))
            
            try:
                out_df = add_logrank_info(out_df,df,motif_name,cp_fig=False)
            except:
                print(data_type,cancertype,motif_name)
        
        
        out_df.to_csv('{}/{}_TF_targets_RP_GT1.csv'.format(outdir,cancertype))
        
        figname = '{}/top{}_{}_clinical_summary.pdf'.format(outdir,data_type,cancertype)
        # box_plot(box_vals,xticklabels,figname,cancertype,data_type)    
        stack_bar(box_vals,xticklabels,figname,cancertype,data_type)    
        
    



            

  


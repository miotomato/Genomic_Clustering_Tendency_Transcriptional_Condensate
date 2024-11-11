import os,sys,argparse,glob,re,bisect
import numpy as np
import pandas as pd
from collections import Counter
import operator
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
#matplotlib.rcParams['agg.path.chunksize'] = 10000
matplotlib.rcParams['font.size']=16
import seaborn as sns
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid' : False})
import scipy
import scipy.optimize
sns.set_style("ticks")
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
from scipy import stats
from scipy.stats import gamma
from sklearn.linear_model import LinearRegression



hg38_chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',
             'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',
             'chr18','chr19','chr20','chr21','chr22','chrX','chrY']




def regression_plot(df,sub_cols,basename,outdir):
    
    X = df[sub_cols]
    Y = df['TFBS_CP']
    reg=LinearRegression()     #initiating linearregression
    reg.fit(X,Y)
    Y_pred = reg.predict(X)
    Intercept=reg.intercept_
    Coefficients=reg.coef_
    print("Coefficients: ", Coefficients)
    print("Intercept: ", Intercept)
    out_df_coef = pd.DataFrame()
    out_df_coef.loc['Intercept',basename] = Intercept
    for ii in np.arange(len(sub_cols)):
        out_df_coef.loc['Coef {}'.format(sub_cols[ii]),basename] = Coefficients[ii]
    
    out_df_data = pd.DataFrame()
    out_df_data = pd.concat([Y,pd.DataFrame(Y_pred,index = Y.index,columns = ['predicted TFBS CP'])],axis=1)
    # from sklearn.linear_model import Lasso
    # lasso = Lasso(alpha=.1, max_iter=10000)
    # lasso.fit(X,Y)
    # Coefficients=lasso.coef_
    # print("Coefficients: ", Coefficients)
    # # Predict the test data
    # y_pred = lasso.predict(X)
    
    
    plt.figure(figsize=(2.6,2.6))
    x,y = Y,Y_pred
    # x,y = Y,df['TFBS_SE']
    plt.scatter(x,y,c='k',s=3)
    # ==== linear regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)       
    x_sort = np.sort(x)
    plt.plot(x_sort,x_sort*slope+intercept,c = 'grey',ls='--',lw=.9)
    plt.text(.85,1.05,'$R$ = {:.2f}'.format(r_value),fontsize=12,transform=plt.axes().transAxes)
    plt.xlabel('True TFBS CP')
    plt.ylabel('Predicted TFBS CP')
    # plt.title(cellType)
    # plt.axhline(y=0,color='k',lw=.5,ls='--')
    # plt.axvline(x=0,color='k',lw=.5,ls='--')
    # figname = '{}/{}_vs_{}_{}.pdf'.format(outdir,xlabel,ylabel,cellType)
    plt.savefig('{}/{}.pdf'.format(outdir,basename),bbox_inches='tight',pad_inches=0.02,transparent=True)
    plt.show()
    plt.close()
    out_df_coef.to_csv('{}/{}.csv'.format(outdir,basename))
    out_df_data.to_csv('{}/{}_data.csv'.format(outdir,basename))
    






indir = 'f1_CP_AICAP_expr_correlation'
outdir = 'f3_CP_regression'
os.makedirs(outdir,exist_ok=True)


df = pd.read_csv('{}/data_all.csv'.format(indir),index_col=0)
cellTypes = list(set([ii.split('_')[0] for ii in df.index]))
columns = ['TFMS_CP', 'TFMS_alpha', 'TFBS_CP', 'TFBS_ZCP', 'TFBS_CP_avgRank',
       'TFBS_SE', 'Expr_TPM', 'Expr_log2_TPM_Plus1', 'negLog2_IDR']


basename = 'all_data_regression_with_SE'
sub_cols = ['negLog2_IDR','TFMS_CP','Expr_log2_TPM_Plus1','TFBS_SE']
regression_plot(df,sub_cols,basename,outdir)

basename = 'all_data_regression_without_SE'
sub_cols = ['negLog2_IDR','TFMS_CP','Expr_log2_TPM_Plus1']
regression_plot(df,sub_cols,basename,outdir)

basename = 'data_with_SE_regression_with_SE'
df = df[df['TFBS_SE']!=0]
sub_cols = ['negLog2_IDR','TFMS_CP','Expr_log2_TPM_Plus1','TFBS_SE']
regression_plot(df,sub_cols,basename,outdir)



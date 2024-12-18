import os
import numpy as np
import pandas as pd
import json
import re
import bisect
from lifelines.statistics import logrank_test
from lifelines import KaplanMeierFitter

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size'] = 16
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"

import seaborn as sns
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid': False})
sns.set_style("ticks")

# === Function Definitions ===

def fdr_adj_p(pvalues, p_index):
    """Calculate FDR-adjusted p-values."""
    df = pvalues.to_frame()
    n, k = len(pvalues), len(pvalues)
    minimum = 1

    for index in df.sort_values(by=p_index, ascending=False).index:
        pvalue = df.loc[index, p_index]
        fdr = n * pvalue / k
        minimum = min(minimum, fdr)
        df.loc[index, 'fdr'] = minimum
        k -= 1

    return df['fdr']


def return_survival_df(region_df, region, cut_off=0.5):
    """Prepare survival DataFrame using top/bottom cutoff for openness scores."""
    df = pd.DataFrame(index=region_df.index)
    df['time_max'] = region_df[['days_to_death', 'days_to_last_follow_up']].fillna(0).astype(float).max(axis=1)
    df['time_sum'] = region_df[['days_to_death', 'days_to_last_follow_up']].fillna(0).astype(float).sum(axis=1)
    df['time'] = df['time_max']
    df = df[df['time'] > 0]

    # Use top and bottom cutoffs
    a_index = df.index[:max(int(len(df.index) * cut_off), 0)]
    b_index = df.index[-max(int(len(df.index) * cut_off), 0):]
    df = df.loc[a_index.union(b_index)]

    df.loc[a_index, 'group'] = 'treat'
    df.loc[b_index, 'group'] = 'ctrl'
    df['status'] = region_df.loc[df.index]['vital_status']
    df['score'] = region_df.loc[df.index][region]
    df.loc[df['status'] == 'Dead', 'status'] = 1
    df.loc[df['status'] == 'Alive', 'status'] = 0

    return df


def survival_for_two(df, treat, ctrl, legends, title, figname):
    """Perform log-rank test and plot survival curves for two groups."""
    ix = df['group'] == treat
    t1, e1 = df.loc[ix]['time'], df.loc[ix]['status']
    t2, e2 = df.loc[~ix]['time'], df.loc[~ix]['status']

    results = logrank_test(t1, t2, event_observed_A=e1, event_observed_B=e2)
    pvalue = results.p_value

    if pvalue < 0.001:
        plt.figure(figsize=(2.6, 2.6))
        ax = plt.subplot(111)

        kmf_treat = KaplanMeierFitter()
        kmf_treat.fit(t1, e1).plot(ax=ax, show_censors=True, label=legends[0],
                                   censor_styles={'ms': 12, 'marker': '+'}, ci_show=False, color='red', ls='-')

        kmf_ctrl = KaplanMeierFitter()
        kmf_ctrl.fit(t2, e2).plot(ax=ax, show_censors=True, label=legends[1],
                                  censor_styles={'ms': 12, 'marker': '+'}, ci_show=False, color='k', ls='--')

        ax.legend(loc='lower left', borderaxespad=0, handletextpad=0.2, labelspacing=0.2, handlelength=1, frameon=False)

        plt.text(df['time'].max() * 0.45, 0.45, f'p={pvalue:.1e}', fontsize=16, ha='center')
        plt.ylim([0.0, 1.0])
        plt.title(title, fontsize=16)
        plt.xlabel('Days', fontsize=16)
        plt.ylabel('Survival probability', fontsize=16)
        plt.savefig(figname, bbox_inches='tight', pad_inches=0.1, dpi=600, transparent=True)
        plt.close()

    return results

# === Main Processing ===

# Output directories
outdir = 'f2_caseID_each_peak_vs_clinical'
os.makedirs(outdir, exist_ok=True)
os.makedirs(f"{outdir}/fig", exist_ok=True)

# Load data mappings
name_match = pd.read_excel('data/TCGA/TCGA-ATAC_SE_cancerType_match.xlsx', index_col=0).dropna()
filtered_df = pd.read_excel('data/TCGA/TCGA-ATAC_clustered_samples.xlsx', index_col=0)

# Process each cancer type
for cancertype in ['BRCA', 'COAD']:
    cancertype_SE = name_match.loc[cancertype, 'SE']
    filtered_id = filtered_df[filtered_df['cohort'] == cancertype]['case_id']

    # Load TCGA peak data
    sig_file = 'data/TCGA/mynorm_TCGA-ATAC_PanCan_Log2_QuantileNorm_Counts_plus5.caseID.avg.txt'
    sig_df = pd.read_csv(sig_file, sep='\t', index_col=3).iloc[:, 6:]
    sig_case_id = [col.split('_')[1] for col in sig_df.columns]
    sig_df.columns = sig_case_id
    sig_case_id = [case for case in sig_case_id if case in filtered_id.values]

    # Load clinical data
    clinical_file = f"data/TCGA/clinical.project-TCGA-{cancertype}.2022-03-20.json"
    with open(clinical_file) as clinical_inf:
        clinical_list = json.load(clinical_inf)

    # Extract clinical information
    clinical_df = pd.DataFrame()
    for clinical_case in clinical_list:
        case_id = clinical_case['case_id']
        if case_id in sig_case_id:
            try:
                clinical_df.loc[case_id, 'vital_status'] = clinical_case['demographic']['vital_status']
                if clinical_df.loc[case_id, 'vital_status'] == 'Alive':
                    clinical_df.loc[case_id, 'days_to_last_follow_up'] = clinical_case['diagnoses'][0]['days_to_last_follow_up']
                elif clinical_df.loc[case_id, 'vital_status'] == 'Dead':
                    clinical_df.loc[case_id, 'days_to_death'] = clinical_case['demographic']['days_to_death']
            except Exception as e:
                print(f"Error processing case {case_id}: {e}")

    clinical_df.to_csv(f"{outdir}/{cancertype}_clinical_info.csv")

    # Analyze each ATAC-seq peak
    region_outdf = pd.DataFrame()
    for region in sig_df.index:
        sig_score = sig_df.loc[region][clinical_df.index]
        region_df = pd.concat([clinical_df, sig_score], axis=1).sort_values(by=[region], ascending=False)

        df = return_survival_df(region_df, region, cut_off=0.5)
        figname = f"{outdir}/fig/{cancertype}_{region}_survival_50th.pdf"
        legends = ['more open', 'less open']
        results = survival_for_two(df, 'treat', 'ctrl', legends, region, figname)

        if results.p_value < 0.001:
            df.to_csv(f"{outdir}/fig/{cancertype}_{region}_survival_50th.csv")

        # Save region statistics
        ix = df['group'] == 'treat'
        region_outdf.loc[region, 'treat avg atac-score'] = df.loc[ix, 'score'].mean().round(2)
        region_outdf.loc[region, 'ctrl avg atac-score'] = df.loc[~ix, 'score'].mean().round(2)
        region_outdf.loc[region, 'treat time'] = df.loc[ix, 'time'].mean().round(2)
        region_outdf.loc[region, 'ctrl time'] = df.loc[~ix, 'time'].mean().round(2)
        region_outdf.loc[region, 'log rank p'] = results.p_value

    region_outdf['fdr'] = fdr_adj_p(region_outdf['log rank p'], 'log rank p')
    region_outdf.sort_values(by='log rank p', inplace=True)
    region_outdf.to_csv(f"{outdir}/{cancertype}_logrank_info.csv")

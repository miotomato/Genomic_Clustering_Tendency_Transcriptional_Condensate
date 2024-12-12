import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import warnings

# == Configure Plotting Parameters
warnings.filterwarnings("ignore")
matplotlib.rcParams['font.size'] = 13
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"

sns.set(font_scale=1)
sns.set_style("whitegrid", {'axes.grid': False})
sns.set_style("ticks")

# == Directories
indir = '../../f1_TF_cluster_potential/f3_clustered_TFBS/f4_cobinding_TFBS_venn3'
outdir = 'f3_DCI_overlap_coBinding_TFBS'
os.makedirs(outdir, exist_ok=True)

# == Load Selected Factors for Each Cell Type
selected_factors = {}
tfbs_cp_dir = '../../f1_TF_cluster_potential/f2_cor_CP_SE_AICAP/f9_per_CT_TFBS_CP_cor_zscore_CP_with_motif_SE/TFBS_CP/'

for ct in ['MCF-7', 'HCT-116']:
    # Load TFBS CP ranking data for each cell type
    df = pd.read_csv(f"{tfbs_cp_dir}/_CP_TFBS_nonBlackList_vs_TFMS_{ct}.csv", index_col=0)
    selected_factors[f"{ct} top_TFBSCP"] = df['TFBS CP rank'].sort_values().iloc[:3].index
    selected_factors[f"{ct} top_zscored_TFBSCP"] = df['avg rank'].sort_values().iloc[:3].index

# == DCI Directory and Parameters
dci_dir = '../../../f11_TF_condensates_KS_test/f3_public_data/f1_hct116_hic_RAD21_auxin/f1_bart3d/'
genomicDistances = [100000, 200000, 500000]
reps = ['rep1', 'rep2', 'all_reps']

# == Loop Over Treatment Flags, Cell Types, and Factor Types
for treat_flag in ['percentile_T_ExtendMerge']:
    for ct in ['MCF-7', 'HCT-116']:
        for factorType in ['top_TFBSCP', 'top_zscored_TFBSCP']:
            # Create output subdirectory for the current cell type and factor type
            subdir = f"{ct}_{factorType}"
            os.makedirs(f"{outdir}/{subdir}", exist_ok=True)

            # Get selected TF factors and genomic annotations
            factors = selected_factors[f"{ct} {factorType}"]
            prenames = np.append(factors, ['hg38_exons', 'hg38_introns', 'hg38_4k_promoter_geneID'])

            # Load the cobinding data
            cobinding_file = f"{indir}/{subdir}/{treat_flag}_{ct}_cobinding.csv"
            outdf = pd.read_csv(cobinding_file)

            # Create data subsets for each combination of factors
            dfs = [
                outdf,
                outdf[(outdf[prenames[0]] != 0)],
                outdf[(outdf[prenames[1]] != 0)],
                outdf[(outdf[prenames[2]] != 0)],
                outdf[(outdf[prenames[0]] != 0) & (outdf[prenames[1]] != 0)],
                outdf[(outdf[prenames[0]] != 0) & (outdf[prenames[1]] != 0) & (outdf[prenames[2]] != 0)],
            ]

            xticklabels = [
                'Union',
                '-'.join([prenames[0]]),
                '-'.join([prenames[1]]),
                '-'.join([prenames[2]]),
                '-'.join(prenames[:2]),
                '-'.join(prenames[:3]),
            ]

            # Save BED files for each subset and perform DCI overlap analysis
            for ii, label in enumerate(xticklabels):
                outbed = f"{outdir}/{subdir}/{treat_flag}_{ct}_{label}.bed"
                dfs[ii].to_csv(outbed, index=False, sep='\t', header=None)

                # Intersect with DCI data for different genomic distances and replicates
                for genomicDis in genomicDistances:
                    for rep in reps:
                        dci_file = f"{dci_dir}/RAD21_6hr_auxin_over_NT_{rep}_dis{int(genomicDis / 1000)}k_differential_score.bed"
                        overlapped_bed = f"{outdir}/{subdir}/DCI_RAD21_KO_{rep}_dis{int(genomicDis / 1000)}k_overlapped_{treat_flag}_{ct}_{label}.bed"
                        commandLine = f"bedtools intersect \\\n-a {dci_file} \\\n-b {outbed} \\\n-wa -u > {overlapped_bed}\n"
                        print(commandLine)
                        os.system(commandLine)

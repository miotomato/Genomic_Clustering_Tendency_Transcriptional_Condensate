import os
import glob
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# == Plotting Configuration
plt.rcParams['font.size'] = 14
plt.rcParams["font.sans-serif"] = ["Arial"]
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid': False})
sns.set_style("ticks")

# == Constants
HG38_CHROMS = [
    'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
    'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
    'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY'
]

# == Load Chromosome Sizes
chrom_size = pd.read_csv(
    'data/Genome/ucsc/hg38/hg38.chrom.sizes',
    sep='\t',
    index_col=0,
    header=None
)
chrom_size.columns = ['len']

# == Output Directories
outdir = 'CP_TFMS_vs_random_plot_by_chrom'
os.makedirs(f"{outdir}/_fig", exist_ok=True)
os.makedirs(f"{outdir}/_csv", exist_ok=True)

# == Process TF Motif Files
infiles = glob.glob('data_TFMS_jarspar/*')
infiles = [i for i in infiles if not re.search('~', i)]
infiles.sort()

# Initialize DataFrame to store results
df_cp = pd.DataFrame()

for infile in infiles:
    try:
        print(f"Processing: {infile}")
        outname = os.path.basename(infile).split('.tsv')[0]
        df_all = pd.read_csv(infile, sep='\t', header=None)

        for chrom in HG38_CHROMS:
            # Filter data for the current chromosome
            df = df_all[df_all[0] == chrom]

            # Skip if there are no entries for this chromosome
            if df.empty:
                continue

            # Calculate distances and their log-transformed values
            col = df.columns[-1]
            df['dis'] = np.abs(df[col])
            values_log = np.log10(df['dis'].clip(lower=1))

            # Generate random samples from an exponential distribution
            nums = df.shape[0]
            mean_val = chrom_size.loc[chrom, 'len'] / (nums + 1)
            sample_vals = np.random.exponential(mean_val, nums)
            sample_vals_log = np.log10(sample_vals)

            # Perform KS tests (less and greater alternatives)
            sl, pl = stats.ks_2samp(sample_vals_log, values_log, alternative='less')
            sg, pg = stats.ks_2samp(sample_vals_log, values_log, alternative='greater')

            # Record the larger KS statistic with the appropriate sign
            if sl > sg:
                df_cp.loc[outname, f'CP {chrom}'] = sl
            else:
                df_cp.loc[outname, f'CP {chrom}'] = -sg

    except Exception as e:
        print(f"Error processing file {infile}: {e}")

# == Save Results
df_cp.index.name = 'TF'
df_cp.to_csv(f"{outdir}/_csv/data_CP_TFMS_vs_random_by_chrom.csv")

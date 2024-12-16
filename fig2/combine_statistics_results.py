
import os
import glob
import re
import pandas as pd
import numpy as np
from scipy import stats
import warnings

# Suppress warnings for cleaner output
warnings.filterwarnings("ignore")

def return_final_ks_statistics(df):
    """
    Computes signed KS statistics based on the 'less' and 'greater' tests.

    Parameters:
    - df: DataFrame containing KS statistics.

    Returns:
    - df: Updated DataFrame with signed KS statistics.
    """
    for ii in df.index:
        if df.loc[ii, 'log10-dis ks_2samp-s less'] > df.loc[ii, 'log10-dis ks_2samp-s greater']:
            df.loc[ii, 'log10-dis ks_2samp-s signed'] = df.loc[ii, 'log10-dis ks_2samp-s less']
            df.loc[ii, 'log10-dis ks_2samp-p signed'] = df.loc[ii, 'log10-dis ks_2samp-p less']
        else:
            df.loc[ii, 'log10-dis ks_2samp-s signed'] = -1 * df.loc[ii, 'log10-dis ks_2samp-s greater']
            df.loc[ii, 'log10-dis ks_2samp-p signed'] = df.loc[ii, 'log10-dis ks_2samp-p greater']
    return df

# ==== Main script starts here

# Directories containing the input data
indirs = [
    'CP_TFBS_nonBlackList_vs_TFMS',
    'CP_TFBS_nonBlackList_overlap_motif_vs_TFMS',
    'CP_TFBS_nonBlackList_NOT_overlap_motif_vs_TFMS'
]

# Load QC data
df_qc = pd.read_csv('data/cistrome/cistrome2019_selected_QC.csv', index_col=0)

print(f'Processing directories: {indirs}')

# Iterate through each directory
for indir in indirs:
    # Collect all CSV files in the current directory
    infiles = glob.glob(f'{indir}/_csv_cp/*csv')
    infiles.sort()

    # DataFrame to store combined statistics
    df_out = pd.DataFrame()

    # Process each input file
    for infile in infiles:
        print(f'Processing file: {infile}')
        basename = os.path.basename(infile).split('.csv')[0]

        # Skip files that don't match the expected naming pattern
        if len(basename.split('_')) != 3:
            continue

        # Extract the Cistrome ID
        dcID = int(basename.split('_')[-1])
        PeaksFCAbove10 = df_qc.loc[dcID].PeaksFoldChangeAbove10

        # Read the input file
        df = pd.read_csv(infile, index_col=0)

        # Compute signed KS statistics
        df = return_final_ks_statistics(df)

        # Store the median values and the standard deviation of KS statistics
        s_values = df['log10-dis ks_2samp-s signed']
        df_out[basename] = df.median()
        df_out.loc['log10-dis ks_2samp-s signed std', basename] = s_values.std()
        df_out.loc['PeaksFCAbove10', basename] = PeaksFCAbove10

    # Transpose the DataFrame for easier analysis
    df_out = df_out.T

    # Filter rows based on conditions
    df_out = df_out[(df_out['PeaksFCAbove10'] > 500) & (df_out['#TFMS'] > 0)]

    # Save the combined statistics to a CSV file
    output_csv_path = f'{indir}/per_data_{indir}.csv'
    df_out.to_csv(output_csv_path)
    print(f'Results saved to {output_csv_path}')

    # ==== Combine results per TF
    df_out2 = pd.DataFrame(columns=df_out.columns)
    tfs = set([i.split('_')[1] for i in df_out.index])

    for tf in sorted(tfs):
        # Filter data corresponding to the current TF
        data_index = [i for i in df_out.index if re.search(f'_{tf}_', i)]
        df = df_out.loc[data_index]

        # Compute the median values for the current TF
        df_out2.loc[tf] = df.median()

    # Save the combined results per TF
    combined_tf_csv_path = f'{indir}/combined_per_TF_{indir}.csv'
    df_out2.to_csv(combined_tf_csv_path)
    print(f'Combined TF results saved to {combined_tf_csv_path}')

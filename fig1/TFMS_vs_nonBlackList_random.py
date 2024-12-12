
import os
import pandas as pd
import numpy as np
from scipy import stats

# List of human chromosomes (excluding Y)
hg38_chroms = [
    'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
    'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
    'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX'
]

# Genome file path and blacklist file path
genomeFile = 'data/hg38_clean.chrom.nonChrY.sizes'
blackList = 'data/blackList_data/hg38-blacklist.v2.bed'

def bedtools_shuffle_closest(infile, outdir, outname):
    """
    Shuffles the input BED file, sorts it, and finds the closest features.
    
    Parameters:
    - infile: Path to the input BED file.
    - outdir: Directory where output files will be saved.
    - outname: Base name for output files.

    Returns:
    - sort_file: Path to the sorted BED file.
    - closest_file: Path to the closest feature BED file.
    """
    sample_file = f'{outdir}/_bed/{outname}.shuffle.bed'
    sort_file = f'{outdir}/_bed/{outname}.shuffle.sort.bed'
    closest_file = f'{outdir}/_bed/{outname}.shuffle.closest.bed'

    # Shuffle the input BED file excluding blacklist regions
    commandLine = f'bedtools shuffle -i {infile} -g {genomeFile} -excl {blackList} -seed 0 | cut -f1-3 > {sample_file}'
    os.system(commandLine)

    # Sort the shuffled BED file
    commandLine = f'bedtools sort -i {sample_file} > {sort_file}'
    os.system(commandLine)

    # Find the closest features
    commandLine = f'bedtools closest -a {sort_file} -b {sort_file} -D ref -fd -io -t first > {closest_file}'
    os.system(commandLine)

    # Remove the intermediate shuffled file
    os.remove(sample_file)

    return sort_file, closest_file

def return_distance_distribution(infile):
    """
    Calculates the distribution of distances from a BED file.

    Parameters:
    - infile: Path to the BED file containing distance information.

    Returns:
    - values: Absolute distance values.
    - values_log: Log10-transformed distance values.
    """
    df = pd.read_csv(infile, sep='\t', header=None, low_memory=False)
    df = df[df[0].isin(hg38_chroms)]
    col = df.columns[-1]
    values = np.abs(df[col])
    values_log = np.log10(values.clip(1))
    return values, values_log

def compr_sites_distribution(df_cp, infile, outdir):
    """
    Compares site distributions and generates shuffled controls.

    Parameters:
    - df_cp: DataFrame (unused in current function).
    - infile: Path to the input BED file.
    - outdir: Directory where output files will be saved.

    Prints the name of the input file, the number of entries, and the proportion of entries with distance = 1.
    """
    outname = os.path.basename(infile).split('.tsv')[0]
    df = pd.read_csv(infile, sep='\t', header=None)
    df = df[df[0].isin(hg38_chroms)]
    col = df.columns[-1]
    df['dis'] = np.abs(df[col])
    values = df['dis']
    print(f'{outname}\t{len(values)}\t{sum(values == 1)}\t{sum(values == 1) / len(values)}')

    # Generate shuffled controls and distance distribution
    sort_file, closest_file = bedtools_shuffle_closest(infile, outdir, outname)
    sample_vals, sample_vals_log = return_distance_distribution(closest_file)

import glob
import re

# Define the output directory
outdir = 'output_directory'  # Replace with the desired path

# == Check the TF motifs by processing multiple input files
infiles = glob.glob('f0_bedtools_closest/data_TFMS_jarspar_nonBlackList/*')
infiles = [i for i in infiles if not re.search('~', i)]
infiles.sort()
print(f'Input files: {infiles}')

# Initialize an empty DataFrame for storing results
df_cp = pd.DataFrame()

# Process each input file and generate distribution comparisons
for infile in infiles:
    try:
        df_cp = compr_sites_distribution(df_cp, infile, outdir)
    except Exception as e:
        print(f'Error processing file {infile}: {e}')  # Handle errors for empty or invalid files

# Save the results to a CSV file
df_cp.index.name = 'TF'
df_cp = df_cp.sort_values(by='log10-dis Wilcoxon-rank-sum-s', ascending=False)
output_csv_path = f'{outdir}/data_CP_TFMS_vs_random.csv'
df_cp.to_csv(output_csv_path)
print(f'Results saved to {output_csv_path}')


import os
import sys
import subprocess
import pandas as pd
import numpy as np
from scipy import stats

# List of human chromosomes (excluding Y)
hg38_chroms = [
    'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
    'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
    'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX'
]

# Paths to genome and blacklist files
genome_file = 'data/hg38_clean.chrom.nonChrY.sizes'
blackList = 'data/blackList_data/hg38-blacklist.v2.bed'

def run_bedtools_intersect(df, bfile, outdir, prename):
    """
    Intersects a BED file with a super-enhancer file and counts overlapping regions.
    """
    motif_file = f'{outdir}/{prename}.bed'
    df.to_csv(motif_file, index=False, sep='\t', header=None)
    
    cl = f'bedtools intersect -a {motif_file} -b {bfile} -u -wa | wc -l'
    wa_count = subprocess.check_output(cl, shell=True).decode(sys.stdout.encoding).strip()
    
    return int(wa_count), motif_file

def run_bedtools_shuffle_intersect(ii, motif_file, se_file, outdir, prename):
    """
    Shuffles a BED file, sorts it, finds the closest features, and intersects with super-enhancers.
    """
    motif_shuffle = f'{outdir}/{prename}_shuffle_{ii}.bed'
    sort_file = f'{outdir}/{prename}_shuffle_{ii}.sort.bed'
    closest_file = f'{outdir}/{prename}_shuffle_{ii}.closest.bed'

    # Shuffle the BED file excluding blacklist regions
    cl = f'bedtools shuffle -i {motif_file} -g {genome_file} -excl {blackList} | cut -f1-3 > {motif_shuffle}'
    os.system(cl)

    # Sort the shuffled BED file
    cl = f'bedtools sort -i {motif_shuffle} > {sort_file}'
    os.system(cl)

    # Find the closest features
    cl = f'bedtools closest -a {sort_file} -b {sort_file} -D ref -fd -io -t first > {closest_file}'
    os.system(cl)

    # Intersect with super-enhancer file
    cl = f'bedtools intersect -a {motif_shuffle} -b {se_file} -u -wa | wc -l'
    wa_count = subprocess.check_output(cl, shell=True).decode(sys.stdout.encoding).strip()

    # Read distances
    df_closest = pd.read_csv(closest_file, sep='\t', header=None, low_memory=False)
    df_closest = df_closest[df_closest[0].isin(hg38_chroms)]
    col = df_closest.columns[-1]
    values_closest = np.abs(df_closest[col])
    values_log_closest = np.log10(values_closest.clip(1))

    # Clean up intermediate files
    os.remove(motif_shuffle)
    os.remove(sort_file)
    os.remove(closest_file)

    return int(wa_count), values_closest, values_log_closest

def main(infile):
    """
    Main function to compute enrichment of TFMS at super-enhancers and save results.
    """
    outdir = 'CP_TFMS_vs_random_nonBlackList'
    se_file = 'data/SE_hg38/all_hg38_SE.bed'  # Update this path as needed

    # Ensure output directory exists
    os.makedirs(f'{outdir}/_csv', exist_ok=True)

    # Initialize result DataFrame
    df_out_tmp = pd.DataFrame()

    prename = os.path.basename(infile).split('_Exponential_Pvalue')[0]

    # Read input BED file
    df = pd.read_csv(infile)
    total = df.shape[0]

    # Count total TFMS overlapped with SE
    total_overlapped, motif_file = run_bedtools_intersect(df, se_file, outdir, prename)
    true_values = df['dis']
    values_log = np.log10(true_values.clip(1))

    # Perform 100 random shuffles and calculate statistics
    for ii in range(100):
        random_overlapped, sample_vals, sample_vals_log = run_bedtools_shuffle_intersect(ii, motif_file, se_file, outdir, prename)

        # Fisher's exact test
        s, p = stats.fisher_exact([[total_overlapped, total - total_overlapped],
                                   [random_overlapped, total - random_overlapped]])
        df_out_tmp.loc[ii, 'random total'] = total
        df_out_tmp.loc[ii, 'random SE overlapped'] = random_overlapped
        df_out_tmp.loc[ii, 'enrich-at-SE-fisher-exact-s'] = round(s, 2)
        df_out_tmp.loc[ii, 'enrich-at-SE-fisher-exact-p'] = '{:.2e}'.format(p)

    # Add total counts to the DataFrame
    df_out_tmp['motif total'] = total
    df_out_tmp['motif SE overlapped'] = total_overlapped

    # Save the result DataFrame to a CSV file
    output_csv_path = f'{outdir}/data_TFMS_enrich_at_SE.csv'
    df_out_tmp.to_csv(output_csv_path)
    print(f'Results saved to {output_csv_path}')


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'infile', metavar = '<file>')
    
    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main(args.infile)

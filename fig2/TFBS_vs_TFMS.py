
import os
import sys
import argparse
import glob
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# Set plotting styles
plt.rcParams['font.size'] = 14
sns.set(font_scale=1.1)
sns.set_style("whitegrid", {'axes.grid': False})
sns.set_style("ticks")
plt.rcParams["font.sans-serif"] = ["Arial"]

# Import custom modules
from module import compr_cluster_potential, enrichment_odds_ratio


def main(dcID):
    """
    Main function to compare TFBS types with TFMS, compute cluster potential, and enrichment.

    Parameters:
    - dcID: Cistrome ID used to identify the factor and cell type.
    """

    # Define TFBS types to process
    tfbs_types = [
        'TFBS_nonBlackList',
        'TFBS_nonBlackList_overlap_motif',
        'TFBS_nonBlackList_NOT_overlap_motif'
    ]

    # Loop through each TFBS type
    for tfbs_type in tfbs_types:
        outdir = f'CP_{tfbs_type}_vs_TFMS'
        os.makedirs(f'{outdir}/_fig', exist_ok=True)
        os.makedirs(f'{outdir}/_csv_sample', exist_ok=True)
        os.makedirs(f'{outdir}/_csv_cp', exist_ok=True)

        # ==== Load cistrome information
        df_qc = pd.read_csv('../../../f12_KS_test_Rename/data/cistrome/cistrome2019_selected_QC.csv', index_col=0)
        factor = df_qc.loc[dcID].Factor
        celltype = df_qc.loc[dcID].Cell_line

        # ==== Define paths to peak, SE, and motif files
        peak_dir = f'../f0_bedtools_closest/data_{tfbs_type}'
        peak_file = f'{peak_dir}/{celltype}_{factor}_{dcID}.tsv'

        celltype_SE_dir = '../../../f12_KS_test_Rename/data/SE_hg38'
        celltype_SE_file = f'{celltype_SE_dir}/{celltype}.bed'

        motif_dir = "../../data/motif/motif_fimo_jaspar_mid_nonBlackList"
        motif_files = glob.glob(f'{motif_dir}/{factor}_*')
        motif_files = [i for i in motif_files if not re.search('~', i)]
        assert len(motif_files) == 1, "Error: Multiple or no motif files found."
        motif_file = motif_files[0]

        # ==== Define output name and labels
        outname = f'{celltype}_{factor}_{dcID}'
        labels = [tfbs_type, 'TFMS']
        colors = ['tab:purple', 'tab:green']

        # ==== Initialize an empty DataFrame for results
        df_out = pd.DataFrame()

        # ==== Perform 100 iterations of comparison and enrichment calculation
        for ii in range(100):
            outname_ii = f'{outname}_sample{ii}'

            # Compute cluster potential
            df_out, sort_file = compr_cluster_potential(df_out, peak_file, motif_file, outdir, outname_ii, labels, colors)

            # Compute enrichment odds ratio
            df_out = enrichment_odds_ratio(df_out, peak_file, sort_file, celltype_SE_file, outname_ii, labels)

            # Remove temporary sorted file
            os.remove(sort_file)

        # ==== Save the final results to a CSV file
        df_out.index.name = 'TF'
        output_csv_path = f'{outdir}/_csv_cp/{outname}.csv'
        df_out.to_csv(output_csv_path)
        print(f'Results saved to {output_csv_path}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compare TFBS types with TFMS and compute enrichment.')
    parser.add_argument('-i', '--id', type=int, dest='id', required=True, help='Cistrome ID')

    args = parser.parse_args()

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)

    main(args.id)

import os
import pandas as pd
import numpy as np

# List of TFBS types to analyze
tfbs_prenames = ['TFBS_nonBlackList', 'TFBS_nonBlackList_overlap_motif', 'TFBS_nonBlackList_NOT_overlap_motif']
# List of percentile thresholds to use for distance filtering
percentile_thres = ['percentile-5', 'percentile-1']

# Iterate over the first percentile threshold (percentile-5)
for percentile_thre in percentile_thres[:1]:
    # Create output directory for the current threshold
    outdir = f'f1_clustered_TFBS_{percentile_thre}'
    os.makedirs(outdir, exist_ok=True)
    
    # Iterate over the first TFBS type (TFBS_nonBlackList)
    for tfbs_prename in tfbs_prenames[:1]:
        # Construct the input file name based on the current TFBS type
        cp_prename = f'CP_{tfbs_prename}_vs_TFMS'
        input_file = f'../f1_TFMS_TFBS_CP/{cp_prename}/per_data_{cp_prename}.csv'
        
        # Read the input data into a DataFrame, dropping rows with all NaN values
        df = pd.read_csv(input_file, index_col=0)
        df = df.dropna(how='all')
        
        # Filter rows based on the '#TFBS_nonBlackList' column
        df = df[df['#TFBS_nonBlackList'] > 2000]
        
        # Process each row in the filtered DataFrame
        for index in df.index:
            # Extract cell type, transcription factor, and data ID from the index
            ct, tf, data_id = index.split('_')
            
            # Only process data for specific cell types
            if ct not in ['MCF-7', 'HCT-116', 'HeLa', 'LNCaP', 'U87', 'HepG2']:
                continue
            
            # Get the distance threshold for clustering from the DataFrame
            dis_thre = df.loc[index, f'TFMS {percentile_thre}'].astype(int)
            
            # Construct the path to the TFBS data file
            data_file = f'data_TFBS_both_side_nonBlackList/{index}.tsv'
            
            # Skip if the file doesn't exist
            if not os.path.isfile(data_file):
                print(f"File not found: {index}")
                continue
            
            # Read the TFBS data into a DataFrame
            data_df = pd.read_csv(data_file, sep='\t', header=None)
            dis_col = data_df.columns[-1]
            
            # Filter data into clustered (T) and non-clustered (C) based on the distance threshold
            data_df_T = data_df[data_df[dis_col] <= dis_thre]
            data_df_C = data_df[data_df[dis_col] > dis_thre]
            
            # Create output directory for the current cell type
            os.makedirs(f"{outdir}/{ct}", exist_ok=True)
            outname = f"{ct}_{tf}_{data_id}"
            
            # Define file paths for the output BED files
            se_file = f'data/SE_hg38/{ct}.bed'
            bed_T_file = f"{outdir}/{ct}/{outname}_T.bed"
            bed_T_SE = f"{outdir}/{ct}/{outname}_T_on_SE.bed"
            bed_C_file = f"{outdir}/{ct}/{outname}_C.bed"
            bed_C_SE = f"{outdir}/{ct}/{outname}_C_on_SE.bed"
            merge_file = f"{outdir}/{ct}/{outname}_T_ExtendMerge.bed"
            merge_SE = f"{outdir}/{ct}/{outname}_T_ExtendMerge_on_SE.bed"
            
            # Save the clustered (T) and non-clustered (C) TFBS to BED files
            data_df_T.to_csv(bed_T_file, sep='\t', header=None, index=False)
            data_df_C.to_csv(bed_C_file, sep='\t', header=None, index=False)
            
            # Merge nearby TFBS clusters and save the result
            merge_command = f"bedtools merge -d {dis_thre} -i {bed_T_file} > {merge_file}"
            print(merge_command)
            os.system(merge_command)
            
            # Find overlaps between TFBS clusters and super-enhancers (SE)
            intersect_T_SE = f"bedtools intersect -a {bed_T_file} -b {se_file} -wa -u > {bed_T_SE}"
            print(intersect_T_SE)
            os.system(intersect_T_SE)
            
            intersect_C_SE = f"bedtools intersect -a {bed_C_file} -b {se_file} -wa -u > {bed_C_SE}"
            print(intersect_C_SE)
            os.system(intersect_C_SE)
            
            intersect_merge_SE = f"bedtools intersect -a {merge_file} -b {se_file} -wa -u > {merge_SE}"
            print(intersect_merge_SE)
            os.system(intersect_merge_SE)
            
            # Update the DataFrame with the new TFBS count
            df.loc[index, '#TFBS_new'] = data_df.shape[0]
        
        # Save the updated DataFrame to a new CSV file
        output_csv = f"{outdir}/{cp_prename}_new.csv"
        df.dropna().to_csv(output_csv)

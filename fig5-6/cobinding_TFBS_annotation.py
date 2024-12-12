import os
import pandas as pd

# == Input and Output Directories
indir = 'f2_bedfiles_merged'
outdir = 'f3_coBinding_merge'
os.makedirs(outdir, exist_ok=True)

# == Load Selected Factors for Each Cell Type
selected_factors = {}
tfbs_cp_dir = '../f2_cor_CP_SE_AICAP/f9_per_CT_TFBS_CP_cor_zscore_CP_with_motif_SE/TFBS_CP/'

for ct in ['MCF-7', 'HCT-116']:
    # Load TFBS CP ranking data for each cell type
    df = pd.read_csv(f"{tfbs_cp_dir}/_CP_TFBS_nonBlackList_vs_TFMS_{ct}.csv", index_col=0)
    selected_factors[f"{ct} top_TFBSCP"] = df['TFBS CP rank'].sort_values().iloc[:3].index
    selected_factors[f"{ct} top_zscored_TFBSCP"] = df['avg rank'].sort_values().iloc[:3].index

# == Process Each Treatment Flag, Cell Type, and Factor Type
for treat_flag in ['percentile_T', 'percentile_T_ExtendMerge']:
    for ct in ['MCF-7', 'HCT-116']:
        for factorType in ['top_TFBSCP', 'top_zscored_TFBSCP']:
            factors = selected_factors[f"{ct} {factorType}"]
            subdir = f"{ct}_{factorType}"
            os.makedirs(f"{outdir}/{subdir}", exist_ok=True)

            # Create paths for concatenated, sorted, and merged BED files
            bedfiles = [f"{indir}/{ct}/{ct}_{factor}_{treat_flag}.merge.bed" for factor in factors]
            catFile = f"{outdir}/{subdir}/{treat_flag}.cat.bed"
            sortFile = f"{outdir}/{subdir}/{treat_flag}.sort.bed"
            mergeFile = f"{outdir}/{subdir}/{treat_flag}.merge.bed"

            # Step 1: Concatenate BED files
            commandLine = f"\ncat {' \\\n'.join(bedfiles)} > {catFile}"
            print(commandLine)

            # Step 2: Sort the concatenated file
            commandLine = f"bedtools sort -i {catFile} > {sortFile}"
            print(commandLine)

            # Step 3: Merge overlapping regions
            commandLine = f"bedtools merge -i {sortFile} > {mergeFile}\n"
            print(commandLine)

            # Step 4: Intersect merged file with individual factor files
            for factor in factors:
                afile = mergeFile
                bfile = f"{indir}/{ct}/{ct}_{factor}_{treat_flag}.merge.bed"
                outfile = f"{outdir}/{subdir}/{treat_flag}.merge.{factor}.overlapped.bed"
                commandLine = f"bedtools intersect -a {afile} -b {bfile} -wa -c > {outfile}"
                print(commandLine)

            # Step 5: Intersect merged file with annotation files (exons, introns, promoters)
            for prename in ['hg38_exons', 'hg38_introns', 'hg38_4k_promoter_geneID']:
                afile = mergeFile
                bfile = f"/standard/vol190/zanglab/zw5j/data/geneID_annotation/hg38/{prename}.bed"
                outfile = f"{outdir}/{subdir}/{treat_flag}.merge.{prename}.overlapped.bed"
                commandLine = f"bedtools intersect -a {afile} -b {bfile} -wa -c > {outfile}"
                print(commandLine)

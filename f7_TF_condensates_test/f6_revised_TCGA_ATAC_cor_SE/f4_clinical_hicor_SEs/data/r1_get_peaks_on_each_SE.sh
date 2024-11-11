# ==== BRCA
se_file='../../data/SE_hg38/MCF-7.bed'
atac_file='../../f3_clinical_hicor_SE_overlapped_peaks/f1_clinical_at_each_peak/f1_ATAC_overlap_SE_caseID/BRCA_ATAC_overlap_MCF-7_SE_caseID.bed'

bedtools intersect -a $se_file -b $atac_file -wa -wb > BRCA_SE_with_overlapped_peaks.bed


# ==== COAD
se_file='../../data/SE_hg38/HCT-116.bed'
atac_file='../../f3_clinical_hicor_SE_overlapped_peaks/f1_clinical_at_each_peak/f1_ATAC_overlap_SE_caseID/COAD_ATAC_overlap_HCT-116_SE_caseID.bed'

bedtools intersect -a $se_file -b $atac_file -wa -wb > COAD_SE_with_overlapped_peaks.bed

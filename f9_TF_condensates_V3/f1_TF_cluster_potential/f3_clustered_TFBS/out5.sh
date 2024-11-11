bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/MCF-7_top_TFBSCP/percentile_T_MCF-7_Union.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/MCF-7_top_TFBSCP/atac_overlapped_percentile_T_MCF-7_Union.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/MCF-7_top_TFBSCP/percentile_T_MCF-7_NR2F2.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/MCF-7_top_TFBSCP/atac_overlapped_percentile_T_MCF-7_NR2F2.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/MCF-7_top_TFBSCP/percentile_T_MCF-7_MAX.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/MCF-7_top_TFBSCP/atac_overlapped_percentile_T_MCF-7_MAX.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/MCF-7_top_TFBSCP/percentile_T_MCF-7_ERG.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/MCF-7_top_TFBSCP/atac_overlapped_percentile_T_MCF-7_ERG.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/MCF-7_top_TFBSCP/percentile_T_MCF-7_NR2F2-MAX.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/MCF-7_top_TFBSCP/atac_overlapped_percentile_T_MCF-7_NR2F2-MAX.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/MCF-7_top_TFBSCP/percentile_T_MCF-7_NR2F2-MAX-ERG.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/MCF-7_top_TFBSCP/atac_overlapped_percentile_T_MCF-7_NR2F2-MAX-ERG.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/MCF-7_top_zscored_TFBSCP/percentile_T_MCF-7_Union.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/MCF-7_top_zscored_TFBSCP/atac_overlapped_percentile_T_MCF-7_Union.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/MCF-7_top_zscored_TFBSCP/percentile_T_MCF-7_ERG.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/MCF-7_top_zscored_TFBSCP/atac_overlapped_percentile_T_MCF-7_ERG.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/MCF-7_top_zscored_TFBSCP/percentile_T_MCF-7_E2F1.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/MCF-7_top_zscored_TFBSCP/atac_overlapped_percentile_T_MCF-7_E2F1.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/MCF-7_top_zscored_TFBSCP/percentile_T_MCF-7_MYC.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/MCF-7_top_zscored_TFBSCP/atac_overlapped_percentile_T_MCF-7_MYC.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/MCF-7_top_zscored_TFBSCP/percentile_T_MCF-7_ERG-E2F1.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/MCF-7_top_zscored_TFBSCP/atac_overlapped_percentile_T_MCF-7_ERG-E2F1.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/MCF-7_top_zscored_TFBSCP/percentile_T_MCF-7_ERG-E2F1-MYC.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/MCF-7_top_zscored_TFBSCP/atac_overlapped_percentile_T_MCF-7_ERG-E2F1-MYC.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/HCT-116_top_TFBSCP/percentile_T_HCT-116_Union.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/HCT-116_top_TFBSCP/atac_overlapped_percentile_T_HCT-116_Union.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/HCT-116_top_TFBSCP/percentile_T_HCT-116_JUND.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/HCT-116_top_TFBSCP/atac_overlapped_percentile_T_HCT-116_JUND.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/HCT-116_top_TFBSCP/percentile_T_HCT-116_CEBPB.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/HCT-116_top_TFBSCP/atac_overlapped_percentile_T_HCT-116_CEBPB.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/HCT-116_top_TFBSCP/percentile_T_HCT-116_MAX.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/HCT-116_top_TFBSCP/atac_overlapped_percentile_T_HCT-116_MAX.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/HCT-116_top_TFBSCP/percentile_T_HCT-116_JUND-CEBPB.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/HCT-116_top_TFBSCP/atac_overlapped_percentile_T_HCT-116_JUND-CEBPB.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/HCT-116_top_TFBSCP/percentile_T_HCT-116_JUND-CEBPB-MAX.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/HCT-116_top_TFBSCP/atac_overlapped_percentile_T_HCT-116_JUND-CEBPB-MAX.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/HCT-116_top_zscored_TFBSCP/percentile_T_HCT-116_Union.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/HCT-116_top_zscored_TFBSCP/atac_overlapped_percentile_T_HCT-116_Union.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/HCT-116_top_zscored_TFBSCP/percentile_T_HCT-116_SRF.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/HCT-116_top_zscored_TFBSCP/atac_overlapped_percentile_T_HCT-116_SRF.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/HCT-116_top_zscored_TFBSCP/percentile_T_HCT-116_JUND.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/HCT-116_top_zscored_TFBSCP/atac_overlapped_percentile_T_HCT-116_JUND.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/HCT-116_top_zscored_TFBSCP/percentile_T_HCT-116_CEBPB.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/HCT-116_top_zscored_TFBSCP/atac_overlapped_percentile_T_HCT-116_CEBPB.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/HCT-116_top_zscored_TFBSCP/percentile_T_HCT-116_SRF-JUND.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/HCT-116_top_zscored_TFBSCP/atac_overlapped_percentile_T_HCT-116_SRF-JUND.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/HCT-116_top_zscored_TFBSCP/percentile_T_HCT-116_SRF-JUND-CEBPB.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/HCT-116_top_zscored_TFBSCP/atac_overlapped_percentile_T_HCT-116_SRF-JUND-CEBPB.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/MCF-7_top_TFBSCP/percentile_T_ExtendMerge_MCF-7_Union.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/MCF-7_top_TFBSCP/atac_overlapped_percentile_T_ExtendMerge_MCF-7_Union.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/MCF-7_top_TFBSCP/percentile_T_ExtendMerge_MCF-7_NR2F2.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/MCF-7_top_TFBSCP/atac_overlapped_percentile_T_ExtendMerge_MCF-7_NR2F2.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/MCF-7_top_TFBSCP/percentile_T_ExtendMerge_MCF-7_MAX.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/MCF-7_top_TFBSCP/atac_overlapped_percentile_T_ExtendMerge_MCF-7_MAX.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/MCF-7_top_TFBSCP/percentile_T_ExtendMerge_MCF-7_ERG.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/MCF-7_top_TFBSCP/atac_overlapped_percentile_T_ExtendMerge_MCF-7_ERG.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/MCF-7_top_TFBSCP/percentile_T_ExtendMerge_MCF-7_NR2F2-MAX.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/MCF-7_top_TFBSCP/atac_overlapped_percentile_T_ExtendMerge_MCF-7_NR2F2-MAX.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/MCF-7_top_TFBSCP/percentile_T_ExtendMerge_MCF-7_NR2F2-MAX-ERG.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/MCF-7_top_TFBSCP/atac_overlapped_percentile_T_ExtendMerge_MCF-7_NR2F2-MAX-ERG.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/MCF-7_top_zscored_TFBSCP/percentile_T_ExtendMerge_MCF-7_Union.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/MCF-7_top_zscored_TFBSCP/atac_overlapped_percentile_T_ExtendMerge_MCF-7_Union.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/MCF-7_top_zscored_TFBSCP/percentile_T_ExtendMerge_MCF-7_ERG.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/MCF-7_top_zscored_TFBSCP/atac_overlapped_percentile_T_ExtendMerge_MCF-7_ERG.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/MCF-7_top_zscored_TFBSCP/percentile_T_ExtendMerge_MCF-7_E2F1.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/MCF-7_top_zscored_TFBSCP/atac_overlapped_percentile_T_ExtendMerge_MCF-7_E2F1.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/MCF-7_top_zscored_TFBSCP/percentile_T_ExtendMerge_MCF-7_MYC.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/MCF-7_top_zscored_TFBSCP/atac_overlapped_percentile_T_ExtendMerge_MCF-7_MYC.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/MCF-7_top_zscored_TFBSCP/percentile_T_ExtendMerge_MCF-7_ERG-E2F1.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/MCF-7_top_zscored_TFBSCP/atac_overlapped_percentile_T_ExtendMerge_MCF-7_ERG-E2F1.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/MCF-7_top_zscored_TFBSCP/percentile_T_ExtendMerge_MCF-7_ERG-E2F1-MYC.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/MCF-7_top_zscored_TFBSCP/atac_overlapped_percentile_T_ExtendMerge_MCF-7_ERG-E2F1-MYC.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/HCT-116_top_TFBSCP/percentile_T_ExtendMerge_HCT-116_Union.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/HCT-116_top_TFBSCP/atac_overlapped_percentile_T_ExtendMerge_HCT-116_Union.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/HCT-116_top_TFBSCP/percentile_T_ExtendMerge_HCT-116_JUND.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/HCT-116_top_TFBSCP/atac_overlapped_percentile_T_ExtendMerge_HCT-116_JUND.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/HCT-116_top_TFBSCP/percentile_T_ExtendMerge_HCT-116_CEBPB.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/HCT-116_top_TFBSCP/atac_overlapped_percentile_T_ExtendMerge_HCT-116_CEBPB.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/HCT-116_top_TFBSCP/percentile_T_ExtendMerge_HCT-116_MAX.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/HCT-116_top_TFBSCP/atac_overlapped_percentile_T_ExtendMerge_HCT-116_MAX.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/HCT-116_top_TFBSCP/percentile_T_ExtendMerge_HCT-116_JUND-CEBPB.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/HCT-116_top_TFBSCP/atac_overlapped_percentile_T_ExtendMerge_HCT-116_JUND-CEBPB.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/HCT-116_top_TFBSCP/percentile_T_ExtendMerge_HCT-116_JUND-CEBPB-MAX.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/HCT-116_top_TFBSCP/atac_overlapped_percentile_T_ExtendMerge_HCT-116_JUND-CEBPB-MAX.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/HCT-116_top_zscored_TFBSCP/percentile_T_ExtendMerge_HCT-116_Union.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/HCT-116_top_zscored_TFBSCP/atac_overlapped_percentile_T_ExtendMerge_HCT-116_Union.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/HCT-116_top_zscored_TFBSCP/percentile_T_ExtendMerge_HCT-116_SRF.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/HCT-116_top_zscored_TFBSCP/atac_overlapped_percentile_T_ExtendMerge_HCT-116_SRF.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/HCT-116_top_zscored_TFBSCP/percentile_T_ExtendMerge_HCT-116_JUND.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/HCT-116_top_zscored_TFBSCP/atac_overlapped_percentile_T_ExtendMerge_HCT-116_JUND.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/HCT-116_top_zscored_TFBSCP/percentile_T_ExtendMerge_HCT-116_CEBPB.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/HCT-116_top_zscored_TFBSCP/atac_overlapped_percentile_T_ExtendMerge_HCT-116_CEBPB.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/HCT-116_top_zscored_TFBSCP/percentile_T_ExtendMerge_HCT-116_SRF-JUND.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/HCT-116_top_zscored_TFBSCP/atac_overlapped_percentile_T_ExtendMerge_HCT-116_SRF-JUND.bed

bedtools intersect \
-a ../../data/TCGA/tcga_atac.bed \
-b f5_atac_overlap_coBinding_TFBS/HCT-116_top_zscored_TFBSCP/percentile_T_ExtendMerge_HCT-116_SRF-JUND-CEBPB.bed \
-wa -u > f5_atac_overlap_coBinding_TFBS/HCT-116_top_zscored_TFBSCP/atac_overlapped_percentile_T_ExtendMerge_HCT-116_SRF-JUND-CEBPB.bed


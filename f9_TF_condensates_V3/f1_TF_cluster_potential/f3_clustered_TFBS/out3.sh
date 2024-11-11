
cat f2_bedfiles_merged/MCF-7/MCF-7_NR2F2_percentile_T.merge.bed \
f2_bedfiles_merged/MCF-7/MCF-7_MAX_percentile_T.merge.bed \
f2_bedfiles_merged/MCF-7/MCF-7_ERG_percentile_T.merge.bed > f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T.cat.bed

bedtools sort -i f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T.cat.bed > f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T.sort.bed
bedtools merge -i f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T.sort.bed > f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T.merge.bed

bedtools intersect -a f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T.merge.bed -b f2_bedfiles_merged/MCF-7/MCF-7_NR2F2_percentile_T.merge.bed -wa -c > f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T.merge.NR2F2.overlapped.bed
bedtools intersect -a f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T.merge.bed -b f2_bedfiles_merged/MCF-7/MCF-7_MAX_percentile_T.merge.bed -wa -c > f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T.merge.MAX.overlapped.bed
bedtools intersect -a f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T.merge.bed -b f2_bedfiles_merged/MCF-7/MCF-7_ERG_percentile_T.merge.bed -wa -c > f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T.merge.ERG.overlapped.bed
bedtools intersect -a f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_exons.bed -wa -c > f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T.merge.hg38_exons.overlapped.bed
bedtools intersect -a f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_introns.bed -wa -c > f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T.merge.hg38_introns.overlapped.bed
bedtools intersect -a f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_4k_promoter_geneID.bed -wa -c > f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T.merge.hg38_4k_promoter_geneID.overlapped.bed

cat f2_bedfiles_merged/MCF-7/MCF-7_ERG_percentile_T.merge.bed \
f2_bedfiles_merged/MCF-7/MCF-7_E2F1_percentile_T.merge.bed \
f2_bedfiles_merged/MCF-7/MCF-7_MYC_percentile_T.merge.bed > f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T.cat.bed

bedtools sort -i f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T.cat.bed > f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T.sort.bed
bedtools merge -i f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T.sort.bed > f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T.merge.bed

bedtools intersect -a f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T.merge.bed -b f2_bedfiles_merged/MCF-7/MCF-7_ERG_percentile_T.merge.bed -wa -c > f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T.merge.ERG.overlapped.bed
bedtools intersect -a f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T.merge.bed -b f2_bedfiles_merged/MCF-7/MCF-7_E2F1_percentile_T.merge.bed -wa -c > f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T.merge.E2F1.overlapped.bed
bedtools intersect -a f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T.merge.bed -b f2_bedfiles_merged/MCF-7/MCF-7_MYC_percentile_T.merge.bed -wa -c > f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T.merge.MYC.overlapped.bed
bedtools intersect -a f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_exons.bed -wa -c > f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T.merge.hg38_exons.overlapped.bed
bedtools intersect -a f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_introns.bed -wa -c > f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T.merge.hg38_introns.overlapped.bed
bedtools intersect -a f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_4k_promoter_geneID.bed -wa -c > f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T.merge.hg38_4k_promoter_geneID.overlapped.bed

cat f2_bedfiles_merged/HCT-116/HCT-116_JUND_percentile_T.merge.bed \
f2_bedfiles_merged/HCT-116/HCT-116_CEBPB_percentile_T.merge.bed \
f2_bedfiles_merged/HCT-116/HCT-116_MAX_percentile_T.merge.bed > f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T.cat.bed

bedtools sort -i f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T.cat.bed > f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T.sort.bed
bedtools merge -i f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T.sort.bed > f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T.merge.bed

bedtools intersect -a f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T.merge.bed -b f2_bedfiles_merged/HCT-116/HCT-116_JUND_percentile_T.merge.bed -wa -c > f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T.merge.JUND.overlapped.bed
bedtools intersect -a f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T.merge.bed -b f2_bedfiles_merged/HCT-116/HCT-116_CEBPB_percentile_T.merge.bed -wa -c > f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T.merge.CEBPB.overlapped.bed
bedtools intersect -a f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T.merge.bed -b f2_bedfiles_merged/HCT-116/HCT-116_MAX_percentile_T.merge.bed -wa -c > f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T.merge.MAX.overlapped.bed
bedtools intersect -a f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_exons.bed -wa -c > f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T.merge.hg38_exons.overlapped.bed
bedtools intersect -a f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_introns.bed -wa -c > f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T.merge.hg38_introns.overlapped.bed
bedtools intersect -a f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_4k_promoter_geneID.bed -wa -c > f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T.merge.hg38_4k_promoter_geneID.overlapped.bed

cat f2_bedfiles_merged/HCT-116/HCT-116_SRF_percentile_T.merge.bed \
f2_bedfiles_merged/HCT-116/HCT-116_JUND_percentile_T.merge.bed \
f2_bedfiles_merged/HCT-116/HCT-116_CEBPB_percentile_T.merge.bed > f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T.cat.bed

bedtools sort -i f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T.cat.bed > f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T.sort.bed
bedtools merge -i f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T.sort.bed > f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T.merge.bed

bedtools intersect -a f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T.merge.bed -b f2_bedfiles_merged/HCT-116/HCT-116_SRF_percentile_T.merge.bed -wa -c > f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T.merge.SRF.overlapped.bed
bedtools intersect -a f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T.merge.bed -b f2_bedfiles_merged/HCT-116/HCT-116_JUND_percentile_T.merge.bed -wa -c > f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T.merge.JUND.overlapped.bed
bedtools intersect -a f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T.merge.bed -b f2_bedfiles_merged/HCT-116/HCT-116_CEBPB_percentile_T.merge.bed -wa -c > f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T.merge.CEBPB.overlapped.bed
bedtools intersect -a f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_exons.bed -wa -c > f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T.merge.hg38_exons.overlapped.bed
bedtools intersect -a f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_introns.bed -wa -c > f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T.merge.hg38_introns.overlapped.bed
bedtools intersect -a f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_4k_promoter_geneID.bed -wa -c > f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T.merge.hg38_4k_promoter_geneID.overlapped.bed

cat f2_bedfiles_merged/MCF-7/MCF-7_NR2F2_percentile_T_ExtendMerge.merge.bed \
f2_bedfiles_merged/MCF-7/MCF-7_MAX_percentile_T_ExtendMerge.merge.bed \
f2_bedfiles_merged/MCF-7/MCF-7_ERG_percentile_T_ExtendMerge.merge.bed > f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T_ExtendMerge.cat.bed

bedtools sort -i f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T_ExtendMerge.cat.bed > f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T_ExtendMerge.sort.bed
bedtools merge -i f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T_ExtendMerge.sort.bed > f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T_ExtendMerge.merge.bed

bedtools intersect -a f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T_ExtendMerge.merge.bed -b f2_bedfiles_merged/MCF-7/MCF-7_NR2F2_percentile_T_ExtendMerge.merge.bed -wa -c > f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T_ExtendMerge.merge.NR2F2.overlapped.bed
bedtools intersect -a f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T_ExtendMerge.merge.bed -b f2_bedfiles_merged/MCF-7/MCF-7_MAX_percentile_T_ExtendMerge.merge.bed -wa -c > f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T_ExtendMerge.merge.MAX.overlapped.bed
bedtools intersect -a f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T_ExtendMerge.merge.bed -b f2_bedfiles_merged/MCF-7/MCF-7_ERG_percentile_T_ExtendMerge.merge.bed -wa -c > f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T_ExtendMerge.merge.ERG.overlapped.bed
bedtools intersect -a f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T_ExtendMerge.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_exons.bed -wa -c > f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T_ExtendMerge.merge.hg38_exons.overlapped.bed
bedtools intersect -a f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T_ExtendMerge.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_introns.bed -wa -c > f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T_ExtendMerge.merge.hg38_introns.overlapped.bed
bedtools intersect -a f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T_ExtendMerge.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_4k_promoter_geneID.bed -wa -c > f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T_ExtendMerge.merge.hg38_4k_promoter_geneID.overlapped.bed

cat f2_bedfiles_merged/MCF-7/MCF-7_ERG_percentile_T_ExtendMerge.merge.bed \
f2_bedfiles_merged/MCF-7/MCF-7_E2F1_percentile_T_ExtendMerge.merge.bed \
f2_bedfiles_merged/MCF-7/MCF-7_MYC_percentile_T_ExtendMerge.merge.bed > f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T_ExtendMerge.cat.bed

bedtools sort -i f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T_ExtendMerge.cat.bed > f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T_ExtendMerge.sort.bed
bedtools merge -i f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T_ExtendMerge.sort.bed > f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.bed

bedtools intersect -a f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.bed -b f2_bedfiles_merged/MCF-7/MCF-7_ERG_percentile_T_ExtendMerge.merge.bed -wa -c > f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.ERG.overlapped.bed
bedtools intersect -a f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.bed -b f2_bedfiles_merged/MCF-7/MCF-7_E2F1_percentile_T_ExtendMerge.merge.bed -wa -c > f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.E2F1.overlapped.bed
bedtools intersect -a f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.bed -b f2_bedfiles_merged/MCF-7/MCF-7_MYC_percentile_T_ExtendMerge.merge.bed -wa -c > f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.MYC.overlapped.bed
bedtools intersect -a f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_exons.bed -wa -c > f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.hg38_exons.overlapped.bed
bedtools intersect -a f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_introns.bed -wa -c > f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.hg38_introns.overlapped.bed
bedtools intersect -a f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_4k_promoter_geneID.bed -wa -c > f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.hg38_4k_promoter_geneID.overlapped.bed

cat f2_bedfiles_merged/HCT-116/HCT-116_JUND_percentile_T_ExtendMerge.merge.bed \
f2_bedfiles_merged/HCT-116/HCT-116_CEBPB_percentile_T_ExtendMerge.merge.bed \
f2_bedfiles_merged/HCT-116/HCT-116_MAX_percentile_T_ExtendMerge.merge.bed > f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T_ExtendMerge.cat.bed

bedtools sort -i f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T_ExtendMerge.cat.bed > f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T_ExtendMerge.sort.bed
bedtools merge -i f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T_ExtendMerge.sort.bed > f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T_ExtendMerge.merge.bed

bedtools intersect -a f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T_ExtendMerge.merge.bed -b f2_bedfiles_merged/HCT-116/HCT-116_JUND_percentile_T_ExtendMerge.merge.bed -wa -c > f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T_ExtendMerge.merge.JUND.overlapped.bed
bedtools intersect -a f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T_ExtendMerge.merge.bed -b f2_bedfiles_merged/HCT-116/HCT-116_CEBPB_percentile_T_ExtendMerge.merge.bed -wa -c > f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T_ExtendMerge.merge.CEBPB.overlapped.bed
bedtools intersect -a f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T_ExtendMerge.merge.bed -b f2_bedfiles_merged/HCT-116/HCT-116_MAX_percentile_T_ExtendMerge.merge.bed -wa -c > f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T_ExtendMerge.merge.MAX.overlapped.bed
bedtools intersect -a f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T_ExtendMerge.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_exons.bed -wa -c > f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T_ExtendMerge.merge.hg38_exons.overlapped.bed
bedtools intersect -a f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T_ExtendMerge.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_introns.bed -wa -c > f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T_ExtendMerge.merge.hg38_introns.overlapped.bed
bedtools intersect -a f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T_ExtendMerge.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_4k_promoter_geneID.bed -wa -c > f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T_ExtendMerge.merge.hg38_4k_promoter_geneID.overlapped.bed

cat f2_bedfiles_merged/HCT-116/HCT-116_SRF_percentile_T_ExtendMerge.merge.bed \
f2_bedfiles_merged/HCT-116/HCT-116_JUND_percentile_T_ExtendMerge.merge.bed \
f2_bedfiles_merged/HCT-116/HCT-116_CEBPB_percentile_T_ExtendMerge.merge.bed > f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T_ExtendMerge.cat.bed

bedtools sort -i f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T_ExtendMerge.cat.bed > f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T_ExtendMerge.sort.bed
bedtools merge -i f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T_ExtendMerge.sort.bed > f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.bed

bedtools intersect -a f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.bed -b f2_bedfiles_merged/HCT-116/HCT-116_SRF_percentile_T_ExtendMerge.merge.bed -wa -c > f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.SRF.overlapped.bed
bedtools intersect -a f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.bed -b f2_bedfiles_merged/HCT-116/HCT-116_JUND_percentile_T_ExtendMerge.merge.bed -wa -c > f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.JUND.overlapped.bed
bedtools intersect -a f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.bed -b f2_bedfiles_merged/HCT-116/HCT-116_CEBPB_percentile_T_ExtendMerge.merge.bed -wa -c > f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.CEBPB.overlapped.bed
bedtools intersect -a f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_exons.bed -wa -c > f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.hg38_exons.overlapped.bed
bedtools intersect -a f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_introns.bed -wa -c > f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.hg38_introns.overlapped.bed
bedtools intersect -a f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_4k_promoter_geneID.bed -wa -c > f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.hg38_4k_promoter_geneID.overlapped.bed

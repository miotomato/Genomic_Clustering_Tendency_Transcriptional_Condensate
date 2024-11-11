
cat f2_bedfiles_merged/MCF-7/MCF-7_ERG_percentile_T.merge.bed \
f2_bedfiles_merged/MCF-7/MCF-7_ELK1_percentile_T.merge.bed \
f2_bedfiles_merged/MCF-7/MCF-7_FOS_percentile_T.merge.bed > f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T.cat.bed

bedtools sort -i f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T.cat.bed > f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T.sort.bed
bedtools merge -i f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T.sort.bed > f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T.merge.bed

bedtools intersect -a f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T.merge.bed -b f2_bedfiles_merged/MCF-7/MCF-7_ERG_percentile_T.merge.bed -wa -c > f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T.merge.ERG.overlapped.bed
bedtools intersect -a f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T.merge.bed -b f2_bedfiles_merged/MCF-7/MCF-7_ELK1_percentile_T.merge.bed -wa -c > f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T.merge.ELK1.overlapped.bed
bedtools intersect -a f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T.merge.bed -b f2_bedfiles_merged/MCF-7/MCF-7_FOS_percentile_T.merge.bed -wa -c > f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T.merge.FOS.overlapped.bed
bedtools intersect -a f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_exons.bed -wa -c > f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T.merge.hg38_exons.overlapped.bed
bedtools intersect -a f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_introns.bed -wa -c > f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T.merge.hg38_introns.overlapped.bed
bedtools intersect -a f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_4k_promoter_geneID.bed -wa -c > f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T.merge.hg38_4k_promoter_geneID.overlapped.bed

cat f2_bedfiles_merged/MCF-7/MCF-7_ERG_percentile_T.merge.bed \
f2_bedfiles_merged/MCF-7/MCF-7_JUND_percentile_T.merge.bed \
f2_bedfiles_merged/MCF-7/MCF-7_ZNF143_percentile_T.merge.bed > f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T.cat.bed

bedtools sort -i f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T.cat.bed > f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T.sort.bed
bedtools merge -i f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T.sort.bed > f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T.merge.bed

bedtools intersect -a f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T.merge.bed -b f2_bedfiles_merged/MCF-7/MCF-7_ERG_percentile_T.merge.bed -wa -c > f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T.merge.ERG.overlapped.bed
bedtools intersect -a f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T.merge.bed -b f2_bedfiles_merged/MCF-7/MCF-7_JUND_percentile_T.merge.bed -wa -c > f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T.merge.JUND.overlapped.bed
bedtools intersect -a f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T.merge.bed -b f2_bedfiles_merged/MCF-7/MCF-7_ZNF143_percentile_T.merge.bed -wa -c > f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T.merge.ZNF143.overlapped.bed
bedtools intersect -a f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_exons.bed -wa -c > f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T.merge.hg38_exons.overlapped.bed
bedtools intersect -a f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_introns.bed -wa -c > f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T.merge.hg38_introns.overlapped.bed
bedtools intersect -a f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_4k_promoter_geneID.bed -wa -c > f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T.merge.hg38_4k_promoter_geneID.overlapped.bed

cat f2_bedfiles_merged/HCT-116/HCT-116_ELF1_percentile_T.merge.bed \
f2_bedfiles_merged/HCT-116/HCT-116_JUND_percentile_T.merge.bed \
f2_bedfiles_merged/HCT-116/HCT-116_SRF_percentile_T.merge.bed > f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T.cat.bed

bedtools sort -i f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T.cat.bed > f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T.sort.bed
bedtools merge -i f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T.sort.bed > f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T.merge.bed

bedtools intersect -a f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T.merge.bed -b f2_bedfiles_merged/HCT-116/HCT-116_ELF1_percentile_T.merge.bed -wa -c > f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T.merge.ELF1.overlapped.bed
bedtools intersect -a f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T.merge.bed -b f2_bedfiles_merged/HCT-116/HCT-116_JUND_percentile_T.merge.bed -wa -c > f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T.merge.JUND.overlapped.bed
bedtools intersect -a f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T.merge.bed -b f2_bedfiles_merged/HCT-116/HCT-116_SRF_percentile_T.merge.bed -wa -c > f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T.merge.SRF.overlapped.bed
bedtools intersect -a f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_exons.bed -wa -c > f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T.merge.hg38_exons.overlapped.bed
bedtools intersect -a f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_introns.bed -wa -c > f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T.merge.hg38_introns.overlapped.bed
bedtools intersect -a f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_4k_promoter_geneID.bed -wa -c > f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T.merge.hg38_4k_promoter_geneID.overlapped.bed

cat f2_bedfiles_merged/HCT-116/HCT-116_JUND_percentile_T.merge.bed \
f2_bedfiles_merged/HCT-116/HCT-116_CEBPB_percentile_T.merge.bed \
f2_bedfiles_merged/HCT-116/HCT-116_YY1_percentile_T.merge.bed > f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T.cat.bed

bedtools sort -i f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T.cat.bed > f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T.sort.bed
bedtools merge -i f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T.sort.bed > f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T.merge.bed

bedtools intersect -a f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T.merge.bed -b f2_bedfiles_merged/HCT-116/HCT-116_JUND_percentile_T.merge.bed -wa -c > f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T.merge.JUND.overlapped.bed
bedtools intersect -a f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T.merge.bed -b f2_bedfiles_merged/HCT-116/HCT-116_CEBPB_percentile_T.merge.bed -wa -c > f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T.merge.CEBPB.overlapped.bed
bedtools intersect -a f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T.merge.bed -b f2_bedfiles_merged/HCT-116/HCT-116_YY1_percentile_T.merge.bed -wa -c > f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T.merge.YY1.overlapped.bed
bedtools intersect -a f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_exons.bed -wa -c > f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T.merge.hg38_exons.overlapped.bed
bedtools intersect -a f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_introns.bed -wa -c > f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T.merge.hg38_introns.overlapped.bed
bedtools intersect -a f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_4k_promoter_geneID.bed -wa -c > f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T.merge.hg38_4k_promoter_geneID.overlapped.bed

cat f2_bedfiles_merged/MCF-7/MCF-7_ERG_percentile_T_ExtendMerge.merge.bed \
f2_bedfiles_merged/MCF-7/MCF-7_ELK1_percentile_T_ExtendMerge.merge.bed \
f2_bedfiles_merged/MCF-7/MCF-7_FOS_percentile_T_ExtendMerge.merge.bed > f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T_ExtendMerge.cat.bed

bedtools sort -i f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T_ExtendMerge.cat.bed > f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T_ExtendMerge.sort.bed
bedtools merge -i f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T_ExtendMerge.sort.bed > f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T_ExtendMerge.merge.bed

bedtools intersect -a f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T_ExtendMerge.merge.bed -b f2_bedfiles_merged/MCF-7/MCF-7_ERG_percentile_T_ExtendMerge.merge.bed -wa -c > f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T_ExtendMerge.merge.ERG.overlapped.bed
bedtools intersect -a f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T_ExtendMerge.merge.bed -b f2_bedfiles_merged/MCF-7/MCF-7_ELK1_percentile_T_ExtendMerge.merge.bed -wa -c > f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T_ExtendMerge.merge.ELK1.overlapped.bed
bedtools intersect -a f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T_ExtendMerge.merge.bed -b f2_bedfiles_merged/MCF-7/MCF-7_FOS_percentile_T_ExtendMerge.merge.bed -wa -c > f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T_ExtendMerge.merge.FOS.overlapped.bed
bedtools intersect -a f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T_ExtendMerge.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_exons.bed -wa -c > f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T_ExtendMerge.merge.hg38_exons.overlapped.bed
bedtools intersect -a f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T_ExtendMerge.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_introns.bed -wa -c > f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T_ExtendMerge.merge.hg38_introns.overlapped.bed
bedtools intersect -a f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T_ExtendMerge.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_4k_promoter_geneID.bed -wa -c > f3_coBinding_merge/MCF-7_top_TFBSCP/percentile_T_ExtendMerge.merge.hg38_4k_promoter_geneID.overlapped.bed

cat f2_bedfiles_merged/MCF-7/MCF-7_ERG_percentile_T_ExtendMerge.merge.bed \
f2_bedfiles_merged/MCF-7/MCF-7_JUND_percentile_T_ExtendMerge.merge.bed \
f2_bedfiles_merged/MCF-7/MCF-7_ZNF143_percentile_T_ExtendMerge.merge.bed > f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T_ExtendMerge.cat.bed

bedtools sort -i f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T_ExtendMerge.cat.bed > f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T_ExtendMerge.sort.bed
bedtools merge -i f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T_ExtendMerge.sort.bed > f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.bed

bedtools intersect -a f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.bed -b f2_bedfiles_merged/MCF-7/MCF-7_ERG_percentile_T_ExtendMerge.merge.bed -wa -c > f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.ERG.overlapped.bed
bedtools intersect -a f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.bed -b f2_bedfiles_merged/MCF-7/MCF-7_JUND_percentile_T_ExtendMerge.merge.bed -wa -c > f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.JUND.overlapped.bed
bedtools intersect -a f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.bed -b f2_bedfiles_merged/MCF-7/MCF-7_ZNF143_percentile_T_ExtendMerge.merge.bed -wa -c > f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.ZNF143.overlapped.bed
bedtools intersect -a f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_exons.bed -wa -c > f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.hg38_exons.overlapped.bed
bedtools intersect -a f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_introns.bed -wa -c > f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.hg38_introns.overlapped.bed
bedtools intersect -a f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_4k_promoter_geneID.bed -wa -c > f3_coBinding_merge/MCF-7_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.hg38_4k_promoter_geneID.overlapped.bed

cat f2_bedfiles_merged/HCT-116/HCT-116_ELF1_percentile_T_ExtendMerge.merge.bed \
f2_bedfiles_merged/HCT-116/HCT-116_JUND_percentile_T_ExtendMerge.merge.bed \
f2_bedfiles_merged/HCT-116/HCT-116_SRF_percentile_T_ExtendMerge.merge.bed > f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T_ExtendMerge.cat.bed

bedtools sort -i f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T_ExtendMerge.cat.bed > f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T_ExtendMerge.sort.bed
bedtools merge -i f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T_ExtendMerge.sort.bed > f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T_ExtendMerge.merge.bed

bedtools intersect -a f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T_ExtendMerge.merge.bed -b f2_bedfiles_merged/HCT-116/HCT-116_ELF1_percentile_T_ExtendMerge.merge.bed -wa -c > f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T_ExtendMerge.merge.ELF1.overlapped.bed
bedtools intersect -a f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T_ExtendMerge.merge.bed -b f2_bedfiles_merged/HCT-116/HCT-116_JUND_percentile_T_ExtendMerge.merge.bed -wa -c > f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T_ExtendMerge.merge.JUND.overlapped.bed
bedtools intersect -a f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T_ExtendMerge.merge.bed -b f2_bedfiles_merged/HCT-116/HCT-116_SRF_percentile_T_ExtendMerge.merge.bed -wa -c > f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T_ExtendMerge.merge.SRF.overlapped.bed
bedtools intersect -a f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T_ExtendMerge.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_exons.bed -wa -c > f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T_ExtendMerge.merge.hg38_exons.overlapped.bed
bedtools intersect -a f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T_ExtendMerge.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_introns.bed -wa -c > f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T_ExtendMerge.merge.hg38_introns.overlapped.bed
bedtools intersect -a f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T_ExtendMerge.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_4k_promoter_geneID.bed -wa -c > f3_coBinding_merge/HCT-116_top_TFBSCP/percentile_T_ExtendMerge.merge.hg38_4k_promoter_geneID.overlapped.bed

cat f2_bedfiles_merged/HCT-116/HCT-116_JUND_percentile_T_ExtendMerge.merge.bed \
f2_bedfiles_merged/HCT-116/HCT-116_CEBPB_percentile_T_ExtendMerge.merge.bed \
f2_bedfiles_merged/HCT-116/HCT-116_YY1_percentile_T_ExtendMerge.merge.bed > f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T_ExtendMerge.cat.bed

bedtools sort -i f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T_ExtendMerge.cat.bed > f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T_ExtendMerge.sort.bed
bedtools merge -i f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T_ExtendMerge.sort.bed > f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.bed

bedtools intersect -a f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.bed -b f2_bedfiles_merged/HCT-116/HCT-116_JUND_percentile_T_ExtendMerge.merge.bed -wa -c > f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.JUND.overlapped.bed
bedtools intersect -a f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.bed -b f2_bedfiles_merged/HCT-116/HCT-116_CEBPB_percentile_T_ExtendMerge.merge.bed -wa -c > f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.CEBPB.overlapped.bed
bedtools intersect -a f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.bed -b f2_bedfiles_merged/HCT-116/HCT-116_YY1_percentile_T_ExtendMerge.merge.bed -wa -c > f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.YY1.overlapped.bed
bedtools intersect -a f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_exons.bed -wa -c > f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.hg38_exons.overlapped.bed
bedtools intersect -a f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_introns.bed -wa -c > f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.hg38_introns.overlapped.bed
bedtools intersect -a f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.bed -b /nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/hg38_4k_promoter_geneID.bed -wa -c > f3_coBinding_merge/HCT-116_top_zscored_TFBSCP/percentile_T_ExtendMerge.merge.hg38_4k_promoter_geneID.overlapped.bed

immune Tregs Adenocarcinoma
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Tregs/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Tregs/Adenocarcinoma_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Tregs/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Tregs/Adenocarcinoma_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Tregs/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Tregs/Adenocarcinoma_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Tregs/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Tregs/Adenocarcinoma_SE_NOT_overlapped.bed

immune Tregs Normal
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Tregs/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Tregs/Normal_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Tregs/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Tregs/Normal_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Tregs/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Tregs/Normal_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Tregs/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Tregs/Normal_SE_NOT_overlapped.bed

immune Tregs Polyp
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Tregs/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Tregs/Polyp_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Tregs/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Tregs/Polyp_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Tregs/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Tregs/Polyp_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Tregs/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Tregs/Polyp_SE_NOT_overlapped.bed

immune Tregs Unaffected
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Tregs/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Tregs/Unaffected_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Tregs/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Tregs/Unaffected_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Tregs/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Tregs/Unaffected_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Tregs/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Tregs/Unaffected_SE_NOT_overlapped.bed

immune Macrophages Adenocarcinoma
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Macrophages/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Macrophages/Adenocarcinoma_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Macrophages/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Macrophages/Adenocarcinoma_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Macrophages/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Macrophages/Adenocarcinoma_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Macrophages/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Macrophages/Adenocarcinoma_SE_NOT_overlapped.bed

immune Macrophages Normal
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Macrophages/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Macrophages/Normal_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Macrophages/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Macrophages/Normal_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Macrophages/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Macrophages/Normal_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Macrophages/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Macrophages/Normal_SE_NOT_overlapped.bed

immune Macrophages Polyp
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Macrophages/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Macrophages/Polyp_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Macrophages/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Macrophages/Polyp_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Macrophages/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Macrophages/Polyp_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Macrophages/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Macrophages/Polyp_SE_NOT_overlapped.bed

immune Macrophages Unaffected
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Macrophages/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Macrophages/Unaffected_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Macrophages/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Macrophages/Unaffected_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Macrophages/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Macrophages/Unaffected_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Macrophages/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Macrophages/Unaffected_SE_NOT_overlapped.bed

immune GC Adenocarcinoma
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/GC/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/GC/Adenocarcinoma_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/GC/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/GC/Adenocarcinoma_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/GC/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/GC/Adenocarcinoma_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/GC/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/GC/Adenocarcinoma_SE_NOT_overlapped.bed

immune GC Normal
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/GC/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/GC/Normal_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/GC/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/GC/Normal_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/GC/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/GC/Normal_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/GC/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/GC/Normal_SE_NOT_overlapped.bed

immune GC Polyp
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/GC/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/GC/Polyp_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/GC/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/GC/Polyp_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/GC/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/GC/Polyp_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/GC/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/GC/Polyp_SE_NOT_overlapped.bed

immune GC Unaffected
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/GC/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/GC/Unaffected_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/GC/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/GC/Unaffected_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/GC/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/GC/Unaffected_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/GC/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/GC/Unaffected_SE_NOT_overlapped.bed

immune NK Adenocarcinoma
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/NK/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/NK/Adenocarcinoma_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/NK/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/NK/Adenocarcinoma_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/NK/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/NK/Adenocarcinoma_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/NK/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/NK/Adenocarcinoma_SE_NOT_overlapped.bed

immune NK Normal
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/NK/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/NK/Normal_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/NK/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/NK/Normal_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/NK/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/NK/Normal_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/NK/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/NK/Normal_SE_NOT_overlapped.bed

immune NK Polyp
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/NK/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/NK/Polyp_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/NK/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/NK/Polyp_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/NK/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/NK/Polyp_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/NK/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/NK/Polyp_SE_NOT_overlapped.bed

immune NK Unaffected
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/NK/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/NK/Unaffected_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/NK/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/NK/Unaffected_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/NK/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/NK/Unaffected_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/NK/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/NK/Unaffected_SE_NOT_overlapped.bed

immune DC Adenocarcinoma
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/DC/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/DC/Adenocarcinoma_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/DC/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/DC/Adenocarcinoma_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/DC/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/DC/Adenocarcinoma_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/DC/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/DC/Adenocarcinoma_SE_NOT_overlapped.bed

immune DC Normal
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/DC/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/DC/Normal_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/DC/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/DC/Normal_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/DC/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/DC/Normal_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/DC/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/DC/Normal_SE_NOT_overlapped.bed

immune DC Polyp
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/DC/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/DC/Polyp_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/DC/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/DC/Polyp_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/DC/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/DC/Polyp_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/DC/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/DC/Polyp_SE_NOT_overlapped.bed

immune DC Unaffected
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/DC/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/DC/Unaffected_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/DC/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/DC/Unaffected_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/DC/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/DC/Unaffected_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/DC/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/DC/Unaffected_SE_NOT_overlapped.bed

immune Mast2 Adenocarcinoma
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Mast2/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Mast2/Adenocarcinoma_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Mast2/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Mast2/Adenocarcinoma_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Mast2/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Mast2/Adenocarcinoma_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Mast2/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Mast2/Adenocarcinoma_SE_NOT_overlapped.bed

immune Mast2 Normal
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Mast2/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Mast2/Normal_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Mast2/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Mast2/Normal_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Mast2/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Mast2/Normal_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Mast2/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Mast2/Normal_SE_NOT_overlapped.bed

immune Mast2 Polyp
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Mast2/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Mast2/Polyp_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Mast2/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Mast2/Polyp_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Mast2/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Mast2/Polyp_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Mast2/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Mast2/Polyp_SE_NOT_overlapped.bed

immune Mast2 Unaffected
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Mast2/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Mast2/Unaffected_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Mast2/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Mast2/Unaffected_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Mast2/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Mast2/Unaffected_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Mast2/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Mast2/Unaffected_SE_NOT_overlapped.bed

immune CD4+ Tfh PD1+ Adenocarcinoma
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+_Tfh_PD1+/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/CD4+_Tfh_PD1+/Adenocarcinoma_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+_Tfh_PD1+/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/CD4+_Tfh_PD1+/Adenocarcinoma_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+_Tfh_PD1+/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/CD4+_Tfh_PD1+/Adenocarcinoma_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+_Tfh_PD1+/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/CD4+_Tfh_PD1+/Adenocarcinoma_SE_NOT_overlapped.bed

immune CD4+ Tfh PD1+ Normal
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+_Tfh_PD1+/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/CD4+_Tfh_PD1+/Normal_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+_Tfh_PD1+/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/CD4+_Tfh_PD1+/Normal_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+_Tfh_PD1+/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/CD4+_Tfh_PD1+/Normal_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+_Tfh_PD1+/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/CD4+_Tfh_PD1+/Normal_SE_NOT_overlapped.bed

immune CD4+ Tfh PD1+ Polyp
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+_Tfh_PD1+/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/CD4+_Tfh_PD1+/Polyp_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+_Tfh_PD1+/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/CD4+_Tfh_PD1+/Polyp_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+_Tfh_PD1+/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/CD4+_Tfh_PD1+/Polyp_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+_Tfh_PD1+/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/CD4+_Tfh_PD1+/Polyp_SE_NOT_overlapped.bed

immune CD4+ Tfh PD1+ Unaffected
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+_Tfh_PD1+/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/CD4+_Tfh_PD1+/Unaffected_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+_Tfh_PD1+/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/CD4+_Tfh_PD1+/Unaffected_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+_Tfh_PD1+/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/CD4+_Tfh_PD1+/Unaffected_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+_Tfh_PD1+/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/CD4+_Tfh_PD1+/Unaffected_SE_NOT_overlapped.bed

immune Plasma Adenocarcinoma
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Plasma/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Plasma/Adenocarcinoma_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Plasma/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Plasma/Adenocarcinoma_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Plasma/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Plasma/Adenocarcinoma_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Plasma/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Plasma/Adenocarcinoma_SE_NOT_overlapped.bed

immune Plasma Normal
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Plasma/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Plasma/Normal_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Plasma/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Plasma/Normal_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Plasma/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Plasma/Normal_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Plasma/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Plasma/Normal_SE_NOT_overlapped.bed

immune Plasma Polyp
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Plasma/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Plasma/Polyp_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Plasma/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Plasma/Polyp_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Plasma/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Plasma/Polyp_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Plasma/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Plasma/Polyp_SE_NOT_overlapped.bed

immune Plasma Unaffected
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Plasma/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Plasma/Unaffected_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Plasma/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Plasma/Unaffected_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Plasma/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Plasma/Unaffected_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Plasma/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Plasma/Unaffected_SE_NOT_overlapped.bed

immune CD4+ Activated Adenocarcinoma
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+_Activated/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/CD4+_Activated/Adenocarcinoma_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+_Activated/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/CD4+_Activated/Adenocarcinoma_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+_Activated/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/CD4+_Activated/Adenocarcinoma_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+_Activated/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/CD4+_Activated/Adenocarcinoma_SE_NOT_overlapped.bed

immune CD4+ Activated Normal
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+_Activated/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/CD4+_Activated/Normal_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+_Activated/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/CD4+_Activated/Normal_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+_Activated/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/CD4+_Activated/Normal_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+_Activated/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/CD4+_Activated/Normal_SE_NOT_overlapped.bed

immune CD4+ Activated Polyp
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+_Activated/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/CD4+_Activated/Polyp_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+_Activated/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/CD4+_Activated/Polyp_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+_Activated/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/CD4+_Activated/Polyp_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+_Activated/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/CD4+_Activated/Polyp_SE_NOT_overlapped.bed

immune CD4+ Activated Unaffected
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+_Activated/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/CD4+_Activated/Unaffected_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+_Activated/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/CD4+_Activated/Unaffected_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+_Activated/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/CD4+_Activated/Unaffected_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+_Activated/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/CD4+_Activated/Unaffected_SE_NOT_overlapped.bed

immune ILCs Adenocarcinoma
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/ILCs/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/ILCs/Adenocarcinoma_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/ILCs/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/ILCs/Adenocarcinoma_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/ILCs/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/ILCs/Adenocarcinoma_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/ILCs/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/ILCs/Adenocarcinoma_SE_NOT_overlapped.bed

immune ILCs Normal
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/ILCs/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/ILCs/Normal_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/ILCs/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/ILCs/Normal_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/ILCs/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/ILCs/Normal_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/ILCs/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/ILCs/Normal_SE_NOT_overlapped.bed

immune ILCs Polyp
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/ILCs/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/ILCs/Polyp_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/ILCs/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/ILCs/Polyp_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/ILCs/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/ILCs/Polyp_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/ILCs/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/ILCs/Polyp_SE_NOT_overlapped.bed

immune ILCs Unaffected
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/ILCs/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/ILCs/Unaffected_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/ILCs/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/ILCs/Unaffected_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/ILCs/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/ILCs/Unaffected_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/ILCs/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/ILCs/Unaffected_SE_NOT_overlapped.bed

immune CD8+ Adenocarcinoma
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD8+/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/CD8+/Adenocarcinoma_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD8+/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/CD8+/Adenocarcinoma_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD8+/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/CD8+/Adenocarcinoma_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD8+/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/CD8+/Adenocarcinoma_SE_NOT_overlapped.bed

immune CD8+ Normal
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD8+/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/CD8+/Normal_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD8+/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/CD8+/Normal_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD8+/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/CD8+/Normal_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD8+/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/CD8+/Normal_SE_NOT_overlapped.bed

immune CD8+ Polyp
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD8+/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/CD8+/Polyp_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD8+/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/CD8+/Polyp_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD8+/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/CD8+/Polyp_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD8+/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/CD8+/Polyp_SE_NOT_overlapped.bed

immune CD8+ Unaffected
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD8+/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/CD8+/Unaffected_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD8+/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/CD8+/Unaffected_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD8+/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/CD8+/Unaffected_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD8+/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/CD8+/Unaffected_SE_NOT_overlapped.bed

immune Exhausted T cells Adenocarcinoma
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Exhausted_T_cells/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Exhausted_T_cells/Adenocarcinoma_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Exhausted_T_cells/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Exhausted_T_cells/Adenocarcinoma_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Exhausted_T_cells/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Exhausted_T_cells/Adenocarcinoma_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Exhausted_T_cells/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Exhausted_T_cells/Adenocarcinoma_SE_NOT_overlapped.bed

immune Exhausted T cells Normal
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Exhausted_T_cells/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Exhausted_T_cells/Normal_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Exhausted_T_cells/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Exhausted_T_cells/Normal_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Exhausted_T_cells/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Exhausted_T_cells/Normal_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Exhausted_T_cells/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Exhausted_T_cells/Normal_SE_NOT_overlapped.bed

immune Exhausted T cells Polyp
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Exhausted_T_cells/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Exhausted_T_cells/Polyp_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Exhausted_T_cells/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Exhausted_T_cells/Polyp_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Exhausted_T_cells/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Exhausted_T_cells/Polyp_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Exhausted_T_cells/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Exhausted_T_cells/Polyp_SE_NOT_overlapped.bed

immune Exhausted T cells Unaffected
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Exhausted_T_cells/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Exhausted_T_cells/Unaffected_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Exhausted_T_cells/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Exhausted_T_cells/Unaffected_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Exhausted_T_cells/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Exhausted_T_cells/Unaffected_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Exhausted_T_cells/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Exhausted_T_cells/Unaffected_SE_NOT_overlapped.bed

immune Naive T Adenocarcinoma
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Naive_T/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Naive_T/Adenocarcinoma_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Naive_T/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Naive_T/Adenocarcinoma_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Naive_T/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Naive_T/Adenocarcinoma_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Naive_T/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Naive_T/Adenocarcinoma_SE_NOT_overlapped.bed

immune Naive T Polyp
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Naive_T/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Naive_T/Polyp_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Naive_T/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Naive_T/Polyp_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Naive_T/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Naive_T/Polyp_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Naive_T/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Naive_T/Polyp_SE_NOT_overlapped.bed

immune Naive T Unaffected
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Naive_T/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Naive_T/Unaffected_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Naive_T/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Naive_T/Unaffected_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Naive_T/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Naive_T/Unaffected_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Naive_T/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Naive_T/Unaffected_SE_NOT_overlapped.bed

immune Naive B Adenocarcinoma
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Naive_B/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Naive_B/Adenocarcinoma_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Naive_B/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Naive_B/Adenocarcinoma_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Naive_B/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Naive_B/Adenocarcinoma_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Naive_B/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Naive_B/Adenocarcinoma_SE_NOT_overlapped.bed

immune Naive B Normal
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Naive_B/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Naive_B/Normal_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Naive_B/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Naive_B/Normal_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Naive_B/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Naive_B/Normal_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Naive_B/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Naive_B/Normal_SE_NOT_overlapped.bed

immune Naive B Polyp
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Naive_B/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Naive_B/Polyp_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Naive_B/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Naive_B/Polyp_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Naive_B/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Naive_B/Polyp_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Naive_B/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Naive_B/Polyp_SE_NOT_overlapped.bed

immune Naive B Unaffected
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Naive_B/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Naive_B/Unaffected_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Naive_B/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Naive_B/Unaffected_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Naive_B/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Naive_B/Unaffected_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Naive_B/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Naive_B/Unaffected_SE_NOT_overlapped.bed

immune CD4+ Adenocarcinoma
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/CD4+/Adenocarcinoma_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/CD4+/Adenocarcinoma_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/CD4+/Adenocarcinoma_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/CD4+/Adenocarcinoma_SE_NOT_overlapped.bed

immune CD4+ Normal
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/CD4+/Normal_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/CD4+/Normal_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/CD4+/Normal_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/CD4+/Normal_SE_NOT_overlapped.bed

immune CD4+ Polyp
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/CD4+/Polyp_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/CD4+/Polyp_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/CD4+/Polyp_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/CD4+/Polyp_SE_NOT_overlapped.bed

immune CD4+ Unaffected
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/CD4+/Unaffected_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/CD4+/Unaffected_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/CD4+/Unaffected_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/CD4+/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/CD4+/Unaffected_SE_NOT_overlapped.bed

immune Mast Adenocarcinoma
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Mast/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Mast/Adenocarcinoma_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Mast/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Mast/Adenocarcinoma_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Mast/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Mast/Adenocarcinoma_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Mast/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Mast/Adenocarcinoma_SE_NOT_overlapped.bed

immune Mast Normal
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Mast/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Mast/Normal_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Mast/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Mast/Normal_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Mast/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Mast/Normal_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Mast/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Mast/Normal_SE_NOT_overlapped.bed

immune Mast Polyp
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Mast/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Mast/Polyp_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Mast/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Mast/Polyp_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Mast/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Mast/Polyp_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Mast/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Mast/Polyp_SE_NOT_overlapped.bed

immune Mast Unaffected
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Mast/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Mast/Unaffected_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Mast/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Mast/Unaffected_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Mast/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Mast/Unaffected_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Mast/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Mast/Unaffected_SE_NOT_overlapped.bed

immune Memory B Adenocarcinoma
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Memory_B/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Memory_B/Adenocarcinoma_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Memory_B/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Memory_B/Adenocarcinoma_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Memory_B/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Memory_B/Adenocarcinoma_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Memory_B/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Memory_B/Adenocarcinoma_SE_NOT_overlapped.bed

immune Memory B Normal
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Memory_B/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Memory_B/Normal_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Memory_B/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Memory_B/Normal_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Memory_B/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Memory_B/Normal_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Memory_B/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Memory_B/Normal_SE_NOT_overlapped.bed

immune Memory B Polyp
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Memory_B/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Memory_B/Polyp_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Memory_B/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Memory_B/Polyp_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Memory_B/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Memory_B/Polyp_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Memory_B/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Memory_B/Polyp_SE_NOT_overlapped.bed

immune Memory B Unaffected
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Memory_B/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Memory_B/Unaffected_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Memory_B/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Memory_B/Unaffected_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Memory_B/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/immune/Memory_B/Unaffected_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/immune/Memory_B/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/immune/Memory_B/Unaffected_SE_NOT_overlapped.bed

stromal Cancer Associated Fibroblasts Adenocarcinoma
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Cancer_Associated_Fibroblasts/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Cancer_Associated_Fibroblasts/Adenocarcinoma_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Cancer_Associated_Fibroblasts/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Cancer_Associated_Fibroblasts/Adenocarcinoma_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Cancer_Associated_Fibroblasts/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Cancer_Associated_Fibroblasts/Adenocarcinoma_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Cancer_Associated_Fibroblasts/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Cancer_Associated_Fibroblasts/Adenocarcinoma_SE_NOT_overlapped.bed

stromal Cancer Associated Fibroblasts Normal
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Cancer_Associated_Fibroblasts/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Cancer_Associated_Fibroblasts/Normal_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Cancer_Associated_Fibroblasts/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Cancer_Associated_Fibroblasts/Normal_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Cancer_Associated_Fibroblasts/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Cancer_Associated_Fibroblasts/Normal_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Cancer_Associated_Fibroblasts/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Cancer_Associated_Fibroblasts/Normal_SE_NOT_overlapped.bed

stromal Cancer Associated Fibroblasts Polyp
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Cancer_Associated_Fibroblasts/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Cancer_Associated_Fibroblasts/Polyp_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Cancer_Associated_Fibroblasts/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Cancer_Associated_Fibroblasts/Polyp_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Cancer_Associated_Fibroblasts/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Cancer_Associated_Fibroblasts/Polyp_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Cancer_Associated_Fibroblasts/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Cancer_Associated_Fibroblasts/Polyp_SE_NOT_overlapped.bed

stromal Cancer Associated Fibroblasts Unaffected
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Cancer_Associated_Fibroblasts/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Cancer_Associated_Fibroblasts/Unaffected_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Cancer_Associated_Fibroblasts/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Cancer_Associated_Fibroblasts/Unaffected_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Cancer_Associated_Fibroblasts/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Cancer_Associated_Fibroblasts/Unaffected_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Cancer_Associated_Fibroblasts/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Cancer_Associated_Fibroblasts/Unaffected_SE_NOT_overlapped.bed

stromal preCAF Adenocarcinoma
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/preCAF/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/preCAF/Adenocarcinoma_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/preCAF/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/preCAF/Adenocarcinoma_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/preCAF/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/preCAF/Adenocarcinoma_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/preCAF/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/preCAF/Adenocarcinoma_SE_NOT_overlapped.bed

stromal preCAF Normal
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/preCAF/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/preCAF/Normal_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/preCAF/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/preCAF/Normal_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/preCAF/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/preCAF/Normal_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/preCAF/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/preCAF/Normal_SE_NOT_overlapped.bed

stromal preCAF Polyp
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/preCAF/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/preCAF/Polyp_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/preCAF/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/preCAF/Polyp_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/preCAF/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/preCAF/Polyp_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/preCAF/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/preCAF/Polyp_SE_NOT_overlapped.bed

stromal preCAF Unaffected
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/preCAF/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/preCAF/Unaffected_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/preCAF/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/preCAF/Unaffected_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/preCAF/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/preCAF/Unaffected_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/preCAF/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/preCAF/Unaffected_SE_NOT_overlapped.bed

stromal Crypt Fibroblasts 3 Adenocarcinoma
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_3/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_3/Adenocarcinoma_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_3/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_3/Adenocarcinoma_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_3/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_3/Adenocarcinoma_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_3/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_3/Adenocarcinoma_SE_NOT_overlapped.bed

stromal Crypt Fibroblasts 3 Normal
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_3/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_3/Normal_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_3/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_3/Normal_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_3/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_3/Normal_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_3/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_3/Normal_SE_NOT_overlapped.bed

stromal Crypt Fibroblasts 3 Polyp
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_3/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_3/Polyp_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_3/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_3/Polyp_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_3/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_3/Polyp_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_3/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_3/Polyp_SE_NOT_overlapped.bed

stromal Crypt Fibroblasts 3 Unaffected
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_3/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_3/Unaffected_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_3/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_3/Unaffected_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_3/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_3/Unaffected_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_3/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_3/Unaffected_SE_NOT_overlapped.bed

stromal Unknown Normal
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Unknown/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Unknown/Normal_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Unknown/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Unknown/Normal_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Unknown/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Unknown/Normal_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Unknown/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Unknown/Normal_SE_NOT_overlapped.bed

stromal Crypt Fibroblasts 2 Adenocarcinoma
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_2/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_2/Adenocarcinoma_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_2/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_2/Adenocarcinoma_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_2/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_2/Adenocarcinoma_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_2/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_2/Adenocarcinoma_SE_NOT_overlapped.bed

stromal Crypt Fibroblasts 2 Normal
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_2/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_2/Normal_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_2/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_2/Normal_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_2/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_2/Normal_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_2/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_2/Normal_SE_NOT_overlapped.bed

stromal Crypt Fibroblasts 2 Polyp
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_2/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_2/Polyp_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_2/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_2/Polyp_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_2/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_2/Polyp_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_2/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_2/Polyp_SE_NOT_overlapped.bed

stromal Crypt Fibroblasts 2 Unaffected
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_2/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_2/Unaffected_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_2/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_2/Unaffected_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_2/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_2/Unaffected_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_2/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_2/Unaffected_SE_NOT_overlapped.bed

stromal Glia Adenocarcinoma
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Glia/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Glia/Adenocarcinoma_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Glia/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Glia/Adenocarcinoma_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Glia/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Glia/Adenocarcinoma_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Glia/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Glia/Adenocarcinoma_SE_NOT_overlapped.bed

stromal Glia Normal
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Glia/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Glia/Normal_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Glia/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Glia/Normal_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Glia/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Glia/Normal_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Glia/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Glia/Normal_SE_NOT_overlapped.bed

stromal Glia Polyp
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Glia/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Glia/Polyp_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Glia/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Glia/Polyp_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Glia/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Glia/Polyp_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Glia/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Glia/Polyp_SE_NOT_overlapped.bed

stromal Glia Unaffected
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Glia/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Glia/Unaffected_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Glia/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Glia/Unaffected_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Glia/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Glia/Unaffected_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Glia/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Glia/Unaffected_SE_NOT_overlapped.bed

stromal Endothelial Adenocarcinoma
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Endothelial/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Endothelial/Adenocarcinoma_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Endothelial/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Endothelial/Adenocarcinoma_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Endothelial/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Endothelial/Adenocarcinoma_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Endothelial/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Endothelial/Adenocarcinoma_SE_NOT_overlapped.bed

stromal Endothelial Normal
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Endothelial/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Endothelial/Normal_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Endothelial/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Endothelial/Normal_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Endothelial/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Endothelial/Normal_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Endothelial/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Endothelial/Normal_SE_NOT_overlapped.bed

stromal Endothelial Polyp
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Endothelial/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Endothelial/Polyp_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Endothelial/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Endothelial/Polyp_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Endothelial/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Endothelial/Polyp_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Endothelial/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Endothelial/Polyp_SE_NOT_overlapped.bed

stromal Endothelial Unaffected
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Endothelial/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Endothelial/Unaffected_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Endothelial/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Endothelial/Unaffected_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Endothelial/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Endothelial/Unaffected_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Endothelial/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Endothelial/Unaffected_SE_NOT_overlapped.bed

stromal Crypt Fibroblasts 1 Adenocarcinoma
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_1/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_1/Adenocarcinoma_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_1/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_1/Adenocarcinoma_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_1/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_1/Adenocarcinoma_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_1/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_1/Adenocarcinoma_SE_NOT_overlapped.bed

stromal Crypt Fibroblasts 1 Normal
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_1/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_1/Normal_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_1/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_1/Normal_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_1/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_1/Normal_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_1/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_1/Normal_SE_NOT_overlapped.bed

stromal Crypt Fibroblasts 1 Polyp
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_1/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_1/Polyp_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_1/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_1/Polyp_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_1/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_1/Polyp_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_1/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_1/Polyp_SE_NOT_overlapped.bed

stromal Crypt Fibroblasts 1 Unaffected
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_1/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_1/Unaffected_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_1/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_1/Unaffected_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_1/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_1/Unaffected_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_1/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_1/Unaffected_SE_NOT_overlapped.bed

stromal Myofibroblasts/Smooth Muscle 1 Adenocarcinoma
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Myofibroblasts_Smooth_Muscle_1/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Myofibroblasts_Smooth_Muscle_1/Adenocarcinoma_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Myofibroblasts_Smooth_Muscle_1/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Myofibroblasts_Smooth_Muscle_1/Adenocarcinoma_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Myofibroblasts_Smooth_Muscle_1/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Myofibroblasts_Smooth_Muscle_1/Adenocarcinoma_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Myofibroblasts_Smooth_Muscle_1/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Myofibroblasts_Smooth_Muscle_1/Adenocarcinoma_SE_NOT_overlapped.bed

stromal Myofibroblasts/Smooth Muscle 1 Normal
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Myofibroblasts_Smooth_Muscle_1/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Myofibroblasts_Smooth_Muscle_1/Normal_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Myofibroblasts_Smooth_Muscle_1/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Myofibroblasts_Smooth_Muscle_1/Normal_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Myofibroblasts_Smooth_Muscle_1/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Myofibroblasts_Smooth_Muscle_1/Normal_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Myofibroblasts_Smooth_Muscle_1/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Myofibroblasts_Smooth_Muscle_1/Normal_SE_NOT_overlapped.bed

stromal Myofibroblasts/Smooth Muscle 1 Polyp
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Myofibroblasts_Smooth_Muscle_1/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Myofibroblasts_Smooth_Muscle_1/Polyp_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Myofibroblasts_Smooth_Muscle_1/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Myofibroblasts_Smooth_Muscle_1/Polyp_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Myofibroblasts_Smooth_Muscle_1/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Myofibroblasts_Smooth_Muscle_1/Polyp_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Myofibroblasts_Smooth_Muscle_1/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Myofibroblasts_Smooth_Muscle_1/Polyp_SE_NOT_overlapped.bed

stromal Myofibroblasts/Smooth Muscle 1 Unaffected
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Myofibroblasts_Smooth_Muscle_1/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Myofibroblasts_Smooth_Muscle_1/Unaffected_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Myofibroblasts_Smooth_Muscle_1/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Myofibroblasts_Smooth_Muscle_1/Unaffected_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Myofibroblasts_Smooth_Muscle_1/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Myofibroblasts_Smooth_Muscle_1/Unaffected_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Myofibroblasts_Smooth_Muscle_1/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Myofibroblasts_Smooth_Muscle_1/Unaffected_SE_NOT_overlapped.bed

stromal Lymphatic Endothelial Cells Adenocarcinoma
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Lymphatic_Endothelial_Cells/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Lymphatic_Endothelial_Cells/Adenocarcinoma_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Lymphatic_Endothelial_Cells/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Lymphatic_Endothelial_Cells/Adenocarcinoma_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Lymphatic_Endothelial_Cells/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Lymphatic_Endothelial_Cells/Adenocarcinoma_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Lymphatic_Endothelial_Cells/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Lymphatic_Endothelial_Cells/Adenocarcinoma_SE_NOT_overlapped.bed

stromal Lymphatic Endothelial Cells Normal
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Lymphatic_Endothelial_Cells/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Lymphatic_Endothelial_Cells/Normal_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Lymphatic_Endothelial_Cells/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Lymphatic_Endothelial_Cells/Normal_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Lymphatic_Endothelial_Cells/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Lymphatic_Endothelial_Cells/Normal_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Lymphatic_Endothelial_Cells/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Lymphatic_Endothelial_Cells/Normal_SE_NOT_overlapped.bed

stromal Lymphatic Endothelial Cells Polyp
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Lymphatic_Endothelial_Cells/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Lymphatic_Endothelial_Cells/Polyp_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Lymphatic_Endothelial_Cells/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Lymphatic_Endothelial_Cells/Polyp_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Lymphatic_Endothelial_Cells/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Lymphatic_Endothelial_Cells/Polyp_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Lymphatic_Endothelial_Cells/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Lymphatic_Endothelial_Cells/Polyp_SE_NOT_overlapped.bed

stromal Lymphatic Endothelial Cells Unaffected
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Lymphatic_Endothelial_Cells/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Lymphatic_Endothelial_Cells/Unaffected_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Lymphatic_Endothelial_Cells/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Lymphatic_Endothelial_Cells/Unaffected_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Lymphatic_Endothelial_Cells/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Lymphatic_Endothelial_Cells/Unaffected_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Lymphatic_Endothelial_Cells/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Lymphatic_Endothelial_Cells/Unaffected_SE_NOT_overlapped.bed

stromal Myofibroblasts/Smooth Muscle 3 Normal
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Myofibroblasts_Smooth_Muscle_3/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Myofibroblasts_Smooth_Muscle_3/Normal_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Myofibroblasts_Smooth_Muscle_3/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Myofibroblasts_Smooth_Muscle_3/Normal_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Myofibroblasts_Smooth_Muscle_3/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Myofibroblasts_Smooth_Muscle_3/Normal_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Myofibroblasts_Smooth_Muscle_3/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Myofibroblasts_Smooth_Muscle_3/Normal_SE_NOT_overlapped.bed

stromal Myofibroblasts/Smooth Muscle 3 Polyp
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Myofibroblasts_Smooth_Muscle_3/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Myofibroblasts_Smooth_Muscle_3/Polyp_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Myofibroblasts_Smooth_Muscle_3/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Myofibroblasts_Smooth_Muscle_3/Polyp_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Myofibroblasts_Smooth_Muscle_3/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Myofibroblasts_Smooth_Muscle_3/Polyp_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Myofibroblasts_Smooth_Muscle_3/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Myofibroblasts_Smooth_Muscle_3/Polyp_SE_NOT_overlapped.bed

stromal Myofibroblasts/Smooth Muscle 3 Unaffected
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Myofibroblasts_Smooth_Muscle_3/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Myofibroblasts_Smooth_Muscle_3/Unaffected_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Myofibroblasts_Smooth_Muscle_3/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Myofibroblasts_Smooth_Muscle_3/Unaffected_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Myofibroblasts_Smooth_Muscle_3/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Myofibroblasts_Smooth_Muscle_3/Unaffected_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Myofibroblasts_Smooth_Muscle_3/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Myofibroblasts_Smooth_Muscle_3/Unaffected_SE_NOT_overlapped.bed

stromal Crypt Fibroblasts 4 Adenocarcinoma
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_4/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_4/Adenocarcinoma_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_4/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_4/Adenocarcinoma_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_4/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_4/Adenocarcinoma_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_4/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_4/Adenocarcinoma_SE_NOT_overlapped.bed

stromal Crypt Fibroblasts 4 Normal
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_4/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_4/Normal_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_4/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_4/Normal_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_4/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_4/Normal_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_4/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_4/Normal_SE_NOT_overlapped.bed

stromal Crypt Fibroblasts 4 Polyp
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_4/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_4/Polyp_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_4/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_4/Polyp_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_4/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_4/Polyp_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_4/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_4/Polyp_SE_NOT_overlapped.bed

stromal Crypt Fibroblasts 4 Unaffected
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_4/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_4/Unaffected_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_4/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_4/Unaffected_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_4/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_4/Unaffected_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Crypt_Fibroblasts_4/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Crypt_Fibroblasts_4/Unaffected_SE_NOT_overlapped.bed

stromal Pericytes Adenocarcinoma
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Pericytes/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Pericytes/Adenocarcinoma_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Pericytes/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Pericytes/Adenocarcinoma_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Pericytes/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Pericytes/Adenocarcinoma_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Pericytes/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Pericytes/Adenocarcinoma_SE_NOT_overlapped.bed

stromal Pericytes Normal
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Pericytes/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Pericytes/Normal_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Pericytes/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Pericytes/Normal_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Pericytes/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Pericytes/Normal_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Pericytes/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Pericytes/Normal_SE_NOT_overlapped.bed

stromal Pericytes Polyp
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Pericytes/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Pericytes/Polyp_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Pericytes/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Pericytes/Polyp_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Pericytes/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Pericytes/Polyp_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Pericytes/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Pericytes/Polyp_SE_NOT_overlapped.bed

stromal Pericytes Unaffected
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Pericytes/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Pericytes/Unaffected_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Pericytes/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Pericytes/Unaffected_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Pericytes/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Pericytes/Unaffected_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Pericytes/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Pericytes/Unaffected_SE_NOT_overlapped.bed

stromal Myofibroblasts/Smooth Muscle 2 Adenocarcinoma
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Myofibroblasts_Smooth_Muscle_2/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Myofibroblasts_Smooth_Muscle_2/Adenocarcinoma_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Myofibroblasts_Smooth_Muscle_2/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Myofibroblasts_Smooth_Muscle_2/Adenocarcinoma_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Myofibroblasts_Smooth_Muscle_2/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Myofibroblasts_Smooth_Muscle_2/Adenocarcinoma_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Myofibroblasts_Smooth_Muscle_2/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Myofibroblasts_Smooth_Muscle_2/Adenocarcinoma_SE_NOT_overlapped.bed

stromal Myofibroblasts/Smooth Muscle 2 Normal
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Myofibroblasts_Smooth_Muscle_2/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Myofibroblasts_Smooth_Muscle_2/Normal_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Myofibroblasts_Smooth_Muscle_2/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Myofibroblasts_Smooth_Muscle_2/Normal_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Myofibroblasts_Smooth_Muscle_2/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Myofibroblasts_Smooth_Muscle_2/Normal_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Myofibroblasts_Smooth_Muscle_2/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Myofibroblasts_Smooth_Muscle_2/Normal_SE_NOT_overlapped.bed

stromal Myofibroblasts/Smooth Muscle 2 Polyp
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Myofibroblasts_Smooth_Muscle_2/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Myofibroblasts_Smooth_Muscle_2/Polyp_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Myofibroblasts_Smooth_Muscle_2/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Myofibroblasts_Smooth_Muscle_2/Polyp_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Myofibroblasts_Smooth_Muscle_2/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Myofibroblasts_Smooth_Muscle_2/Polyp_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Myofibroblasts_Smooth_Muscle_2/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Myofibroblasts_Smooth_Muscle_2/Polyp_SE_NOT_overlapped.bed

stromal Myofibroblasts/Smooth Muscle 2 Unaffected
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Myofibroblasts_Smooth_Muscle_2/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Myofibroblasts_Smooth_Muscle_2/Unaffected_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Myofibroblasts_Smooth_Muscle_2/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Myofibroblasts_Smooth_Muscle_2/Unaffected_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Myofibroblasts_Smooth_Muscle_2/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Myofibroblasts_Smooth_Muscle_2/Unaffected_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Myofibroblasts_Smooth_Muscle_2/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Myofibroblasts_Smooth_Muscle_2/Unaffected_SE_NOT_overlapped.bed

stromal Villus Fibroblasts WNT5B+ Adenocarcinoma
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Villus_Fibroblasts_WNT5B+/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Villus_Fibroblasts_WNT5B+/Adenocarcinoma_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Villus_Fibroblasts_WNT5B+/Adenocarcinoma.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Villus_Fibroblasts_WNT5B+/Adenocarcinoma_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Villus_Fibroblasts_WNT5B+/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Villus_Fibroblasts_WNT5B+/Adenocarcinoma_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Villus_Fibroblasts_WNT5B+/Adenocarcinoma.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Villus_Fibroblasts_WNT5B+/Adenocarcinoma_SE_NOT_overlapped.bed

stromal Villus Fibroblasts WNT5B+ Normal
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Villus_Fibroblasts_WNT5B+/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Villus_Fibroblasts_WNT5B+/Normal_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Villus_Fibroblasts_WNT5B+/Normal.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Villus_Fibroblasts_WNT5B+/Normal_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Villus_Fibroblasts_WNT5B+/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Villus_Fibroblasts_WNT5B+/Normal_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Villus_Fibroblasts_WNT5B+/Normal.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Villus_Fibroblasts_WNT5B+/Normal_SE_NOT_overlapped.bed

stromal Villus Fibroblasts WNT5B+ Polyp
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Villus_Fibroblasts_WNT5B+/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Villus_Fibroblasts_WNT5B+/Polyp_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Villus_Fibroblasts_WNT5B+/Polyp.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Villus_Fibroblasts_WNT5B+/Polyp_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Villus_Fibroblasts_WNT5B+/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Villus_Fibroblasts_WNT5B+/Polyp_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Villus_Fibroblasts_WNT5B+/Polyp.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Villus_Fibroblasts_WNT5B+/Polyp_SE_NOT_overlapped.bed

stromal Villus Fibroblasts WNT5B+ Unaffected
bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Villus_Fibroblasts_WNT5B+/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Villus_Fibroblasts_WNT5B+/Unaffected_co_binding_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Villus_Fibroblasts_WNT5B+/Unaffected.bed \
-b f1_continuum_changes_at_coBinding/data_HCT-116_SRF_JUND_CEBPB_coBinding.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Villus_Fibroblasts_WNT5B+/Unaffected_co_binding_NOT_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Villus_Fibroblasts_WNT5B+/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -u > f1_continuum_changes_at_coBinding/stromal/Villus_Fibroblasts_WNT5B+/Unaffected_SE_overlapped.bed

bedtools intersect -a /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/f3_public_data/data/greenLeaf_NG2022//processed_cell_Regions_MERGE_by_GrossPathology/stromal/Villus_Fibroblasts_WNT5B+/Unaffected.bed \
-b /nv/vol190/zanglab/zw5j/since2019_projects/phase_separation_FEpiTR//f9_TF_condensates_V3/data/SE_hg38/HCT-116.bed \
-wa -v > f1_continuum_changes_at_coBinding/stromal/Villus_Fibroblasts_WNT5B+/Unaffected_SE_NOT_overlapped.bed


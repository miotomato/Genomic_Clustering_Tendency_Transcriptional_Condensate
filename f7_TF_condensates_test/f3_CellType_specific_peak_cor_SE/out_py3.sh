cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/72784_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/47270_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/45306_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/39730_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/908_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/57094_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/53321_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H3K4me3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H3K4me3.cat.bed > f3_union_peaks_across_cellTypes/H3K4me3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H3K4me3.sort.bed > f3_union_peaks_across_cellTypes/H3K4me3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/72784_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H3K4me3_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/47270_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H3K4me3_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/45306_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H3K4me3_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/39730_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H3K4me3_merge_intersect_A549.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/908_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H3K4me3_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/57094_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H3K4me3_merge_intersect_HCT-116.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/53321_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H3K4me3_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/72784_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_H3K4me3_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/47270_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_H3K4me3_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/45306_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_H3K4me3_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/39730_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_H3K4me3_A549.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/908_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_H3K4me3_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/57094_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_H3K4me3_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/53321_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_H3K4me3_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/33859_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45985_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46256_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46061_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46199_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46201_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/73296_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/POLR2A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/POLR2A.cat.bed > f3_union_peaks_across_cellTypes/POLR2A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/POLR2A.sort.bed > f3_union_peaks_across_cellTypes/POLR2A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/33859_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/POLR2A_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45985_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/POLR2A_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46256_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/POLR2A_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46061_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/POLR2A_merge_intersect_A549.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46199_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/POLR2A_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46201_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/POLR2A_merge_intersect_HCT-116.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/73296_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/POLR2A_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/33859_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_POLR2A_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45985_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_POLR2A_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46256_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_POLR2A_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46061_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_POLR2A_A549.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46199_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_POLR2A_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46201_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_POLR2A_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/73296_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_POLR2A_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/35400_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/33141_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/45299_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/39026_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/6531_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/62110_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/62017_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H3K27ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H3K27ac.cat.bed > f3_union_peaks_across_cellTypes/H3K27ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H3K27ac.sort.bed > f3_union_peaks_across_cellTypes/H3K27ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/35400_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H3K27ac_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/33141_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H3K27ac_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/45299_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H3K27ac_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/39026_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H3K27ac_merge_intersect_A549.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/6531_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H3K27ac_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/62110_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H3K27ac_merge_intersect_HCT-116.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/62017_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H3K27ac_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/35400_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_H3K27ac_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/33141_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_H3K27ac_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/45299_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_H3K27ac_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/39026_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_H3K27ac_A549.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/6531_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_H3K27ac_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/62110_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_H3K27ac_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/62017_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_H3K27ac_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45899_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/83736_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/1380_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45558_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/50987_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/42149_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/8087_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/CTCF.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/CTCF.cat.bed > f3_union_peaks_across_cellTypes/CTCF.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/CTCF.sort.bed > f3_union_peaks_across_cellTypes/CTCF.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45899_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/CTCF_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/83736_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/CTCF_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/1380_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/CTCF_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45558_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/CTCF_merge_intersect_A549.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/50987_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/CTCF_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/42149_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/CTCF_merge_intersect_HCT-116.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/8087_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/CTCF_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45899_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_CTCF_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/83736_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_CTCF_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/1380_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_CTCF_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45558_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_CTCF_A549.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/50987_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_CTCF_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/42149_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_CTCF_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/8087_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_CTCF_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/45382_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/75112_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/45307_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/39041_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/897_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/42910_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/45288_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H3K4me2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H3K4me2.cat.bed > f3_union_peaks_across_cellTypes/H3K4me2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H3K4me2.sort.bed > f3_union_peaks_across_cellTypes/H3K4me2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/45382_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H3K4me2_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/75112_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H3K4me2_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/45307_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H3K4me2_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/39041_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H3K4me2_merge_intersect_A549.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/897_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H3K4me2_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/42910_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H3K4me2_merge_intersect_HCT-116.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/45288_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H3K4me2_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/45382_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_H3K4me2_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/75112_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_H3K4me2_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/45307_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_H3K4me2_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/39041_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_H3K4me2_A549.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/897_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_H3K4me2_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/42910_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_H3K4me2_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/45288_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_H3K4me2_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45963_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/33155_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46225_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46056_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46165_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/42907_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45728_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/EP300.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/EP300.cat.bed > f3_union_peaks_across_cellTypes/EP300.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/EP300.sort.bed > f3_union_peaks_across_cellTypes/EP300.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45963_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/EP300_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/33155_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/EP300_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46225_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/EP300_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46056_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/EP300_merge_intersect_A549.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46165_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/EP300_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/42907_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/EP300_merge_intersect_HCT-116.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45728_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/EP300_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45963_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_EP300_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/33155_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_EP300_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46225_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_EP300_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46056_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_EP300_A549.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46165_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_EP300_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/42907_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_EP300_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45728_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_EP300_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45894_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/38193_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45827_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45560_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45691_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/70806_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45745_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MYC.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MYC.cat.bed > f3_union_peaks_across_cellTypes/MYC.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MYC.sort.bed > f3_union_peaks_across_cellTypes/MYC.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45894_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MYC_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/38193_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MYC_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45827_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MYC_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45560_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MYC_merge_intersect_A549.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45691_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MYC_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/70806_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MYC_merge_intersect_HCT-116.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45745_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MYC_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45894_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_MYC_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/38193_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_MYC_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45827_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_MYC_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45560_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_MYC_A549.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45691_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_MYC_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/70806_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_MYC_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45745_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_MYC_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/45406_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/86292_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/45305_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/39729_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/45227_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/87280_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/67816_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H3K9ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H3K9ac.cat.bed > f3_union_peaks_across_cellTypes/H3K9ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H3K9ac.sort.bed > f3_union_peaks_across_cellTypes/H3K9ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/45406_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H3K9ac_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/86292_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H3K9ac_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/45305_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H3K9ac_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/39729_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H3K9ac_merge_intersect_A549.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/45227_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H3K9ac_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/87280_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H3K9ac_merge_intersect_HCT-116.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/67816_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H3K9ac_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/45406_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_H3K9ac_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/86292_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_H3K9ac_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/45305_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_H3K9ac_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/39729_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_H3K9ac_A549.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/45227_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_H3K9ac_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/87280_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_H3K9ac_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_hm/67816_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_H3K9ac_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45923_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46316_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45814_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45559_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45681_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46206_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45762_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/CEBPB.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/CEBPB.cat.bed > f3_union_peaks_across_cellTypes/CEBPB.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/CEBPB.sort.bed > f3_union_peaks_across_cellTypes/CEBPB.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/CEBPB.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45923_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/CEBPB_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/CEBPB.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46316_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/CEBPB_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/CEBPB.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45814_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/CEBPB_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/CEBPB.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45559_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/CEBPB_merge_intersect_A549.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/CEBPB.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45681_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/CEBPB_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/CEBPB.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46206_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/CEBPB_merge_intersect_HCT-116.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/CEBPB.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45762_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/CEBPB_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45923_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_CEBPB_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46316_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_CEBPB_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45814_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_CEBPB_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45559_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_CEBPB_A549.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45681_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_CEBPB_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46206_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_CEBPB_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45762_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_CEBPB_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46275_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/2387_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46252_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46084_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46164_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46207_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45757_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/RAD21.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/RAD21.cat.bed > f3_union_peaks_across_cellTypes/RAD21.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/RAD21.sort.bed > f3_union_peaks_across_cellTypes/RAD21.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/RAD21.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46275_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/RAD21_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/RAD21.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/2387_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/RAD21_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/RAD21.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46252_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/RAD21_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/RAD21.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46084_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/RAD21_merge_intersect_A549.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/RAD21.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46164_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/RAD21_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/RAD21.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46207_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/RAD21_merge_intersect_HCT-116.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/RAD21.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45757_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/RAD21_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46275_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_RAD21_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/2387_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_RAD21_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46252_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_RAD21_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46084_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_RAD21_A549.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46164_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_RAD21_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46207_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_RAD21_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45757_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_RAD21_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46286_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46313_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46257_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46069_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46172_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46203_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46222_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/REST.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/REST.cat.bed > f3_union_peaks_across_cellTypes/REST.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/REST.sort.bed > f3_union_peaks_across_cellTypes/REST.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/REST.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46286_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/REST_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/REST.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46313_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/REST_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/REST.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46257_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/REST_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/REST.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46069_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/REST_merge_intersect_A549.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/REST.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46172_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/REST_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/REST.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46203_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/REST_merge_intersect_HCT-116.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/REST.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46222_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/REST_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46286_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_REST_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46313_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_REST_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46257_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_REST_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46069_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_REST_A549.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46172_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_REST_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46203_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_REST_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46222_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_REST_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46292_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46324_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46244_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46070_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46180_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46216_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/38119_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MAX.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MAX.cat.bed > f3_union_peaks_across_cellTypes/MAX.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MAX.sort.bed > f3_union_peaks_across_cellTypes/MAX.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MAX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46292_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MAX_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MAX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46324_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MAX_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MAX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46244_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MAX_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MAX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46070_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MAX_merge_intersect_A549.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MAX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46180_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MAX_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MAX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46216_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MAX_merge_intersect_HCT-116.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MAX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/38119_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MAX_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46292_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_MAX_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46324_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_MAX_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46244_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_MAX_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46070_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_MAX_A549.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46180_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_MAX_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46216_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_MAX_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/38119_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_MAX_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46268_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46314_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/56045_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46080_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46185_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46219_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GABPA.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GABPA.cat.bed > f3_union_peaks_across_cellTypes/GABPA.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GABPA.sort.bed > f3_union_peaks_across_cellTypes/GABPA.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GABPA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46268_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GABPA_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GABPA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46314_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GABPA_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GABPA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/56045_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GABPA_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GABPA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46080_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GABPA_merge_intersect_A549.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GABPA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46185_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GABPA_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GABPA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46219_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GABPA_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46268_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_GABPA_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46314_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_GABPA_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/56045_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_GABPA_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46080_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_GABPA_A549.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46185_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_GABPA_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46219_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_GABPA_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46282_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46320_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46258_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46077_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46196_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46221_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/TAF1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/TAF1.cat.bed > f3_union_peaks_across_cellTypes/TAF1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/TAF1.sort.bed > f3_union_peaks_across_cellTypes/TAF1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/TAF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46282_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/TAF1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/TAF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46320_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/TAF1_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/TAF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46258_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/TAF1_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/TAF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46077_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/TAF1_merge_intersect_A549.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/TAF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46196_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/TAF1_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/TAF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46221_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/TAF1_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46282_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_TAF1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46320_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_TAF1_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46258_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_TAF1_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46077_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_TAF1_A549.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46196_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_TAF1_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46221_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_TAF1_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46303_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46226_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46065_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46170_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46204_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/74510_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/YY1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/YY1.cat.bed > f3_union_peaks_across_cellTypes/YY1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/YY1.sort.bed > f3_union_peaks_across_cellTypes/YY1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/YY1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46303_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/YY1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/YY1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46226_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/YY1_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/YY1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46065_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/YY1_merge_intersect_A549.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/YY1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46170_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/YY1_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/YY1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46204_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/YY1_merge_intersect_HCT-116.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/YY1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/74510_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/YY1_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46303_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_YY1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46226_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_YY1_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46065_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_YY1_A549.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46170_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_YY1_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46204_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_YY1_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/74510_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_YY1_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46284_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46315_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46238_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46075_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46189_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46210_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/TEAD4.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/TEAD4.cat.bed > f3_union_peaks_across_cellTypes/TEAD4.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/TEAD4.sort.bed > f3_union_peaks_across_cellTypes/TEAD4.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/TEAD4.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46284_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/TEAD4_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/TEAD4.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46315_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/TEAD4_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/TEAD4.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46238_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/TEAD4_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/TEAD4.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46075_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/TEAD4_merge_intersect_A549.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/TEAD4.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46189_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/TEAD4_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/TEAD4.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46210_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/TEAD4_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46284_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_TEAD4_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46315_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_TEAD4_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46238_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_TEAD4_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46075_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_TEAD4_A549.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46189_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_TEAD4_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46210_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_TEAD4_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45880_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46308_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45828_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45696_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46213_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45743_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/JUND.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/JUND.cat.bed > f3_union_peaks_across_cellTypes/JUND.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/JUND.sort.bed > f3_union_peaks_across_cellTypes/JUND.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/JUND.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45880_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/JUND_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/JUND.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46308_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/JUND_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/JUND.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45828_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/JUND_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/JUND.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45696_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/JUND_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/JUND.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46213_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/JUND_merge_intersect_HCT-116.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/JUND.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45743_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/JUND_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45880_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_JUND_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46308_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_JUND_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45828_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_JUND_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45696_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_JUND_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/46213_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_JUND_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/45743_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_JUND_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/74667_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/49652_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/71877_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/71871_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MED1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MED1.cat.bed > f3_union_peaks_across_cellTypes/MED1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MED1.sort.bed > f3_union_peaks_across_cellTypes/MED1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MED1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/74667_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MED1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MED1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/49652_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MED1_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MED1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/71877_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MED1_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MED1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/71871_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MED1_merge_intersect_A549.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/74667_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_MED1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/49652_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_MED1_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/71877_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_MED1_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019_new/human_factor/71871_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/_MED1_A549.bed


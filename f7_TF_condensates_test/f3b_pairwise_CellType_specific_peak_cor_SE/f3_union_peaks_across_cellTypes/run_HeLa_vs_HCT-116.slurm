#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=8000
#SBATCH -t 4:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f3_union_peaks_across_cellTypes/slurm_HeLa_vs_HCT-116.out

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/73296_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46201_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/POLR2A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/POLR2A.cat.bed > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/POLR2A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/POLR2A.sort.bed > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/POLR2A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/73296_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/POLR2A_merge_intersect_HeLa.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46201_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/POLR2A_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/73296_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/_POLR2A_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46201_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/_POLR2A_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/53321_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57094_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/H3K4me3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/H3K4me3.cat.bed > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/H3K4me3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/H3K4me3.sort.bed > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/H3K4me3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/53321_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/H3K4me3_merge_intersect_HeLa.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57094_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/H3K4me3_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/53321_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/_H3K4me3_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57094_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/_H3K4me3_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62017_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62110_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/H3K27ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/H3K27ac.cat.bed > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/H3K27ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/H3K27ac.sort.bed > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/H3K27ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62017_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/H3K27ac_merge_intersect_HeLa.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62110_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/H3K27ac_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62017_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/_H3K27ac_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62110_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/_H3K27ac_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8087_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/42149_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/CTCF.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/CTCF.cat.bed > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/CTCF.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/CTCF.sort.bed > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/CTCF.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8087_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/CTCF_merge_intersect_HeLa.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/42149_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/CTCF_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8087_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/_CTCF_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/42149_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/_CTCF_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1252_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42910_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/H3K4me2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/H3K4me2.cat.bed > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/H3K4me2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/H3K4me2.sort.bed > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/H3K4me2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1252_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/H3K4me2_merge_intersect_HeLa.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42910_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/H3K4me2_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1252_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/_H3K4me2_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42910_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/_H3K4me2_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56234_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/42907_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/EP300.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/EP300.cat.bed > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/EP300.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/EP300.sort.bed > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/EP300.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56234_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/EP300_merge_intersect_HeLa.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/42907_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/EP300_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56234_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/_EP300_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/42907_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/_EP300_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/67816_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/87280_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/H3K9ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/H3K9ac.cat.bed > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/H3K9ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/H3K9ac.sort.bed > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/H3K9ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/67816_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/H3K9ac_merge_intersect_HeLa.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/87280_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/H3K9ac_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/67816_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/_H3K9ac_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/87280_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/_H3K9ac_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/38115_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/70806_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/MYC.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/MYC.cat.bed > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/MYC.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/MYC.sort.bed > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/MYC.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/38115_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/MYC_merge_intersect_HeLa.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/70806_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/MYC_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/38115_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/_MYC_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/70806_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/_MYC_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74510_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46204_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/YY1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/YY1.cat.bed > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/YY1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/YY1.sort.bed > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/YY1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/YY1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74510_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/YY1_merge_intersect_HeLa.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/YY1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46204_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/YY1_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74510_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/_YY1_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46204_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/_YY1_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/38119_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46216_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/MAX.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/MAX.cat.bed > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/MAX.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/MAX.sort.bed > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/MAX.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/MAX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/38119_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/MAX_merge_intersect_HeLa.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/MAX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46216_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/MAX_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/38119_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/_MAX_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46216_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/_MAX_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1231_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/83065_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/H2AZ.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/H2AZ.cat.bed > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/H2AZ.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/H2AZ.sort.bed > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/H2AZ.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/H2AZ.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1231_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/H2AZ_merge_intersect_HeLa.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/H2AZ.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/83065_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/H2AZ_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1231_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/_H2AZ_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/83065_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/_H2AZ_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/49002_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/55645_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/ELL2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/ELL2.cat.bed > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/ELL2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/ELL2.sort.bed > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/ELL2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/ELL2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/49002_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/ELL2_merge_intersect_HeLa.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/ELL2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/55645_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/ELL2_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/49002_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/_ELL2_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/55645_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/_ELL2_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/73297_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/57091_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/LARP7.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/LARP7.cat.bed > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/LARP7.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/LARP7.sort.bed > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/LARP7.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/LARP7.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/73297_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/LARP7_merge_intersect_HeLa.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/LARP7.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/57091_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/LARP7_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/73297_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/_LARP7_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/57091_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_HCT-116/_LARP7_HCT-116.bed


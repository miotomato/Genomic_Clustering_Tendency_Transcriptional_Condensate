#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=8000
#SBATCH -t 4:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f3_union_peaks_across_cellTypes/slurm_H1_vs_HCT-116.out

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46199_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46201_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/POLR2A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/POLR2A.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/POLR2A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/POLR2A.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/POLR2A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46199_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/POLR2A_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46201_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/POLR2A_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46199_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_POLR2A_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46201_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_POLR2A_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/908_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57094_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/H3K4me3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/H3K4me3.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/H3K4me3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/H3K4me3.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/H3K4me3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/908_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/H3K4me3_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57094_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/H3K4me3_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/908_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_H3K4me3_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57094_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_H3K4me3_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/6531_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62110_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/H3K27ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/H3K27ac.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/H3K27ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/H3K27ac.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/H3K27ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/6531_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/H3K27ac_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62110_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/H3K27ac_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/6531_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_H3K27ac_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62110_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_H3K27ac_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/50987_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/42149_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/CTCF.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/CTCF.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/CTCF.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/CTCF.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/CTCF.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/50987_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/CTCF_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/42149_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/CTCF_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/50987_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_CTCF_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/42149_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_CTCF_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/897_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42910_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/H3K4me2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/H3K4me2.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/H3K4me2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/H3K4me2.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/H3K4me2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/897_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/H3K4me2_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42910_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/H3K4me2_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/897_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_H3K4me2_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42910_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_H3K4me2_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46165_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/42907_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/EP300.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/EP300.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/EP300.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/EP300.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/EP300.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46165_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/EP300_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/42907_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/EP300_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46165_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_EP300_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/42907_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_EP300_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45227_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/87280_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/H3K9ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/H3K9ac.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/H3K9ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/H3K9ac.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/H3K9ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45227_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/H3K9ac_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/87280_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/H3K9ac_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45227_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_H3K9ac_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/87280_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_H3K9ac_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45691_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/70806_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/MYC.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/MYC.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/MYC.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/MYC.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/MYC.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45691_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/MYC_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/70806_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/MYC_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45691_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_MYC_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/70806_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_MYC_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/894_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57095_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/H3K4me1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/H3K4me1.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/H3K4me1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/H3K4me1.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/H3K4me1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/H3K4me1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/894_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/H3K4me1_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/H3K4me1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57095_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/H3K4me1_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/894_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_H3K4me1_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57095_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_H3K4me1_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46172_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46203_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/REST.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/REST.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/REST.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/REST.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/REST.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/REST.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46172_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/REST_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/REST.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46203_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/REST_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46172_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_REST_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46203_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_REST_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46170_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46204_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/YY1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/YY1.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/YY1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/YY1.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/YY1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/YY1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46170_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/YY1_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/YY1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46204_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/YY1_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46170_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_YY1_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46204_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_YY1_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46180_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46216_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/MAX.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/MAX.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/MAX.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/MAX.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/MAX.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/MAX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46180_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/MAX_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/MAX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46216_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/MAX_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46180_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_MAX_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46216_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_MAX_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46164_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46207_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/RAD21.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/RAD21.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/RAD21.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/RAD21.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/RAD21.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/RAD21.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46164_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/RAD21_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/RAD21.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46207_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/RAD21_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46164_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_RAD21_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46207_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_RAD21_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45681_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46206_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/CEBPB.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/CEBPB.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/CEBPB.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/CEBPB.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/CEBPB.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/CEBPB.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45681_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/CEBPB_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/CEBPB.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46206_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/CEBPB_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45681_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_CEBPB_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46206_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_CEBPB_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45696_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46213_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/JUND.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/JUND.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/JUND.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/JUND.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/JUND.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/JUND.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45696_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/JUND_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/JUND.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46213_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/JUND_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45696_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_JUND_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46213_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_JUND_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46189_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46210_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/TEAD4.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/TEAD4.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/TEAD4.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/TEAD4.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/TEAD4.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/TEAD4.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46189_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/TEAD4_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/TEAD4.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46210_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/TEAD4_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46189_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_TEAD4_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46210_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_TEAD4_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46167_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46211_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/SRF.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/SRF.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/SRF.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/SRF.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/SRF.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/SRF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46167_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/SRF_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/SRF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46211_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/SRF_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46167_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_SRF_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46211_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_SRF_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46171_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46205_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/USF1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/USF1.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/USF1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/USF1.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/USF1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/USF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46171_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/USF1_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/USF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46205_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/USF1_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46171_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_USF1_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46205_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_USF1_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45699_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46208_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/SIN3A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/SIN3A.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/SIN3A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/SIN3A.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/SIN3A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/SIN3A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45699_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/SIN3A_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/SIN3A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46208_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/SIN3A_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45699_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_SIN3A_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46208_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_SIN3A_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46190_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46215_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/SP1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/SP1.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/SP1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/SP1.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/SP1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/SP1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46190_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/SP1_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/SP1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46215_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/SP1_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46190_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_SP1_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46215_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_SP1_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46174_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/68464_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/ATF3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/ATF3.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/ATF3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/ATF3.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/ATF3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/ATF3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46174_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/ATF3_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/ATF3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/68464_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/ATF3_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46174_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_ATF3_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/68464_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_ATF3_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46169_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46214_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/EGR1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/EGR1.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/EGR1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/EGR1.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/EGR1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/EGR1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46169_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/EGR1_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/EGR1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46214_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/EGR1_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46169_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_EGR1_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46214_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_EGR1_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46175_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46212_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/FOSL1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/FOSL1.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/FOSL1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_HCT-116/FOSL1.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/FOSL1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/FOSL1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46175_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/FOSL1_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_HCT-116/FOSL1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46212_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_HCT-116/FOSL1_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46175_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_FOSL1_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46212_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_HCT-116/_FOSL1_HCT-116.bed


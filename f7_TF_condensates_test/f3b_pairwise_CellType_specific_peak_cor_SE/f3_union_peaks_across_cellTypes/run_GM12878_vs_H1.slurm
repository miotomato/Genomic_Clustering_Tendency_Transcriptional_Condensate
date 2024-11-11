#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=8000
#SBATCH -t 4:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f3_union_peaks_across_cellTypes/slurm_GM12878_vs_H1.out

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46137_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46199_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/POLR2A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/POLR2A.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/POLR2A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/POLR2A.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/POLR2A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46137_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/POLR2A_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46199_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/POLR2A_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46137_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_POLR2A_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46199_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_POLR2A_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42291_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/908_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K4me3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K4me3.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K4me3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K4me3.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K4me3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42291_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K4me3_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/908_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K4me3_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42291_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_H3K4me3_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/908_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_H3K4me3_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/48285_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/6531_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K27ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K27ac.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K27ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K27ac.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K27ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/48285_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K27ac_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/6531_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K27ac_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/48285_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_H3K27ac_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/6531_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_H3K27ac_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45631_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/50987_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/CTCF.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/CTCF.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/CTCF.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/CTCF.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/CTCF.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45631_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/CTCF_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/50987_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/CTCF_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45631_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_CTCF_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/50987_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_CTCF_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45214_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/897_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K4me2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K4me2.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K4me2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K4me2.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K4me2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45214_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K4me2_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/897_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K4me2_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45214_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_H3K4me2_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/897_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_H3K4me2_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/48272_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46165_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/EP300.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/EP300.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/EP300.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/EP300.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/EP300.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/48272_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/EP300_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46165_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/EP300_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/48272_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_EP300_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46165_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_EP300_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/36652_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45227_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K9ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K9ac.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K9ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K9ac.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K9ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/36652_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K9ac_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45227_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K9ac_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/36652_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_H3K9ac_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45227_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_H3K9ac_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/666_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45691_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/MYC.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/MYC.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/MYC.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/MYC.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/MYC.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/666_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/MYC_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45691_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/MYC_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/666_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_MYC_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45691_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_MYC_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/49297_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/894_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K4me1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K4me1.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K4me1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K4me1.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K4me1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K4me1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/49297_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K4me1_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K4me1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/894_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K4me1_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/49297_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_H3K4me1_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/894_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_H3K4me1_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46113_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46172_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/REST.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/REST.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/REST.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/REST.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/REST.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/REST.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46113_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/REST_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/REST.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46172_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/REST_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46113_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_REST_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46172_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_REST_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46136_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46185_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/GABPA.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/GABPA.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/GABPA.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/GABPA.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/GABPA.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/GABPA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46136_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/GABPA_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/GABPA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46185_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/GABPA_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46136_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_GABPA_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46185_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_GABPA_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/7060_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46170_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/YY1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/YY1.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/YY1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/YY1.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/YY1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/YY1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/7060_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/YY1_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/YY1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46170_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/YY1_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/7060_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_YY1_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46170_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_YY1_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/661_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46180_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/MAX.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/MAX.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/MAX.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/MAX.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/MAX.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/MAX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/661_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/MAX_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/MAX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46180_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/MAX_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/661_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_MAX_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46180_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_MAX_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46150_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46164_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/RAD21.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/RAD21.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/RAD21.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/RAD21.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/RAD21.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/RAD21.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46150_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/RAD21_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/RAD21.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46164_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/RAD21_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46150_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_RAD21_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46164_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_RAD21_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46108_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46196_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/TAF1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/TAF1.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/TAF1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/TAF1.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/TAF1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/TAF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46108_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/TAF1_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/TAF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46196_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/TAF1_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46108_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_TAF1_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46196_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_TAF1_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/660_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45696_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/JUND.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/JUND.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/JUND.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/JUND.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/JUND.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/JUND.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/660_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/JUND_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/JUND.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45696_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/JUND_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/660_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_JUND_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45696_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_JUND_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42282_b_sort_peaks.broadPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/8501_b_sort_peaks.broadPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K27me3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K27me3.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K27me3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K27me3.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K27me3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K27me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42282_b_sort_peaks.broadPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K27me3_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K27me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/8501_b_sort_peaks.broadPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K27me3_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42282_b_sort_peaks.broadPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_H3K27me3_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/8501_b_sort_peaks.broadPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_H3K27me3_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45629_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45688_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/ZNF143.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/ZNF143.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/ZNF143.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/ZNF143.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/ZNF143.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/ZNF143.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45629_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/ZNF143_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/ZNF143.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45688_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/ZNF143_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45629_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_ZNF143_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45688_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_ZNF143_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46120_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46166_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/TCF12.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/TCF12.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/TCF12.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/TCF12.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/TCF12.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/TCF12.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46120_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/TCF12_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/TCF12.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46166_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/TCF12_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46120_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_TCF12_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46166_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_TCF12_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45624_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45697_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/NRF1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/NRF1.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/NRF1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/NRF1.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/NRF1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/NRF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45624_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/NRF1_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/NRF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45697_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/NRF1_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45624_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_NRF1_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45697_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_NRF1_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46148_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46167_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/SRF.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/SRF.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/SRF.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/SRF.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/SRF.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/SRF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46148_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/SRF_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/SRF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46167_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/SRF_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46148_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_SRF_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46167_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_SRF_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46141_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46182_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/CREB1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/CREB1.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/CREB1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/CREB1.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/CREB1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/CREB1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46141_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/CREB1_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/CREB1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46182_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/CREB1_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46141_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_CREB1_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46182_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_CREB1_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46124_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46171_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/USF1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/USF1.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/USF1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/USF1.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/USF1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/USF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46124_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/USF1_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/USF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46171_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/USF1_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46124_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_USF1_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46171_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_USF1_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42284_b_sort_peaks.broadPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/6421_b_sort_peaks.broadPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K36me3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K36me3.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K36me3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K36me3.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K36me3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K36me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42284_b_sort_peaks.broadPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K36me3_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K36me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/6421_b_sort_peaks.broadPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/H3K36me3_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42284_b_sort_peaks.broadPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_H3K36me3_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/6421_b_sort_peaks.broadPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_H3K36me3_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45639_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45699_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/SIN3A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/SIN3A.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/SIN3A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/SIN3A.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/SIN3A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/SIN3A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45639_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/SIN3A_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/SIN3A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45699_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/SIN3A_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45639_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_SIN3A_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45699_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_SIN3A_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46135_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46190_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/SP1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/SP1.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/SP1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/SP1.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/SP1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/SP1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46135_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/SP1_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/SP1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46190_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/SP1_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46135_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_SP1_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46190_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_SP1_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45633_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45695_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/MXI1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/MXI1.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/MXI1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/MXI1.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/MXI1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/MXI1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45633_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/MXI1_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/MXI1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45695_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/MXI1_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45633_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_MXI1_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45695_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_MXI1_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45649_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45692_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/RFX5.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/RFX5.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/RFX5.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/RFX5.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/RFX5.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/RFX5.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45649_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/RFX5_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/RFX5.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45692_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/RFX5_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45649_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_RFX5_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45692_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_RFX5_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46133_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46174_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/ATF3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/ATF3.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/ATF3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/ATF3.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/ATF3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/ATF3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46133_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/ATF3_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/ATF3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46174_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/ATF3_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46133_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_ATF3_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46174_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_ATF3_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45621_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45687_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/MAFK.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/MAFK.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/MAFK.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/MAFK.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/MAFK.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/MAFK.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45621_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/MAFK_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/MAFK.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45687_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/MAFK_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45621_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_MAFK_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45687_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_MAFK_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45618_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45684_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/TBP.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/TBP.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/TBP.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/TBP.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/TBP.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/TBP.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45618_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/TBP_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/TBP.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45684_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/TBP_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45618_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_TBP_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45684_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_TBP_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46126_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46169_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/EGR1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/EGR1.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/EGR1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/EGR1.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/EGR1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/EGR1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46126_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/EGR1_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/EGR1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46169_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/EGR1_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46126_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_EGR1_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46169_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_EGR1_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45611_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45694_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/USF2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/USF2.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/USF2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/USF2.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/USF2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/USF2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45611_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/USF2_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/USF2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45694_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/USF2_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45611_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_USF2_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45694_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_USF2_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45638_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45701_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/CHD2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/CHD2.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/CHD2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/CHD2.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/CHD2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/CHD2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45638_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/CHD2_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/CHD2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45701_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/CHD2_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45638_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_CHD2_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45701_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_CHD2_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45654_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45698_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/BRCA1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/BRCA1.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/BRCA1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/BRCA1.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/BRCA1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/BRCA1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45654_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/BRCA1_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/BRCA1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45698_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/BRCA1_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45654_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_BRCA1_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45698_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_BRCA1_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46142_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46188_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/SIX5.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/SIX5.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/SIX5.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/SIX5.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/SIX5.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/SIX5.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46142_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/SIX5_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/SIX5.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46188_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/SIX5_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46142_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_SIX5_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46188_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_SIX5_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46112_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46193_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/ATF2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/ATF2.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/ATF2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/ATF2.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/ATF2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/ATF2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46112_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/ATF2_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/ATF2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46193_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/ATF2_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46112_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_ATF2_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46193_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_ATF2_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46151_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46178_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/BCL11A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/BCL11A.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/BCL11A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_H1/BCL11A.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_H1/BCL11A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/BCL11A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46151_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/BCL11A_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_H1/BCL11A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46178_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_H1/BCL11A_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46151_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_BCL11A_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46178_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_H1/_BCL11A_H1.bed


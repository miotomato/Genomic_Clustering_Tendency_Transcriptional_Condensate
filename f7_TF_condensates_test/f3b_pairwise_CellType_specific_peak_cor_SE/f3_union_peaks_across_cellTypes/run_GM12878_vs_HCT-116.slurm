#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=8000
#SBATCH -t 4:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f3_union_peaks_across_cellTypes/slurm_GM12878_vs_HCT-116.out

# cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46137_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46201_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/POLR2A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/POLR2A.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/POLR2A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/POLR2A.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/POLR2A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46137_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/POLR2A_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46201_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/POLR2A_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46137_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_POLR2A_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46201_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_POLR2A_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42291_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57094_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H3K4me3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H3K4me3.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H3K4me3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H3K4me3.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H3K4me3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42291_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H3K4me3_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57094_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H3K4me3_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42291_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_H3K4me3_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57094_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_H3K4me3_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/48285_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62110_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H3K27ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H3K27ac.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H3K27ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H3K27ac.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H3K27ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/48285_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H3K27ac_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62110_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H3K27ac_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/48285_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_H3K27ac_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62110_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_H3K27ac_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45631_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/42149_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/CTCF.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/CTCF.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/CTCF.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/CTCF.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/CTCF.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45631_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/CTCF_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/42149_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/CTCF_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45631_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_CTCF_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/42149_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_CTCF_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45214_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42910_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H3K4me2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H3K4me2.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H3K4me2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H3K4me2.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H3K4me2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45214_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H3K4me2_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42910_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H3K4me2_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45214_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_H3K4me2_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42910_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_H3K4me2_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/48272_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/42907_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/EP300.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/EP300.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/EP300.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/EP300.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/EP300.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/48272_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/EP300_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/42907_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/EP300_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/48272_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_EP300_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/42907_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_EP300_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/36652_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/87280_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H3K9ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H3K9ac.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H3K9ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H3K9ac.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H3K9ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/36652_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H3K9ac_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/87280_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H3K9ac_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/36652_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_H3K9ac_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/87280_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_H3K9ac_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/666_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/70806_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/MYC.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/MYC.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/MYC.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/MYC.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/MYC.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/666_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/MYC_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/70806_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/MYC_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/666_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_MYC_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/70806_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_MYC_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/49297_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57095_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H3K4me1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H3K4me1.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H3K4me1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H3K4me1.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H3K4me1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H3K4me1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/49297_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H3K4me1_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H3K4me1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57095_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H3K4me1_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/49297_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_H3K4me1_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57095_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_H3K4me1_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46113_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46203_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/REST.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/REST.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/REST.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/REST.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/REST.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/REST.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46113_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/REST_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/REST.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46203_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/REST_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46113_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_REST_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46203_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_REST_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/7060_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46204_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/YY1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/YY1.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/YY1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/YY1.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/YY1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/YY1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/7060_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/YY1_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/YY1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46204_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/YY1_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/7060_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_YY1_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46204_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_YY1_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/661_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46216_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/MAX.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/MAX.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/MAX.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/MAX.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/MAX.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/MAX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/661_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/MAX_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/MAX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46216_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/MAX_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/661_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_MAX_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46216_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_MAX_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46150_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46207_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/RAD21.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/RAD21.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/RAD21.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/RAD21.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/RAD21.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/RAD21.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46150_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/RAD21_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/RAD21.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46207_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/RAD21_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46150_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_RAD21_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46207_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_RAD21_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45221_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/83065_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H2AZ.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H2AZ.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H2AZ.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H2AZ.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H2AZ.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H2AZ.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45221_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H2AZ_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H2AZ.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/83065_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/H2AZ_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45221_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_H2AZ_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/83065_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_H2AZ_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/660_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46213_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/JUND.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/JUND.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/JUND.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/JUND.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/JUND.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/JUND.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/660_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/JUND_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/JUND.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46213_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/JUND_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/660_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_JUND_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46213_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_JUND_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46148_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46211_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/SRF.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/SRF.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/SRF.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/SRF.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/SRF.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/SRF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46148_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/SRF_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/SRF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46211_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/SRF_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46148_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_SRF_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46211_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_SRF_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46127_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46202_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/ELF1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/ELF1.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/ELF1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/ELF1.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/ELF1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/ELF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46127_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/ELF1_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/ELF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46202_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/ELF1_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46127_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_ELF1_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46202_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_ELF1_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46124_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46205_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/USF1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/USF1.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/USF1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/USF1.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/USF1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/USF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46124_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/USF1_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/USF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46205_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/USF1_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46124_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_USF1_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46205_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_USF1_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45639_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46208_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/SIN3A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/SIN3A.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/SIN3A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/SIN3A.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/SIN3A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/SIN3A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45639_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/SIN3A_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/SIN3A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46208_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/SIN3A_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45639_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_SIN3A_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46208_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_SIN3A_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46135_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46215_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/SP1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/SP1.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/SP1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/SP1.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/SP1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/SP1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46135_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/SP1_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/SP1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46215_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/SP1_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46135_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_SP1_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46215_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_SP1_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/48317_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46200_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/ZBTB33.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/ZBTB33.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/ZBTB33.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/ZBTB33.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/ZBTB33.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/ZBTB33.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/48317_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/ZBTB33_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/ZBTB33.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46200_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/ZBTB33_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/48317_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_ZBTB33_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46200_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_ZBTB33_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46133_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/68464_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/ATF3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/ATF3.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/ATF3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/ATF3.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/ATF3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/ATF3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46133_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/ATF3_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/ATF3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/68464_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/ATF3_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46133_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_ATF3_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/68464_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_ATF3_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46126_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46214_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/EGR1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/EGR1.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/EGR1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/EGR1.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/EGR1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/EGR1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46126_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/EGR1_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/EGR1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46214_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/EGR1_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46126_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_EGR1_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46214_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_EGR1_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46140_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/70804_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/TCF7L1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/TCF7L1.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/TCF7L1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/TCF7L1.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/TCF7L1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/TCF7L1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46140_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/TCF7L1_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/TCF7L1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/70804_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/TCF7L1_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46140_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_TCF7L1_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/70804_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HCT-116/_TCF7L1_HCT-116.bed


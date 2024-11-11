#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=8000
#SBATCH -t 4:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f3_union_peaks_across_cellTypes/slurm_GM12878_vs_HepG2.out

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46137_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46256_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/POLR2A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/POLR2A.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/POLR2A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/POLR2A.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/POLR2A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46137_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/POLR2A_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46256_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/POLR2A_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46137_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_POLR2A_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46256_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_POLR2A_HepG2.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42291_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45306_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/H3K4me3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/H3K4me3.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/H3K4me3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/H3K4me3.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/H3K4me3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42291_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/H3K4me3_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45306_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/H3K4me3_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42291_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_H3K4me3_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45306_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_H3K4me3_HepG2.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/48285_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45299_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/H3K27ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/H3K27ac.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/H3K27ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/H3K27ac.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/H3K27ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/48285_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/H3K27ac_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45299_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/H3K27ac_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/48285_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_H3K27ac_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45299_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_H3K27ac_HepG2.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45631_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1380_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/CTCF.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/CTCF.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/CTCF.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/CTCF.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/CTCF.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45631_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/CTCF_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1380_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/CTCF_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45631_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_CTCF_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1380_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_CTCF_HepG2.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45214_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45307_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/H3K4me2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/H3K4me2.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/H3K4me2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/H3K4me2.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/H3K4me2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45214_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/H3K4me2_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45307_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/H3K4me2_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45214_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_H3K4me2_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45307_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_H3K4me2_HepG2.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/48272_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46225_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/EP300.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/EP300.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/EP300.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/EP300.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/EP300.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/48272_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/EP300_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46225_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/EP300_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/48272_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_EP300_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46225_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_EP300_HepG2.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/36652_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45305_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/H3K9ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/H3K9ac.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/H3K9ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/H3K9ac.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/H3K9ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/36652_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/H3K9ac_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45305_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/H3K9ac_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/36652_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_H3K9ac_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45305_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_H3K9ac_HepG2.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/666_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45827_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/MYC.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/MYC.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/MYC.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/MYC.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/MYC.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/666_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/MYC_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45827_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/MYC_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/666_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_MYC_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45827_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_MYC_HepG2.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/49297_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/51255_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/H3K4me1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/H3K4me1.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/H3K4me1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/H3K4me1.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/H3K4me1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/H3K4me1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/49297_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/H3K4me1_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/H3K4me1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/51255_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/H3K4me1_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/49297_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_H3K4me1_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/51255_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_H3K4me1_HepG2.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46113_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46257_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/REST.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/REST.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/REST.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/REST.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/REST.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/REST.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46113_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/REST_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/REST.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46257_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/REST_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46113_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_REST_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46257_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_REST_HepG2.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46136_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56045_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/GABPA.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/GABPA.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/GABPA.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/GABPA.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/GABPA.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/GABPA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46136_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/GABPA_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/GABPA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56045_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/GABPA_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46136_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_GABPA_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56045_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_GABPA_HepG2.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/7060_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46226_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/YY1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/YY1.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/YY1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/YY1.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/YY1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/YY1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/7060_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/YY1_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/YY1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46226_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/YY1_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/7060_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_YY1_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46226_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_YY1_HepG2.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/661_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46244_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/MAX.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/MAX.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/MAX.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/MAX.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/MAX.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/MAX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/661_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/MAX_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/MAX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46244_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/MAX_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/661_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_MAX_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46244_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_MAX_HepG2.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46150_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46252_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/RAD21.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/RAD21.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/RAD21.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/RAD21.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/RAD21.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/RAD21.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46150_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/RAD21_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/RAD21.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46252_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/RAD21_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46150_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_RAD21_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46252_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_RAD21_HepG2.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46108_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46258_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/TAF1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/TAF1.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/TAF1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/TAF1.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/TAF1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/TAF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46108_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/TAF1_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/TAF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46258_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/TAF1_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46108_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_TAF1_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46258_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_TAF1_HepG2.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45640_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45816_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/SMC3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/SMC3.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/SMC3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/SMC3.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/SMC3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/SMC3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45640_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/SMC3_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/SMC3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45816_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/SMC3_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45640_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_SMC3_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45816_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_SMC3_HepG2.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/660_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45828_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/JUND.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/JUND.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/JUND.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/JUND.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/JUND.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/JUND.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/660_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/JUND_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/JUND.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45828_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/JUND_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/660_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_JUND_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45828_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_JUND_HepG2.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46120_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46228_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/TCF12.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/TCF12.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/TCF12.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/TCF12.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/TCF12.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/TCF12.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46120_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/TCF12_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/TCF12.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46228_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/TCF12_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46120_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_TCF12_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46228_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_TCF12_HepG2.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45624_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64667_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/NRF1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/NRF1.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/NRF1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/NRF1.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/NRF1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/NRF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45624_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/NRF1_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/NRF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64667_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/NRF1_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45624_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_NRF1_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64667_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_NRF1_HepG2.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46148_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46223_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/SRF.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/SRF.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/SRF.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/SRF.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/SRF.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/SRF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46148_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/SRF_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/SRF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46223_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/SRF_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46148_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_SRF_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46223_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_SRF_HepG2.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46127_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46254_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ELF1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ELF1.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ELF1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ELF1.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ELF1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ELF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46127_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ELF1_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ELF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46254_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ELF1_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46127_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_ELF1_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46254_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_ELF1_HepG2.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46141_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56047_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/CREB1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/CREB1.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/CREB1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/CREB1.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/CREB1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/CREB1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46141_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/CREB1_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/CREB1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56047_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/CREB1_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46141_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_CREB1_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56047_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_CREB1_HepG2.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46124_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46248_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/USF1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/USF1.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/USF1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/USF1.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/USF1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/USF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46124_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/USF1_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/USF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46248_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/USF1_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46124_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_USF1_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46248_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_USF1_HepG2.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45639_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46236_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/SIN3A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/SIN3A.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/SIN3A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/SIN3A.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/SIN3A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/SIN3A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45639_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/SIN3A_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/SIN3A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46236_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/SIN3A_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45639_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_SIN3A_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46236_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_SIN3A_HepG2.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46135_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46234_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/SP1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/SP1.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/SP1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/SP1.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/SP1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/SP1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46135_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/SP1_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/SP1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46234_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/SP1_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46135_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_SP1_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46234_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_SP1_HepG2.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45633_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45810_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/MXI1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/MXI1.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/MXI1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/MXI1.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/MXI1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/MXI1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45633_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/MXI1_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/MXI1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45810_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/MXI1_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45633_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_MXI1_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45810_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_MXI1_HepG2.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45649_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45823_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/RFX5.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/RFX5.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/RFX5.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/RFX5.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/RFX5.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/RFX5.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45649_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/RFX5_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/RFX5.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45823_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/RFX5_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45649_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_RFX5_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45823_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_RFX5_HepG2.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/48317_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46229_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ZBTB33.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ZBTB33.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ZBTB33.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ZBTB33.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ZBTB33.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ZBTB33.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/48317_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ZBTB33_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ZBTB33.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46229_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ZBTB33_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/48317_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_ZBTB33_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46229_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_ZBTB33_HepG2.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46133_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/7301_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ATF3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ATF3.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ATF3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ATF3.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ATF3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ATF3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46133_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ATF3_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ATF3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/7301_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ATF3_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46133_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_ATF3_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/7301_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_ATF3_HepG2.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45621_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45815_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/MAFK.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/MAFK.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/MAFK.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/MAFK.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/MAFK.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/MAFK.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45621_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/MAFK_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/MAFK.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45815_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/MAFK_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45621_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_MAFK_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45815_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_MAFK_HepG2.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45618_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45791_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/TBP.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/TBP.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/TBP.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/TBP.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/TBP.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/TBP.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45618_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/TBP_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/TBP.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45791_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/TBP_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45618_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_TBP_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45791_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_TBP_HepG2.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45661_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45788_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/MAZ.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/MAZ.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/MAZ.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/MAZ.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/MAZ.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/MAZ.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45661_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/MAZ_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/MAZ.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45788_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/MAZ_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45661_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_MAZ_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45788_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_MAZ_HepG2.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/42295_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1492_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/STAG1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/STAG1.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/STAG1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/STAG1.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/STAG1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/STAG1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/42295_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/STAG1_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/STAG1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1492_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/STAG1_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/42295_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_STAG1_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1492_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_STAG1_HepG2.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/48277_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1374_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/CEBPA.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/CEBPA.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/CEBPA.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/CEBPA.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/CEBPA.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/CEBPA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/48277_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/CEBPA_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/CEBPA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1374_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/CEBPA_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/48277_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_CEBPA_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1374_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_CEBPA_HepG2.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45632_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45807_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/BHLHE40.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/BHLHE40.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/BHLHE40.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/BHLHE40.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/BHLHE40.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/BHLHE40.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45632_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/BHLHE40_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/BHLHE40.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45807_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/BHLHE40_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45632_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_BHLHE40_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45807_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_BHLHE40_HepG2.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/669_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45800_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/NR2C2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/NR2C2.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/NR2C2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/NR2C2.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/NR2C2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/NR2C2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/669_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/NR2C2_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/NR2C2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45800_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/NR2C2_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/669_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_NR2C2_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45800_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_NR2C2_HepG2.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45611_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74655_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/USF2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/USF2.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/USF2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/USF2.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/USF2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/USF2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45611_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/USF2_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/USF2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74655_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/USF2_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45611_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_USF2_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74655_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_USF2_HepG2.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45638_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45802_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/CHD2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/CHD2.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/CHD2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/CHD2.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/CHD2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/CHD2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45638_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/CHD2_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/CHD2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45802_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/CHD2_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45638_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_CHD2_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45802_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_CHD2_HepG2.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45654_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45809_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/BRCA1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/BRCA1.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/BRCA1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/BRCA1.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/BRCA1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/BRCA1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45654_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/BRCA1_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/BRCA1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45809_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/BRCA1_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45654_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_BRCA1_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45809_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_BRCA1_HepG2.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45652_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45798_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/RCOR1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/RCOR1.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/RCOR1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/RCOR1.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/RCOR1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/RCOR1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45652_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/RCOR1_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/RCOR1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45798_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/RCOR1_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45652_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_RCOR1_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45798_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_RCOR1_HepG2.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62937_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63800_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ZBTB40.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ZBTB40.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ZBTB40.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ZBTB40.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ZBTB40.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ZBTB40.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62937_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ZBTB40_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ZBTB40.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63800_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ZBTB40_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62937_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_ZBTB40_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63800_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_ZBTB40_HepG2.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46111_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46240_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/NFIC.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/NFIC.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/NFIC.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/NFIC.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/NFIC.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/NFIC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46111_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/NFIC_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/NFIC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46240_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/NFIC_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46111_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_NFIC_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46240_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_NFIC_HepG2.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74652_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64657_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ZMYM3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ZMYM3.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ZMYM3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ZMYM3.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ZMYM3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ZMYM3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74652_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ZMYM3_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ZMYM3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64657_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ZMYM3_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74652_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_ZMYM3_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64657_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_ZMYM3_HepG2.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74683_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74659_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ZBTB2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ZBTB2.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ZBTB2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ZBTB2.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ZBTB2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ZBTB2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74683_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ZBTB2_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ZBTB2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74659_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ZBTB2_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74683_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_ZBTB2_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74659_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_ZBTB2_HepG2.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46139_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46245_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ZEB1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ZEB1.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ZEB1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ZEB1.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ZEB1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ZEB1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46139_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ZEB1_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ZEB1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46245_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/ZEB1_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46139_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_ZEB1_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46245_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_ZEB1_HepG2.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45620_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45793_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/SREBF1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/SREBF1.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/SREBF1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/SREBF1.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/SREBF1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/SREBF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45620_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/SREBF1_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/SREBF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45793_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/SREBF1_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45620_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_SREBF1_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45793_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HepG2/_SREBF1_HepG2.bed


#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=8000
#SBATCH -t 4:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f3_union_peaks_across_cellTypes/slurm_HepG2_vs_H1.out

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46256_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46199_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/POLR2A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/POLR2A.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/POLR2A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/POLR2A.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/POLR2A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46256_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/POLR2A_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46199_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/POLR2A_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46256_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_POLR2A_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46199_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_POLR2A_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45306_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/908_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K4me3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K4me3.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K4me3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K4me3.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K4me3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45306_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K4me3_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/908_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K4me3_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45306_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_H3K4me3_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/908_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_H3K4me3_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45299_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/6531_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K27ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K27ac.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K27ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K27ac.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K27ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45299_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K27ac_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/6531_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K27ac_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45299_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_H3K27ac_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/6531_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_H3K27ac_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1380_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/50987_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/CTCF.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/CTCF.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/CTCF.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/CTCF.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/CTCF.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1380_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/CTCF_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/50987_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/CTCF_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1380_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_CTCF_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/50987_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_CTCF_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45307_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/897_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K4me2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K4me2.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K4me2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K4me2.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K4me2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45307_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K4me2_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/897_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K4me2_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45307_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_H3K4me2_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/897_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_H3K4me2_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46225_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46165_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/EP300.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/EP300.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/EP300.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/EP300.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/EP300.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46225_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/EP300_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46165_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/EP300_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46225_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_EP300_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46165_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_EP300_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45305_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45227_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K9ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K9ac.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K9ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K9ac.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K9ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45305_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K9ac_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45227_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K9ac_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45305_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_H3K9ac_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45227_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_H3K9ac_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45827_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45691_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/MYC.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/MYC.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/MYC.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/MYC.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/MYC.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45827_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/MYC_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45691_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/MYC_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45827_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_MYC_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45691_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_MYC_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/51255_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/894_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K4me1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K4me1.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K4me1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K4me1.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K4me1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K4me1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/51255_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K4me1_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K4me1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/894_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K4me1_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/51255_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_H3K4me1_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/894_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_H3K4me1_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46257_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46172_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/REST.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/REST.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/REST.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/REST.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/REST.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/REST.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46257_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/REST_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/REST.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46172_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/REST_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46257_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_REST_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46172_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_REST_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56045_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46185_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/GABPA.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/GABPA.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/GABPA.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/GABPA.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/GABPA.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/GABPA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56045_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/GABPA_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/GABPA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46185_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/GABPA_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56045_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_GABPA_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46185_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_GABPA_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46226_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46170_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/YY1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/YY1.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/YY1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/YY1.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/YY1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/YY1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46226_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/YY1_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/YY1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46170_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/YY1_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46226_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_YY1_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46170_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_YY1_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46244_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46180_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/MAX.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/MAX.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/MAX.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/MAX.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/MAX.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/MAX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46244_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/MAX_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/MAX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46180_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/MAX_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46244_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_MAX_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46180_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_MAX_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46252_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46164_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/RAD21.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/RAD21.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/RAD21.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/RAD21.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/RAD21.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/RAD21.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46252_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/RAD21_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/RAD21.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46164_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/RAD21_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46252_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_RAD21_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46164_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_RAD21_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46258_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46196_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/TAF1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/TAF1.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/TAF1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/TAF1.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/TAF1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/TAF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46258_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/TAF1_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/TAF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46196_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/TAF1_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46258_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_TAF1_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46196_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_TAF1_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45814_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45681_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/CEBPB.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/CEBPB.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/CEBPB.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/CEBPB.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/CEBPB.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/CEBPB.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45814_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/CEBPB_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/CEBPB.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45681_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/CEBPB_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45814_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_CEBPB_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45681_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_CEBPB_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45828_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45696_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/JUND.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/JUND.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/JUND.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/JUND.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/JUND.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/JUND.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45828_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/JUND_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/JUND.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45696_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/JUND_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45828_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_JUND_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45696_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_JUND_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46238_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46189_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/TEAD4.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/TEAD4.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/TEAD4.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/TEAD4.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/TEAD4.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/TEAD4.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46238_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/TEAD4_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/TEAD4.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46189_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/TEAD4_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46238_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_TEAD4_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46189_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_TEAD4_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45300_b_sort_peaks.broadPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/48924_b_sort_peaks.broadPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K9me3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K9me3.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K9me3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K9me3.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K9me3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K9me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45300_b_sort_peaks.broadPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K9me3_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K9me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/48924_b_sort_peaks.broadPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/H3K9me3_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45300_b_sort_peaks.broadPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_H3K9me3_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/48924_b_sort_peaks.broadPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_H3K9me3_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46228_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46166_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/TCF12.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/TCF12.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/TCF12.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/TCF12.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/TCF12.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/TCF12.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46228_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/TCF12_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/TCF12.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46166_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/TCF12_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46228_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_TCF12_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46166_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_TCF12_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64667_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45697_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/NRF1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/NRF1.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/NRF1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/NRF1.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/NRF1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/NRF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64667_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/NRF1_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/NRF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45697_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/NRF1_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64667_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_NRF1_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45697_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_NRF1_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46223_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46167_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/SRF.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/SRF.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/SRF.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/SRF.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/SRF.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/SRF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46223_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/SRF_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/SRF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46167_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/SRF_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46223_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_SRF_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46167_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_SRF_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56047_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46182_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/CREB1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/CREB1.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/CREB1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/CREB1.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/CREB1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/CREB1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56047_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/CREB1_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/CREB1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46182_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/CREB1_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56047_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_CREB1_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46182_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_CREB1_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45830_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45685_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/JUN.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/JUN.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/JUN.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/JUN.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/JUN.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/JUN.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45830_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/JUN_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/JUN.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45685_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/JUN_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45830_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_JUN_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45685_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_JUN_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46248_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46171_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/USF1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/USF1.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/USF1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/USF1.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/USF1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/USF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46248_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/USF1_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/USF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46171_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/USF1_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46248_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_USF1_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46171_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_USF1_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46236_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45699_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/SIN3A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/SIN3A.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/SIN3A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/SIN3A.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/SIN3A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/SIN3A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46236_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/SIN3A_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/SIN3A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45699_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/SIN3A_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46236_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_SIN3A_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45699_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_SIN3A_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46234_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46190_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/SP1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/SP1.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/SP1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/SP1.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/SP1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/SP1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46234_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/SP1_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/SP1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46190_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/SP1_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46234_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_SP1_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46190_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_SP1_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45810_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45695_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/MXI1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/MXI1.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/MXI1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/MXI1.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/MXI1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/MXI1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45810_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/MXI1_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/MXI1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45695_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/MXI1_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45810_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_MXI1_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45695_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_MXI1_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45823_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45692_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/RFX5.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/RFX5.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/RFX5.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/RFX5.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/RFX5.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/RFX5.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45823_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/RFX5_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/RFX5.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45692_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/RFX5_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45823_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_RFX5_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45692_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_RFX5_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46235_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46181_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/HDAC2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/HDAC2.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/HDAC2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/HDAC2.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/HDAC2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/HDAC2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46235_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/HDAC2_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/HDAC2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46181_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/HDAC2_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46235_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_HDAC2_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46181_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_HDAC2_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/7301_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46174_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/ATF3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/ATF3.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/ATF3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/ATF3.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/ATF3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/ATF3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/7301_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/ATF3_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/ATF3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46174_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/ATF3_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/7301_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_ATF3_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46174_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_ATF3_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45815_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45687_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/MAFK.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/MAFK.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/MAFK.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/MAFK.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/MAFK.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/MAFK.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45815_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/MAFK_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/MAFK.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45687_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/MAFK_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45815_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_MAFK_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45687_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_MAFK_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45791_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45684_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/TBP.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/TBP.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/TBP.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/TBP.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/TBP.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/TBP.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45791_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/TBP_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/TBP.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45684_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/TBP_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45791_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_TBP_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45684_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_TBP_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74655_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45694_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/USF2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/USF2.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/USF2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/USF2.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/USF2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/USF2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74655_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/USF2_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/USF2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45694_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/USF2_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74655_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_USF2_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45694_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_USF2_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45802_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45701_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/CHD2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/CHD2.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/CHD2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/CHD2.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/CHD2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/CHD2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45802_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/CHD2_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/CHD2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45701_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/CHD2_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45802_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_CHD2_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45701_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_CHD2_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45809_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45698_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/BRCA1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/BRCA1.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/BRCA1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/BRCA1.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/BRCA1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/BRCA1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45809_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/BRCA1_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/BRCA1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45698_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/BRCA1_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45809_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_BRCA1_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45698_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_BRCA1_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46246_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46176_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/SP2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/SP2.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/SP2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/SP2.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/SP2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/SP2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46246_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/SP2_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/SP2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46176_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/SP2_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46246_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_SP2_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46176_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_SP2_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46242_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46192_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/RXRA.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/RXRA.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/RXRA.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_H1/RXRA.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_H1/RXRA.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/RXRA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46242_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/RXRA_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_H1/RXRA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46192_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_H1/RXRA_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46242_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_RXRA_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46192_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_H1/_RXRA_H1.bed


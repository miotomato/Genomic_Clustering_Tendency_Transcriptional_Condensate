#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=8000
#SBATCH -t 4:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f3_union_peaks_across_cellTypes/slurm_H1_vs_IMR90.out

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46199_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/41280_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/POLR2A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/POLR2A.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/POLR2A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/POLR2A.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/POLR2A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46199_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/POLR2A_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/41280_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/POLR2A_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46199_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_POLR2A_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/41280_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_POLR2A_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/908_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1834_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K4me3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K4me3.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K4me3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K4me3.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K4me3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/908_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K4me3_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1834_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K4me3_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/908_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_H3K4me3_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1834_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_H3K4me3_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/6531_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62068_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K27ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K27ac.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K27ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K27ac.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K27ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/6531_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K27ac_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62068_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K27ac_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/6531_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_H3K27ac_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62068_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_H3K27ac_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/50987_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45865_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/CTCF.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/CTCF.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/CTCF.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/CTCF.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/CTCF.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/50987_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/CTCF_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45865_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/CTCF_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/50987_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_CTCF_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45865_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_CTCF_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/897_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1832_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K4me2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K4me2.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K4me2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K4me2.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K4me2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/897_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K4me2_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1832_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K4me2_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/897_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_H3K4me2_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1832_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_H3K4me2_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46165_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/36878_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/EP300.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/EP300.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/EP300.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/EP300.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/EP300.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46165_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/EP300_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/36878_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/EP300_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46165_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_EP300_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/36878_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_EP300_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45227_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1845_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K9ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K9ac.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K9ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K9ac.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K9ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45227_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K9ac_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1845_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K9ac_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45227_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_H3K9ac_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1845_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_H3K9ac_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/894_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1829_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K4me1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K4me1.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K4me1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K4me1.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K4me1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K4me1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/894_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K4me1_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K4me1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1829_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K4me1_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/894_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_H3K4me1_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1829_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_H3K4me1_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46164_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45862_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/RAD21.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/RAD21.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/RAD21.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/RAD21.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/RAD21.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/RAD21.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46164_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/RAD21_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/RAD21.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45862_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/RAD21_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46164_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_RAD21_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45862_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_RAD21_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45681_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45861_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/CEBPB.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/CEBPB.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/CEBPB.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/CEBPB.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/CEBPB.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/CEBPB.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45681_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/CEBPB_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/CEBPB.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45861_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/CEBPB_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45681_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_CEBPB_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45861_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_CEBPB_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/6421_b_sort_peaks.broadPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/41279_b_sort_peaks.broadPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K36me3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K36me3.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K36me3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K36me3.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K36me3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K36me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/6421_b_sort_peaks.broadPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K36me3_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K36me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/41279_b_sort_peaks.broadPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K36me3_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/6421_b_sort_peaks.broadPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_H3K36me3_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/41279_b_sort_peaks.broadPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_H3K36me3_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45695_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45866_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/MXI1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/MXI1.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/MXI1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/MXI1.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/MXI1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/MXI1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45695_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/MXI1_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/MXI1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45866_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/MXI1_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45695_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_MXI1_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45866_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_MXI1_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45692_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45863_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/RFX5.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/RFX5.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/RFX5.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/RFX5.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/RFX5.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/RFX5.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45692_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/RFX5_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/RFX5.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45863_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/RFX5_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45692_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_RFX5_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45863_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_RFX5_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45687_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45864_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/MAFK.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/MAFK.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/MAFK.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/MAFK.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/MAFK.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/MAFK.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45687_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/MAFK_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/MAFK.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45864_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/MAFK_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45687_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_MAFK_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45864_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_MAFK_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/38503_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1828_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K4ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K4ac.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K4ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K4ac.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K4ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K4ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/38503_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K4ac_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K4ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1828_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K4ac_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/38503_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_H3K4ac_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1828_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_H3K4ac_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/38520_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1817_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K18ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K18ac.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K18ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K18ac.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K18ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K18ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/38520_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K18ac_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K18ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1817_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K18ac_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/38520_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_H3K18ac_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1817_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_H3K18ac_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/932_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1852_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H4K5ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/H4K5ac.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H4K5ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/H4K5ac.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H4K5ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/H4K5ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/932_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H4K5ac_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/H4K5ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1852_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H4K5ac_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/932_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_H4K5ac_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1852_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_H4K5ac_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/38525_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1814_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K14ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K14ac.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K14ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K14ac.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K14ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K14ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/38525_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K14ac_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K14ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1814_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K14ac_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/38525_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_H3K14ac_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1814_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_H3K14ac_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/6462_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1837_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K56ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K56ac.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K56ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K56ac.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K56ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K56ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/6462_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K56ac_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K56ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1837_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H3K56ac_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/6462_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_H3K56ac_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1837_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_H3K56ac_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/38515_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1804_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H2BK120ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/H2BK120ac.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H2BK120ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/H2BK120ac.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H2BK120ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/H2BK120ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/38515_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H2BK120ac_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/H2BK120ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1804_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H2BK120ac_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/38515_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_H2BK120ac_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1804_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_H2BK120ac_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/38498_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/6527_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H2BK5ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/H2BK5ac.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H2BK5ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/H2BK5ac.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H2BK5ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/H2BK5ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/38498_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H2BK5ac_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/H2BK5ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/6527_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H2BK5ac_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/38498_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_H2BK5ac_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/6527_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_H2BK5ac_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/38527_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1856_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H4K8ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/H4K8ac.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H4K8ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/H4K8ac.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H4K8ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/H4K8ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/38527_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H4K8ac_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/H4K8ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1856_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H4K8ac_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/38527_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_H4K8ac_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1856_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_H4K8ac_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/38506_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1857_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H4K91ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/H4K91ac.cat.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H4K91ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/H1_vs_IMR90/H4K91ac.sort.bed > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H4K91ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/H4K91ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/38506_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H4K91ac_merge_intersect_H1.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/H1_vs_IMR90/H4K91ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1857_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/H1_vs_IMR90/H4K91ac_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/38506_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_H4K91ac_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1857_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/H1_vs_IMR90/_H4K91ac_IMR90.bed


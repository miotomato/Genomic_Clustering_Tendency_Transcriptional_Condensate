#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=8000
#SBATCH -t 4:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f3_union_peaks_across_cellTypes/slurm_GM12878_vs_IMR90.out

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46137_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/41280_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/POLR2A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/POLR2A.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/POLR2A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/POLR2A.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/POLR2A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46137_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/POLR2A_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/41280_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/POLR2A_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46137_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/_POLR2A_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/41280_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/_POLR2A_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42291_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1834_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K4me3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K4me3.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K4me3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K4me3.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K4me3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42291_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K4me3_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1834_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K4me3_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42291_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/_H3K4me3_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1834_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/_H3K4me3_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/48285_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62068_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K27ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K27ac.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K27ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K27ac.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K27ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/48285_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K27ac_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62068_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K27ac_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/48285_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/_H3K27ac_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62068_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/_H3K27ac_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45631_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45865_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/CTCF.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/CTCF.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/CTCF.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/CTCF.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/CTCF.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45631_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/CTCF_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45865_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/CTCF_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45631_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/_CTCF_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45865_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/_CTCF_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45214_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1832_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K4me2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K4me2.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K4me2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K4me2.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K4me2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45214_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K4me2_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1832_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K4me2_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45214_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/_H3K4me2_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1832_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/_H3K4me2_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/48272_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/36878_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/EP300.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/EP300.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/EP300.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/EP300.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/EP300.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/48272_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/EP300_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/36878_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/EP300_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/48272_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/_EP300_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/36878_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/_EP300_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/36652_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1845_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K9ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K9ac.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K9ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K9ac.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K9ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/36652_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K9ac_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1845_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K9ac_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/36652_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/_H3K9ac_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1845_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/_H3K9ac_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/49297_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1829_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K4me1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K4me1.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K4me1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K4me1.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K4me1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K4me1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/49297_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K4me1_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K4me1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1829_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K4me1_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/49297_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/_H3K4me1_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1829_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/_H3K4me1_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46150_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45862_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/RAD21.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/RAD21.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/RAD21.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/RAD21.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/RAD21.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/RAD21.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46150_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/RAD21_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/RAD21.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45862_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/RAD21_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46150_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/_RAD21_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45862_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/_RAD21_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45221_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/6471_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H2AZ.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H2AZ.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H2AZ.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H2AZ.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H2AZ.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H2AZ.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45221_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H2AZ_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H2AZ.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/6471_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H2AZ_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45221_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/_H2AZ_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/6471_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/_H2AZ_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45640_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62504_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/SMC3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/SMC3.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/SMC3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/SMC3.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/SMC3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/SMC3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45640_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/SMC3_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/SMC3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62504_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/SMC3_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45640_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/_SMC3_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62504_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/_SMC3_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42284_b_sort_peaks.broadPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/41279_b_sort_peaks.broadPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K36me3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K36me3.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K36me3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K36me3.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K36me3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K36me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42284_b_sort_peaks.broadPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K36me3_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K36me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/41279_b_sort_peaks.broadPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H3K36me3_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42284_b_sort_peaks.broadPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/_H3K36me3_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/41279_b_sort_peaks.broadPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/_H3K36me3_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45633_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45866_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/MXI1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/MXI1.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/MXI1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/MXI1.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/MXI1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/MXI1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45633_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/MXI1_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/MXI1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45866_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/MXI1_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45633_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/_MXI1_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45866_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/_MXI1_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45649_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45863_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/RFX5.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/RFX5.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/RFX5.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/RFX5.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/RFX5.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/RFX5.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45649_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/RFX5_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/RFX5.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45863_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/RFX5_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45649_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/_RFX5_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45863_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/_RFX5_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45621_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45864_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/MAFK.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/MAFK.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/MAFK.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/MAFK.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/MAFK.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/MAFK.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45621_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/MAFK_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/MAFK.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45864_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/MAFK_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45621_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/_MAFK_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45864_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/_MAFK_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45661_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45867_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/MAZ.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/MAZ.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/MAZ.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/MAZ.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/MAZ.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/MAZ.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45661_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/MAZ_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/MAZ.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45867_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/MAZ_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45661_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/_MAZ_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45867_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/_MAZ_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45632_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64946_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/BHLHE40.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/BHLHE40.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/BHLHE40.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/BHLHE40.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/BHLHE40.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/BHLHE40.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45632_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/BHLHE40_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/BHLHE40.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64946_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/BHLHE40_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45632_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/_BHLHE40_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64946_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/_BHLHE40_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45652_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45869_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/RCOR1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/RCOR1.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/RCOR1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/RCOR1.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/RCOR1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/RCOR1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45652_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/RCOR1_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/RCOR1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45869_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/RCOR1_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45652_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/_RCOR1_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45869_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/_RCOR1_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/67994_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1812_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H2BK20ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H2BK20ac.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H2BK20ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H2BK20ac.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H2BK20ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H2BK20ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/67994_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H2BK20ac_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H2BK20ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1812_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/H2BK20ac_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/67994_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/_H2BK20ac_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1812_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_IMR90/_H2BK20ac_IMR90.bed


#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=8000
#SBATCH -t 4:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f3_union_peaks_across_cellTypes/slurm_HeLa_vs_IMR90.out

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/73296_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/41280_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/POLR2A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/POLR2A.cat.bed > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/POLR2A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/POLR2A.sort.bed > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/POLR2A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/73296_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/POLR2A_merge_intersect_HeLa.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/41280_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/POLR2A_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/73296_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/_POLR2A_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/41280_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/_POLR2A_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/53321_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1834_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/H3K4me3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/H3K4me3.cat.bed > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/H3K4me3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/H3K4me3.sort.bed > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/H3K4me3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/53321_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/H3K4me3_merge_intersect_HeLa.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1834_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/H3K4me3_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/53321_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/_H3K4me3_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1834_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/_H3K4me3_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62017_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62068_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/H3K27ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/H3K27ac.cat.bed > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/H3K27ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/H3K27ac.sort.bed > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/H3K27ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62017_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/H3K27ac_merge_intersect_HeLa.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62068_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/H3K27ac_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62017_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/_H3K27ac_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62068_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/_H3K27ac_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8087_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45865_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/CTCF.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/CTCF.cat.bed > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/CTCF.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/CTCF.sort.bed > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/CTCF.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8087_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/CTCF_merge_intersect_HeLa.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45865_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/CTCF_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8087_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/_CTCF_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45865_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/_CTCF_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1252_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1832_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/H3K4me2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/H3K4me2.cat.bed > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/H3K4me2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/H3K4me2.sort.bed > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/H3K4me2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1252_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/H3K4me2_merge_intersect_HeLa.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1832_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/H3K4me2_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1252_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/_H3K4me2_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1832_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/_H3K4me2_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56234_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/36878_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/EP300.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/EP300.cat.bed > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/EP300.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/EP300.sort.bed > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/EP300.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56234_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/EP300_merge_intersect_HeLa.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/36878_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/EP300_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56234_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/_EP300_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/36878_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/_EP300_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/67816_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1845_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/H3K9ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/H3K9ac.cat.bed > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/H3K9ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/H3K9ac.sort.bed > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/H3K9ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/67816_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/H3K9ac_merge_intersect_HeLa.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1845_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/H3K9ac_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/67816_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/_H3K9ac_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1845_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/_H3K9ac_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1231_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/6471_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/H2AZ.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/H2AZ.cat.bed > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/H2AZ.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/H2AZ.sort.bed > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/H2AZ.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/H2AZ.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1231_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/H2AZ_merge_intersect_HeLa.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/H2AZ.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/6471_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/H2AZ_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1231_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/_H2AZ_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/6471_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/_H2AZ_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/87043_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62504_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/SMC3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/SMC3.cat.bed > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/SMC3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/SMC3.sort.bed > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/SMC3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/SMC3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/87043_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/SMC3_merge_intersect_HeLa.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/SMC3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62504_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/SMC3_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/87043_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/_SMC3_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62504_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/_SMC3_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56961_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/36362_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/GTF2B.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/GTF2B.cat.bed > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/GTF2B.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/GTF2B.sort.bed > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/GTF2B.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/GTF2B.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56961_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/GTF2B_merge_intersect_HeLa.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/GTF2B.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/36362_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/GTF2B_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56961_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/_GTF2B_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/36362_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/_GTF2B_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/9216_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/36356_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/POLR2B.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/POLR2B.cat.bed > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/POLR2B.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/POLR2B.sort.bed > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/POLR2B.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/POLR2B.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/9216_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/POLR2B_merge_intersect_HeLa.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/POLR2B.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/36356_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/POLR2B_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/9216_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/_POLR2B_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/36356_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HeLa_vs_IMR90/_POLR2B_IMR90.bed


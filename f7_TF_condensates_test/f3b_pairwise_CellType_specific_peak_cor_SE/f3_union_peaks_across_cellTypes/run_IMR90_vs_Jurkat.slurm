#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=8000
#SBATCH -t 4:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f3_union_peaks_across_cellTypes/slurm_IMR90_vs_Jurkat.out

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/41280_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/5199_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/POLR2A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/POLR2A.cat.bed > f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/POLR2A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/POLR2A.sort.bed > f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/POLR2A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/41280_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/POLR2A_merge_intersect_IMR90.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/5199_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/POLR2A_merge_intersect_Jurkat.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/41280_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/_POLR2A_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/5199_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/_POLR2A_Jurkat.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1834_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/49527_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/H3K4me3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/H3K4me3.cat.bed > f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/H3K4me3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/H3K4me3.sort.bed > f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/H3K4me3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1834_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/H3K4me3_merge_intersect_IMR90.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/49527_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/H3K4me3_merge_intersect_Jurkat.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1834_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/_H3K4me3_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/49527_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/_H3K4me3_Jurkat.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62068_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/54237_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/H3K27ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/H3K27ac.cat.bed > f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/H3K27ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/H3K27ac.sort.bed > f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/H3K27ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62068_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/H3K27ac_merge_intersect_IMR90.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/54237_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/H3K27ac_merge_intersect_Jurkat.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62068_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/_H3K27ac_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/54237_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/_H3K27ac_Jurkat.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45865_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/54151_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/CTCF.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/CTCF.cat.bed > f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/CTCF.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/CTCF.sort.bed > f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/CTCF.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45865_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/CTCF_merge_intersect_IMR90.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/54151_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/CTCF_merge_intersect_Jurkat.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45865_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/_CTCF_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/54151_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/_CTCF_Jurkat.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1832_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/74134_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/H3K4me2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/H3K4me2.cat.bed > f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/H3K4me2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/H3K4me2.sort.bed > f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/H3K4me2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1832_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/H3K4me2_merge_intersect_IMR90.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/74134_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/H3K4me2_merge_intersect_Jurkat.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1832_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/_H3K4me2_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/74134_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/_H3K4me2_Jurkat.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/68711_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/69617_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/BRD4.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/BRD4.cat.bed > f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/BRD4.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/BRD4.sort.bed > f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/BRD4.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/BRD4.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/68711_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/BRD4_merge_intersect_IMR90.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/BRD4.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/69617_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/BRD4_merge_intersect_Jurkat.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/68711_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/_BRD4_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/69617_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/IMR90_vs_Jurkat/_BRD4_Jurkat.bed


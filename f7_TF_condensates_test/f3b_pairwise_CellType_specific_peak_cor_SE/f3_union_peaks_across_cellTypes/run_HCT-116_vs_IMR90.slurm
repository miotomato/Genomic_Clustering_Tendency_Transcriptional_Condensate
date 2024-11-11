#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=8000
#SBATCH -t 4:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f3_union_peaks_across_cellTypes/slurm_HCT-116_vs_IMR90.out

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46201_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/41280_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/POLR2A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/POLR2A.cat.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/POLR2A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/POLR2A.sort.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/POLR2A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46201_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/POLR2A_merge_intersect_HCT-116.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/41280_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/POLR2A_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46201_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/_POLR2A_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/41280_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/_POLR2A_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57094_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1834_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H3K4me3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H3K4me3.cat.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H3K4me3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H3K4me3.sort.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H3K4me3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57094_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H3K4me3_merge_intersect_HCT-116.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1834_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H3K4me3_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57094_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/_H3K4me3_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1834_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/_H3K4me3_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62110_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62068_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H3K27ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H3K27ac.cat.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H3K27ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H3K27ac.sort.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H3K27ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62110_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H3K27ac_merge_intersect_HCT-116.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62068_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H3K27ac_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62110_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/_H3K27ac_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62068_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/_H3K27ac_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/42149_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45865_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/CTCF.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/CTCF.cat.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/CTCF.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/CTCF.sort.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/CTCF.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/42149_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/CTCF_merge_intersect_HCT-116.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45865_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/CTCF_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/42149_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/_CTCF_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45865_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/_CTCF_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42910_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1832_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H3K4me2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H3K4me2.cat.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H3K4me2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H3K4me2.sort.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H3K4me2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42910_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H3K4me2_merge_intersect_HCT-116.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1832_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H3K4me2_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42910_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/_H3K4me2_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1832_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/_H3K4me2_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/42907_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/36878_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/EP300.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/EP300.cat.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/EP300.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/EP300.sort.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/EP300.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/42907_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/EP300_merge_intersect_HCT-116.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/36878_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/EP300_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/42907_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/_EP300_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/36878_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/_EP300_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/87280_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1845_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H3K9ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H3K9ac.cat.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H3K9ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H3K9ac.sort.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H3K9ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/87280_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H3K9ac_merge_intersect_HCT-116.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1845_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H3K9ac_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/87280_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/_H3K9ac_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1845_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/_H3K9ac_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57095_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1829_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H3K4me1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H3K4me1.cat.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H3K4me1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H3K4me1.sort.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H3K4me1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H3K4me1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57095_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H3K4me1_merge_intersect_HCT-116.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H3K4me1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1829_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H3K4me1_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57095_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/_H3K4me1_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1829_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/_H3K4me1_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46207_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45862_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/RAD21.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/RAD21.cat.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/RAD21.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/RAD21.sort.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/RAD21.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/RAD21.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46207_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/RAD21_merge_intersect_HCT-116.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/RAD21.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45862_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/RAD21_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46207_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/_RAD21_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45862_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/_RAD21_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/83065_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/6471_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H2AZ.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H2AZ.cat.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H2AZ.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H2AZ.sort.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H2AZ.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H2AZ.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/83065_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H2AZ_merge_intersect_HCT-116.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H2AZ.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/6471_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/H2AZ_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/83065_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/_H2AZ_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/6471_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/_H2AZ_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46206_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45861_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/CEBPB.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/CEBPB.cat.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/CEBPB.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/CEBPB.sort.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/CEBPB.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/CEBPB.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46206_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/CEBPB_merge_intersect_HCT-116.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/CEBPB.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45861_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/CEBPB_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46206_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/_CEBPB_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45861_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/_CEBPB_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/68783_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/68711_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/BRD4.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/BRD4.cat.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/BRD4.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/BRD4.sort.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/BRD4.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/BRD4.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/68783_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/BRD4_merge_intersect_HCT-116.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/BRD4.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/68711_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/BRD4_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/68783_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/_BRD4_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/68711_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/_BRD4_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/82544_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/52092_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/TP53.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/TP53.cat.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/TP53.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/TP53.sort.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/TP53.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/TP53.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/82544_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/TP53_merge_intersect_HCT-116.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/TP53.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/52092_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/TP53_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/82544_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/_TP53_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/52092_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/_TP53_IMR90.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/34399_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/34392_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/MECP2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/MECP2.cat.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/MECP2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/MECP2.sort.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/MECP2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/MECP2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/34399_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/MECP2_merge_intersect_HCT-116.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/MECP2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/34392_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/MECP2_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/34399_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/_MECP2_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/34392_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HCT-116_vs_IMR90/_MECP2_IMR90.bed


#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=8000
#SBATCH -t 4:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f3_union_peaks_across_cellTypes/slurm_K562_vs_IMR90.out

echo POLR2A 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/33859_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/41280_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/POLR2A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_IMR90/POLR2A.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/POLR2A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_IMR90/POLR2A.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/POLR2A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_IMR90/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/33859_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_IMR90/POLR2A_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_IMR90/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/41280_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_IMR90/POLR2A_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/33859_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_IMR90/_POLR2A_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/41280_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_IMR90/_POLR2A_IMR90.bed

echo H3K4me3 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/72784_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1834_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3K4me3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3K4me3.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3K4me3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3K4me3.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3K4me3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/72784_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3K4me3_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1834_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3K4me3_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/72784_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_IMR90/_H3K4me3_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1834_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_IMR90/_H3K4me3_IMR90.bed

echo H3K27ac 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/35400_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62068_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3K27ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3K27ac.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3K27ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3K27ac.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3K27ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/35400_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3K27ac_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62068_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3K27ac_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/35400_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_IMR90/_H3K27ac_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62068_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_IMR90/_H3K27ac_IMR90.bed

echo CTCF 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45899_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45865_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/CTCF.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_IMR90/CTCF.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/CTCF.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_IMR90/CTCF.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/CTCF.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_IMR90/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45899_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_IMR90/CTCF_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_IMR90/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45865_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_IMR90/CTCF_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45899_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_IMR90/_CTCF_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45865_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_IMR90/_CTCF_IMR90.bed

echo H3K4me2 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45382_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1832_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3K4me2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3K4me2.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3K4me2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3K4me2.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3K4me2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45382_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3K4me2_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1832_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3K4me2_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45382_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_IMR90/_H3K4me2_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1832_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_IMR90/_H3K4me2_IMR90.bed

echo EP300 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45963_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/36878_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/EP300.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_IMR90/EP300.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/EP300.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_IMR90/EP300.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/EP300.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_IMR90/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45963_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_IMR90/EP300_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_IMR90/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/36878_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_IMR90/EP300_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45963_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_IMR90/_EP300_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/36878_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_IMR90/_EP300_IMR90.bed

echo H3K9ac 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45406_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1845_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3K9ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3K9ac.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3K9ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3K9ac.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3K9ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45406_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3K9ac_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1845_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3K9ac_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45406_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_IMR90/_H3K9ac_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1845_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_IMR90/_H3K9ac_IMR90.bed

echo H3K4me1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/8462_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1829_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3K4me1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3K4me1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3K4me1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3K4me1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3K4me1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3K4me1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/8462_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3K4me1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3K4me1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1829_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3K4me1_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/8462_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_IMR90/_H3K4me1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1829_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_IMR90/_H3K4me1_IMR90.bed

echo RAD21 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46275_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45862_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/RAD21.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_IMR90/RAD21.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/RAD21.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_IMR90/RAD21.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/RAD21.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_IMR90/RAD21.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46275_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_IMR90/RAD21_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_IMR90/RAD21.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45862_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_IMR90/RAD21_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46275_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_IMR90/_RAD21_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45862_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_IMR90/_RAD21_IMR90.bed

echo H2AZ 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45416_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/6471_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/H2AZ.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_IMR90/H2AZ.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/H2AZ.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_IMR90/H2AZ.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/H2AZ.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_IMR90/H2AZ.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45416_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_IMR90/H2AZ_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_IMR90/H2AZ.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/6471_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_IMR90/H2AZ_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45416_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_IMR90/_H2AZ_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/6471_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_IMR90/_H2AZ_IMR90.bed

echo CEBPB 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45923_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45861_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/CEBPB.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_IMR90/CEBPB.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/CEBPB.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_IMR90/CEBPB.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/CEBPB.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_IMR90/CEBPB.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45923_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_IMR90/CEBPB_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_IMR90/CEBPB.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45861_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_IMR90/CEBPB_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45923_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_IMR90/_CEBPB_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45861_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_IMR90/_CEBPB_IMR90.bed

echo SMC3 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45914_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62504_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/SMC3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_IMR90/SMC3.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/SMC3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_IMR90/SMC3.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/SMC3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_IMR90/SMC3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45914_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_IMR90/SMC3_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_IMR90/SMC3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62504_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_IMR90/SMC3_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45914_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_IMR90/_SMC3_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62504_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_IMR90/_SMC3_IMR90.bed

echo MXI1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45879_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45866_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/MXI1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_IMR90/MXI1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/MXI1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_IMR90/MXI1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/MXI1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_IMR90/MXI1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45879_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_IMR90/MXI1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_IMR90/MXI1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45866_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_IMR90/MXI1_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45879_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_IMR90/_MXI1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45866_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_IMR90/_MXI1_IMR90.bed

echo RFX5 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45893_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45863_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/RFX5.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_IMR90/RFX5.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/RFX5.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_IMR90/RFX5.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/RFX5.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_IMR90/RFX5.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45893_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_IMR90/RFX5_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_IMR90/RFX5.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45863_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_IMR90/RFX5_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45893_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_IMR90/_RFX5_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45863_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_IMR90/_RFX5_IMR90.bed

echo MAFK 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45954_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45864_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/MAFK.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_IMR90/MAFK.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/MAFK.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_IMR90/MAFK.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/MAFK.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_IMR90/MAFK.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45954_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_IMR90/MAFK_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_IMR90/MAFK.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45864_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_IMR90/MAFK_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45954_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_IMR90/_MAFK_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45864_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_IMR90/_MAFK_IMR90.bed

echo MAZ 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45924_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45867_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/MAZ.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_IMR90/MAZ.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/MAZ.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_IMR90/MAZ.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/MAZ.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_IMR90/MAZ.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45924_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_IMR90/MAZ_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_IMR90/MAZ.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45867_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_IMR90/MAZ_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45924_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_IMR90/_MAZ_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45867_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_IMR90/_MAZ_IMR90.bed

echo BHLHE40 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45916_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64946_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/BHLHE40.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_IMR90/BHLHE40.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/BHLHE40.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_IMR90/BHLHE40.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/BHLHE40.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_IMR90/BHLHE40.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45916_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_IMR90/BHLHE40_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_IMR90/BHLHE40.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64946_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_IMR90/BHLHE40_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45916_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_IMR90/_BHLHE40_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64946_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_IMR90/_BHLHE40_IMR90.bed

echo RCOR1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45903_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45869_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/RCOR1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_IMR90/RCOR1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/RCOR1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_IMR90/RCOR1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/RCOR1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_IMR90/RCOR1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45903_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_IMR90/RCOR1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_IMR90/RCOR1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45869_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_IMR90/RCOR1_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45903_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_IMR90/_RCOR1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45869_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_IMR90/_RCOR1_IMR90.bed

echo GTF2B 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45932_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/36362_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/GTF2B.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_IMR90/GTF2B.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/GTF2B.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_IMR90/GTF2B.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/GTF2B.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_IMR90/GTF2B.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45932_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_IMR90/GTF2B_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_IMR90/GTF2B.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/36362_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_IMR90/GTF2B_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45932_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_IMR90/_GTF2B_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/36362_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_IMR90/_GTF2B_IMR90.bed

echo H3 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57312_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/35491_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57312_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/35491_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_IMR90/H3_merge_intersect_IMR90.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57312_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_IMR90/_H3_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/35491_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_IMR90/_H3_IMR90.bed


#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=8000
#SBATCH -t 4:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f3_union_peaks_across_cellTypes/slurm_K562_vs_H1.out

echo POLR2A 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/33859_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46199_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/POLR2A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/POLR2A.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/POLR2A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/POLR2A.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/POLR2A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/33859_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/POLR2A_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46199_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/POLR2A_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/33859_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_POLR2A_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46199_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_POLR2A_H1.bed

echo H3K4me3 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/72784_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/908_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/H3K4me3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/H3K4me3.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/H3K4me3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/H3K4me3.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/H3K4me3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/72784_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/H3K4me3_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/908_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/H3K4me3_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/72784_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_H3K4me3_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/908_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_H3K4me3_H1.bed

echo H3K27ac 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/35400_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/6531_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/H3K27ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/H3K27ac.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/H3K27ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/H3K27ac.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/H3K27ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/35400_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/H3K27ac_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/6531_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/H3K27ac_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/35400_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_H3K27ac_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/6531_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_H3K27ac_H1.bed

echo CTCF 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45899_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/50987_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/CTCF.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/CTCF.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/CTCF.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/CTCF.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/CTCF.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45899_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/CTCF_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/50987_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/CTCF_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45899_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_CTCF_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/50987_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_CTCF_H1.bed

echo H3K4me2 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45382_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/897_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/H3K4me2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/H3K4me2.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/H3K4me2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/H3K4me2.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/H3K4me2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45382_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/H3K4me2_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/897_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/H3K4me2_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45382_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_H3K4me2_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/897_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_H3K4me2_H1.bed

echo EP300 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45963_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46165_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/EP300.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/EP300.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/EP300.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/EP300.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/EP300.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45963_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/EP300_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46165_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/EP300_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45963_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_EP300_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46165_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_EP300_H1.bed

echo H3K9ac 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45406_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45227_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/H3K9ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/H3K9ac.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/H3K9ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/H3K9ac.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/H3K9ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45406_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/H3K9ac_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45227_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/H3K9ac_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45406_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_H3K9ac_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45227_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_H3K9ac_H1.bed

echo MYC 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45894_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45691_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/MYC.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/MYC.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/MYC.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/MYC.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/MYC.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45894_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/MYC_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45691_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/MYC_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45894_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_MYC_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45691_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_MYC_H1.bed

echo H3K4me1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/8462_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/894_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/H3K4me1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/H3K4me1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/H3K4me1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/H3K4me1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/H3K4me1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/H3K4me1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/8462_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/H3K4me1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/H3K4me1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/894_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/H3K4me1_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/8462_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_H3K4me1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/894_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_H3K4me1_H1.bed

echo REST 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46286_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46172_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/REST.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/REST.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/REST.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/REST.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/REST.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/REST.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46286_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/REST_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/REST.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46172_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/REST_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46286_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_REST_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46172_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_REST_H1.bed

echo GABPA 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46268_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46185_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/GABPA.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/GABPA.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/GABPA.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/GABPA.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/GABPA.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/GABPA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46268_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/GABPA_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/GABPA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46185_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/GABPA_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46268_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_GABPA_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46185_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_GABPA_H1.bed

echo YY1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46303_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46170_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/YY1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/YY1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/YY1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/YY1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/YY1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/YY1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46303_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/YY1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/YY1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46170_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/YY1_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46303_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_YY1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46170_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_YY1_H1.bed

echo MAX 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46292_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46180_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/MAX.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/MAX.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/MAX.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/MAX.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/MAX.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/MAX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46292_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/MAX_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/MAX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46180_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/MAX_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46292_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_MAX_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46180_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_MAX_H1.bed

echo RAD21 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46275_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46164_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/RAD21.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/RAD21.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/RAD21.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/RAD21.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/RAD21.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/RAD21.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46275_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/RAD21_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/RAD21.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46164_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/RAD21_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46275_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_RAD21_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46164_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_RAD21_H1.bed

echo TAF1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46282_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46196_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/TAF1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/TAF1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/TAF1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/TAF1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/TAF1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/TAF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46282_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/TAF1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/TAF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46196_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/TAF1_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46282_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_TAF1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46196_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_TAF1_H1.bed

echo CEBPB 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45923_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45681_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/CEBPB.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/CEBPB.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/CEBPB.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/CEBPB.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/CEBPB.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/CEBPB.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45923_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/CEBPB_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/CEBPB.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45681_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/CEBPB_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45923_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_CEBPB_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45681_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_CEBPB_H1.bed

echo JUND 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45880_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45696_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/JUND.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/JUND.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/JUND.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/JUND.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/JUND.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/JUND.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45880_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/JUND_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/JUND.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45696_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/JUND_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45880_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_JUND_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45696_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_JUND_H1.bed

echo H3K27me3 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57482_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/8501_b_sort_peaks.broadPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/H3K27me3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/H3K27me3.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/H3K27me3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/H3K27me3.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/H3K27me3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/H3K27me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57482_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/H3K27me3_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/H3K27me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/8501_b_sort_peaks.broadPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/H3K27me3_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57482_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_H3K27me3_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/8501_b_sort_peaks.broadPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_H3K27me3_H1.bed

echo ZNF143 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45946_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45688_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/ZNF143.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/ZNF143.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/ZNF143.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/ZNF143.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/ZNF143.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/ZNF143.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45946_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/ZNF143_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/ZNF143.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45688_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/ZNF143_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45946_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_ZNF143_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45688_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_ZNF143_H1.bed

echo TEAD4 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46284_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46189_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/TEAD4.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/TEAD4.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/TEAD4.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/TEAD4.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/TEAD4.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/TEAD4.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46284_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/TEAD4_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/TEAD4.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46189_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/TEAD4_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46284_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_TEAD4_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46189_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_TEAD4_H1.bed

echo H3K9me3 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/56869_b_sort_peaks.broadPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/48924_b_sort_peaks.broadPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/H3K9me3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/H3K9me3.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/H3K9me3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/H3K9me3.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/H3K9me3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/H3K9me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/56869_b_sort_peaks.broadPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/H3K9me3_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/H3K9me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/48924_b_sort_peaks.broadPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/H3K9me3_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/56869_b_sort_peaks.broadPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_H3K9me3_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/48924_b_sort_peaks.broadPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_H3K9me3_H1.bed

echo NRF1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74666_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45697_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/NRF1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/NRF1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/NRF1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/NRF1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/NRF1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/NRF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74666_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/NRF1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/NRF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45697_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/NRF1_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74666_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_NRF1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45697_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_NRF1_H1.bed

echo SRF 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46299_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46167_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/SRF.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/SRF.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/SRF.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/SRF.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/SRF.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/SRF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46299_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/SRF_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/SRF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46167_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/SRF_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46299_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_SRF_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46167_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_SRF_H1.bed

echo CREB1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46274_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46182_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/CREB1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/CREB1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/CREB1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/CREB1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/CREB1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/CREB1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46274_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/CREB1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/CREB1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46182_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/CREB1_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46274_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_CREB1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46182_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_CREB1_H1.bed

echo JUN 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2078_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45685_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/JUN.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/JUN.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/JUN.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/JUN.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/JUN.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/JUN.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2078_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/JUN_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/JUN.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45685_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/JUN_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2078_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_JUN_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45685_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_JUN_H1.bed

echo USF1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46279_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46171_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/USF1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/USF1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/USF1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/USF1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/USF1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/USF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46279_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/USF1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/USF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46171_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/USF1_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46279_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_USF1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46171_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_USF1_H1.bed

echo SIN3A 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46277_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45699_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/SIN3A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/SIN3A.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/SIN3A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/SIN3A.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/SIN3A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/SIN3A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46277_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/SIN3A_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/SIN3A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45699_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/SIN3A_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46277_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_SIN3A_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45699_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_SIN3A_H1.bed

echo SP1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46300_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46190_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/SP1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/SP1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/SP1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/SP1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/SP1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/SP1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46300_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/SP1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/SP1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46190_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/SP1_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46300_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_SP1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46190_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_SP1_H1.bed

echo MXI1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45879_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45695_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/MXI1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/MXI1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/MXI1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/MXI1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/MXI1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/MXI1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45879_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/MXI1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/MXI1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45695_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/MXI1_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45879_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_MXI1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45695_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_MXI1_H1.bed

echo RFX5 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45893_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45692_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/RFX5.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/RFX5.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/RFX5.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/RFX5.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/RFX5.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/RFX5.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45893_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/RFX5_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/RFX5.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45692_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/RFX5_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45893_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_RFX5_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45692_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_RFX5_H1.bed

echo HDAC2 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46304_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46181_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/HDAC2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/HDAC2.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/HDAC2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/HDAC2.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/HDAC2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/HDAC2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46304_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/HDAC2_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/HDAC2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46181_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/HDAC2_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46304_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_HDAC2_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46181_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_HDAC2_H1.bed

echo ATF3 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64077_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46174_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/ATF3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/ATF3.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/ATF3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/ATF3.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/ATF3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/ATF3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64077_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/ATF3_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/ATF3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46174_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/ATF3_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64077_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_ATF3_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46174_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_ATF3_H1.bed

echo MAFK 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45954_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45687_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/MAFK.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/MAFK.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/MAFK.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/MAFK.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/MAFK.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/MAFK.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45954_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/MAFK_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/MAFK.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45687_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/MAFK_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45954_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_MAFK_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45687_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_MAFK_H1.bed

echo TBP 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45881_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45684_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/TBP.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/TBP.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/TBP.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/TBP.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/TBP.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/TBP.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45881_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/TBP_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/TBP.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45684_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/TBP_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45881_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_TBP_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45684_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_TBP_H1.bed

echo EGR1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46298_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46169_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/EGR1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/EGR1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/EGR1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/EGR1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/EGR1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/EGR1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46298_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/EGR1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/EGR1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46169_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/EGR1_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46298_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_EGR1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46169_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_EGR1_H1.bed

echo USF2 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63383_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45694_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/USF2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/USF2.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/USF2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/USF2.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/USF2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/USF2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63383_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/USF2_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/USF2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45694_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/USF2_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63383_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_USF2_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45694_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_USF2_H1.bed

echo CHD2 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45905_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45701_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/CHD2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/CHD2.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/CHD2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/CHD2.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/CHD2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/CHD2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45905_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/CHD2_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/CHD2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45701_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/CHD2_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45905_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_CHD2_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45701_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_CHD2_H1.bed

echo SIX5 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46302_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46188_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/SIX5.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/SIX5.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/SIX5.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/SIX5.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/SIX5.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/SIX5.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46302_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/SIX5_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/SIX5.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46188_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/SIX5_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46302_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_SIX5_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46188_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_SIX5_H1.bed

echo FOSL1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63076_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46175_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/FOSL1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/FOSL1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/FOSL1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/FOSL1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/FOSL1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/FOSL1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63076_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/FOSL1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/FOSL1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46175_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/FOSL1_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63076_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_FOSL1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46175_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_FOSL1_H1.bed

echo E2F6 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46297_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46177_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/E2F6.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/E2F6.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/E2F6.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/E2F6.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/E2F6.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/E2F6.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46297_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/E2F6_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/E2F6.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46177_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/E2F6_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46297_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_E2F6_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46177_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_E2F6_H1.bed

echo PHF8 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8790_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8818_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/PHF8.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/PHF8.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/PHF8.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/PHF8.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/PHF8.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/PHF8.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8790_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/PHF8_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/PHF8.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8818_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/PHF8_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8790_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_PHF8_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8818_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_PHF8_H1.bed

echo SP2 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46287_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46176_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/SP2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/SP2.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/SP2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/SP2.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/SP2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/SP2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46287_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/SP2_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/SP2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46176_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/SP2_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46287_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_SP2_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46176_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_SP2_H1.bed

echo TRIM28 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46305_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/48923_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/TRIM28.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/TRIM28.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/TRIM28.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/TRIM28.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/TRIM28.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/TRIM28.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46305_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/TRIM28_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/TRIM28.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/48923_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/TRIM28_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46305_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_TRIM28_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/48923_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_TRIM28_H1.bed

echo BACH1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45900_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45690_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/BACH1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/BACH1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/BACH1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/BACH1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/BACH1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/BACH1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45900_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/BACH1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/BACH1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45690_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/BACH1_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45900_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_BACH1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45690_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_BACH1_H1.bed

echo H3K79me2 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/37968_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/38517_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/H3K79me2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/H3K79me2.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/H3K79me2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/H3K79me2.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/H3K79me2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/H3K79me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/37968_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/H3K79me2_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/H3K79me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/38517_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/H3K79me2_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/37968_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_H3K79me2_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/38517_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_H3K79me2_H1.bed

echo TAF7 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46289_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46195_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/TAF7.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/TAF7.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/TAF7.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/TAF7.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/TAF7.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/TAF7.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46289_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/TAF7_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/TAF7.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46195_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/TAF7_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46289_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_TAF7_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46195_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_TAF7_H1.bed

echo RBBP5 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8795_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8820_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/RBBP5.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/RBBP5.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/RBBP5.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/RBBP5.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/RBBP5.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/RBBP5.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8795_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/RBBP5_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/RBBP5.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8820_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/RBBP5_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8795_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_RBBP5_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8820_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_RBBP5_H1.bed

echo SAP30 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8800_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8821_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/SAP30.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_H1/SAP30.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/SAP30.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_H1/SAP30.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_H1/SAP30.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/SAP30.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8800_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/SAP30_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_H1/SAP30.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8821_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_H1/SAP30_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8800_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_SAP30_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8821_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_H1/_SAP30_H1.bed


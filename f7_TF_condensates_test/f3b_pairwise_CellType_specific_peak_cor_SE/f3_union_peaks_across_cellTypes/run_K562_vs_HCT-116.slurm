#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=8000
#SBATCH -t 4:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f3_union_peaks_across_cellTypes/slurm_K562_vs_HCT-116.out

echo POLR2A 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/33859_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46201_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/POLR2A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/POLR2A.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/POLR2A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/POLR2A.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/POLR2A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/33859_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/POLR2A_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46201_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/POLR2A_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/33859_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_POLR2A_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46201_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_POLR2A_HCT-116.bed

echo H3K4me3 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/72784_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57094_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H3K4me3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H3K4me3.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H3K4me3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H3K4me3.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H3K4me3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/72784_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H3K4me3_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57094_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H3K4me3_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/72784_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_H3K4me3_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57094_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_H3K4me3_HCT-116.bed

echo H3K27ac 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/35400_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62110_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H3K27ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H3K27ac.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H3K27ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H3K27ac.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H3K27ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/35400_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H3K27ac_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62110_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H3K27ac_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/35400_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_H3K27ac_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62110_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_H3K27ac_HCT-116.bed

echo CTCF 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45899_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/42149_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/CTCF.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/CTCF.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/CTCF.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/CTCF.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/CTCF.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45899_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/CTCF_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/42149_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/CTCF_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45899_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_CTCF_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/42149_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_CTCF_HCT-116.bed

echo H3K4me2 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45382_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42910_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H3K4me2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H3K4me2.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H3K4me2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H3K4me2.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H3K4me2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45382_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H3K4me2_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42910_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H3K4me2_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45382_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_H3K4me2_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42910_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_H3K4me2_HCT-116.bed

echo EP300 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45963_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/42907_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/EP300.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/EP300.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/EP300.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/EP300.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/EP300.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45963_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/EP300_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/42907_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/EP300_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45963_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_EP300_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/42907_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_EP300_HCT-116.bed

echo H3K9ac 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45406_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/87280_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H3K9ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H3K9ac.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H3K9ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H3K9ac.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H3K9ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45406_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H3K9ac_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/87280_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H3K9ac_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45406_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_H3K9ac_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/87280_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_H3K9ac_HCT-116.bed

echo MYC 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45894_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/70806_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/MYC.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/MYC.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/MYC.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/MYC.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/MYC.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45894_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/MYC_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/70806_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/MYC_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45894_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_MYC_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/70806_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_MYC_HCT-116.bed

echo H3K4me1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/8462_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57095_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H3K4me1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H3K4me1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H3K4me1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H3K4me1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H3K4me1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H3K4me1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/8462_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H3K4me1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H3K4me1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57095_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H3K4me1_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/8462_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_H3K4me1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57095_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_H3K4me1_HCT-116.bed

echo REST 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46286_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46203_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/REST.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/REST.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/REST.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/REST.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/REST.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/REST.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46286_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/REST_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/REST.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46203_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/REST_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46286_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_REST_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46203_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_REST_HCT-116.bed

echo YY1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46303_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46204_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/YY1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/YY1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/YY1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/YY1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/YY1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/YY1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46303_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/YY1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/YY1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46204_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/YY1_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46303_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_YY1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46204_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_YY1_HCT-116.bed

echo MAX 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46292_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46216_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/MAX.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/MAX.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/MAX.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/MAX.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/MAX.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/MAX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46292_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/MAX_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/MAX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46216_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/MAX_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46292_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_MAX_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46216_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_MAX_HCT-116.bed

echo RAD21 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46275_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46207_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/RAD21.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/RAD21.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/RAD21.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/RAD21.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/RAD21.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/RAD21.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46275_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/RAD21_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/RAD21.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46207_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/RAD21_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46275_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_RAD21_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46207_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_RAD21_HCT-116.bed

echo H2AZ 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45416_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/83065_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H2AZ.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H2AZ.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H2AZ.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H2AZ.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H2AZ.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H2AZ.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45416_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H2AZ_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H2AZ.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/83065_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/H2AZ_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45416_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_H2AZ_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/83065_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_H2AZ_HCT-116.bed

echo CEBPB 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45923_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46206_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/CEBPB.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/CEBPB.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/CEBPB.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/CEBPB.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/CEBPB.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/CEBPB.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45923_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/CEBPB_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/CEBPB.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46206_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/CEBPB_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45923_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_CEBPB_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46206_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_CEBPB_HCT-116.bed

echo JUND 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45880_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46213_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/JUND.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/JUND.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/JUND.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/JUND.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/JUND.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/JUND.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45880_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/JUND_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/JUND.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46213_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/JUND_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45880_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_JUND_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46213_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_JUND_HCT-116.bed

echo TEAD4 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46284_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46210_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/TEAD4.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/TEAD4.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/TEAD4.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/TEAD4.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/TEAD4.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/TEAD4.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46284_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/TEAD4_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/TEAD4.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46210_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/TEAD4_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46284_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_TEAD4_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46210_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_TEAD4_HCT-116.bed

echo SRF 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46299_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46211_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/SRF.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/SRF.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/SRF.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/SRF.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/SRF.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/SRF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46299_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/SRF_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/SRF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46211_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/SRF_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46299_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_SRF_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46211_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_SRF_HCT-116.bed

echo TCF7L2 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8435_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45714_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/TCF7L2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/TCF7L2.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/TCF7L2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/TCF7L2.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/TCF7L2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/TCF7L2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8435_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/TCF7L2_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/TCF7L2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45714_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/TCF7L2_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8435_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_TCF7L2_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45714_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_TCF7L2_HCT-116.bed

echo ELF1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46266_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46202_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/ELF1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/ELF1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/ELF1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/ELF1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/ELF1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/ELF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46266_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/ELF1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/ELF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46202_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/ELF1_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46266_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_ELF1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46202_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_ELF1_HCT-116.bed

echo USF1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46279_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46205_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/USF1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/USF1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/USF1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/USF1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/USF1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/USF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46279_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/USF1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/USF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46205_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/USF1_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46279_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_USF1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46205_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_USF1_HCT-116.bed

echo SIN3A 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46277_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46208_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/SIN3A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/SIN3A.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/SIN3A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/SIN3A.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/SIN3A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/SIN3A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46277_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/SIN3A_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/SIN3A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46208_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/SIN3A_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46277_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_SIN3A_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46208_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_SIN3A_HCT-116.bed

echo SP1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46300_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46215_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/SP1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/SP1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/SP1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/SP1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/SP1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/SP1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46300_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/SP1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/SP1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46215_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/SP1_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46300_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_SP1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46215_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_SP1_HCT-116.bed

echo ZBTB33 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64738_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46200_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/ZBTB33.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/ZBTB33.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/ZBTB33.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/ZBTB33.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/ZBTB33.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/ZBTB33.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64738_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/ZBTB33_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/ZBTB33.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46200_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/ZBTB33_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64738_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_ZBTB33_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46200_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_ZBTB33_HCT-116.bed

echo ATF3 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64077_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/68464_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/ATF3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/ATF3.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/ATF3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/ATF3.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/ATF3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/ATF3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64077_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/ATF3_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/ATF3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/68464_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/ATF3_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64077_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_ATF3_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/68464_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_ATF3_HCT-116.bed

echo EGR1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46298_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46214_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/EGR1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/EGR1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/EGR1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/EGR1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/EGR1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/EGR1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46298_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/EGR1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/EGR1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46214_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/EGR1_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46298_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_EGR1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46214_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_EGR1_HCT-116.bed

echo FOSL1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63076_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46212_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/FOSL1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/FOSL1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/FOSL1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/FOSL1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/FOSL1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/FOSL1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63076_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/FOSL1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/FOSL1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46212_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/FOSL1_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63076_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_FOSL1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46212_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_FOSL1_HCT-116.bed

echo ZFX 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64866_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/81160_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/ZFX.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/ZFX.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/ZFX.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/ZFX.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/ZFX.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/ZFX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64866_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/ZFX_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/ZFX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/81160_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/ZFX_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64866_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_ZFX_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/81160_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_ZFX_HCT-116.bed

echo CBX3 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46291_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46209_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/CBX3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/CBX3.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/CBX3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HCT-116/CBX3.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/CBX3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/CBX3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46291_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/CBX3_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HCT-116/CBX3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46209_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HCT-116/CBX3_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46291_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_CBX3_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46209_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HCT-116/_CBX3_HCT-116.bed


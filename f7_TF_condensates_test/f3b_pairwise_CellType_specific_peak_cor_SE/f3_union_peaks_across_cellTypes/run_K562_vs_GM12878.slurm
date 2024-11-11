#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=8000
#SBATCH -t 4:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f3_union_peaks_across_cellTypes/slurm_K562_vs_GM12878.out

echo POLR2A 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/33859_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46137_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/POLR2A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/POLR2A.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/POLR2A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/POLR2A.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/POLR2A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/33859_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/POLR2A_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46137_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/POLR2A_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/33859_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_POLR2A_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46137_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_POLR2A_GM12878.bed

echo H3K4me3 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/72784_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42291_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K4me3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K4me3.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K4me3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K4me3.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K4me3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/72784_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K4me3_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42291_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K4me3_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/72784_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_H3K4me3_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42291_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_H3K4me3_GM12878.bed

echo H3K27ac 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/35400_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/48285_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K27ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K27ac.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K27ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K27ac.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K27ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/35400_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K27ac_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/48285_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K27ac_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/35400_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_H3K27ac_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/48285_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_H3K27ac_GM12878.bed

echo CTCF 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45899_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45631_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/CTCF.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/CTCF.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/CTCF.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/CTCF.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/CTCF.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45899_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/CTCF_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45631_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/CTCF_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45899_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_CTCF_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45631_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_CTCF_GM12878.bed

echo H3K4me2 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45382_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45214_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K4me2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K4me2.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K4me2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K4me2.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K4me2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45382_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K4me2_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45214_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K4me2_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45382_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_H3K4me2_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45214_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_H3K4me2_GM12878.bed

echo EP300 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45963_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/48272_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/EP300.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/EP300.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/EP300.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/EP300.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/EP300.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45963_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/EP300_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/48272_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/EP300_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45963_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_EP300_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/48272_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_EP300_GM12878.bed

echo H3K9ac 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45406_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/36652_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K9ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K9ac.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K9ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K9ac.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K9ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45406_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K9ac_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/36652_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K9ac_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45406_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_H3K9ac_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/36652_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_H3K9ac_GM12878.bed

echo MYC 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45894_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/666_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/MYC.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/MYC.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/MYC.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/MYC.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/MYC.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45894_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/MYC_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/666_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/MYC_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45894_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_MYC_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/666_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_MYC_GM12878.bed

echo H3K4me1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/8462_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/49297_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K4me1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K4me1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K4me1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K4me1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K4me1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K4me1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/8462_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K4me1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K4me1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/49297_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K4me1_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/8462_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_H3K4me1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/49297_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_H3K4me1_GM12878.bed

echo REST 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46286_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46113_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/REST.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/REST.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/REST.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/REST.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/REST.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/REST.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46286_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/REST_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/REST.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46113_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/REST_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46286_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_REST_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46113_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_REST_GM12878.bed

echo GABPA 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46268_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46136_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/GABPA.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/GABPA.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/GABPA.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/GABPA.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/GABPA.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/GABPA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46268_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/GABPA_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/GABPA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46136_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/GABPA_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46268_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_GABPA_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46136_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_GABPA_GM12878.bed

echo YY1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46303_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/7060_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/YY1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/YY1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/YY1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/YY1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/YY1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/YY1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46303_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/YY1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/YY1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/7060_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/YY1_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46303_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_YY1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/7060_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_YY1_GM12878.bed

echo MAX 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46292_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/661_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/MAX.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/MAX.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/MAX.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/MAX.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/MAX.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/MAX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46292_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/MAX_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/MAX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/661_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/MAX_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46292_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_MAX_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/661_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_MAX_GM12878.bed

echo RAD21 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46275_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46150_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/RAD21.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/RAD21.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/RAD21.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/RAD21.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/RAD21.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/RAD21.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46275_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/RAD21_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/RAD21.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46150_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/RAD21_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46275_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_RAD21_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46150_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_RAD21_GM12878.bed

echo TAF1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46282_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46108_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/TAF1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/TAF1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/TAF1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/TAF1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/TAF1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/TAF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46282_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/TAF1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/TAF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46108_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/TAF1_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46282_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_TAF1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46108_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_TAF1_GM12878.bed

echo H2AZ 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45416_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45221_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/H2AZ.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/H2AZ.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/H2AZ.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/H2AZ.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/H2AZ.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/H2AZ.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45416_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/H2AZ_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/H2AZ.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45221_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/H2AZ_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45416_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_H2AZ_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45221_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_H2AZ_GM12878.bed

echo SMC3 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45914_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45640_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/SMC3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/SMC3.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/SMC3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/SMC3.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/SMC3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/SMC3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45914_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/SMC3_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/SMC3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45640_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/SMC3_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45914_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_SMC3_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45640_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_SMC3_GM12878.bed

echo JUND 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45880_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/660_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/JUND.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/JUND.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/JUND.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/JUND.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/JUND.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/JUND.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45880_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/JUND_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/JUND.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/660_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/JUND_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45880_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_JUND_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/660_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_JUND_GM12878.bed

echo H3K27me3 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57482_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42282_b_sort_peaks.broadPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K27me3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K27me3.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K27me3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K27me3.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K27me3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K27me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57482_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K27me3_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K27me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42282_b_sort_peaks.broadPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/H3K27me3_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57482_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_H3K27me3_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42282_b_sort_peaks.broadPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_H3K27me3_GM12878.bed

echo ZNF143 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45946_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45629_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/ZNF143.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/ZNF143.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/ZNF143.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/ZNF143.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/ZNF143.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/ZNF143.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45946_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/ZNF143_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/ZNF143.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45629_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/ZNF143_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45946_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_ZNF143_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45629_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_ZNF143_GM12878.bed

echo NRF1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74666_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45624_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/NRF1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/NRF1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/NRF1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/NRF1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/NRF1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/NRF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74666_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/NRF1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/NRF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45624_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/NRF1_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74666_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_NRF1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45624_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_NRF1_GM12878.bed

echo SRF 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46299_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46148_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/SRF.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/SRF.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/SRF.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/SRF.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/SRF.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/SRF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46299_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/SRF_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/SRF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46148_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/SRF_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46299_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_SRF_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46148_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_SRF_GM12878.bed

echo SPI1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46272_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46119_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/SPI1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/SPI1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/SPI1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/SPI1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/SPI1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/SPI1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46272_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/SPI1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/SPI1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46119_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/SPI1_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46272_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_SPI1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46119_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_SPI1_GM12878.bed

echo ELF1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46266_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46127_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/ELF1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/ELF1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/ELF1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/ELF1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/ELF1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/ELF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46266_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/ELF1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/ELF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46127_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/ELF1_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46266_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_ELF1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46127_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_ELF1_GM12878.bed

echo CREB1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46274_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46141_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/CREB1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/CREB1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/CREB1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/CREB1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/CREB1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/CREB1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46274_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/CREB1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/CREB1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46141_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/CREB1_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46274_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_CREB1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46141_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_CREB1_GM12878.bed

echo USF1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46279_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46124_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/USF1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/USF1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/USF1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/USF1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/USF1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/USF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46279_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/USF1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/USF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46124_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/USF1_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46279_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_USF1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46124_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_USF1_GM12878.bed

echo SIN3A 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46277_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45639_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/SIN3A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/SIN3A.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/SIN3A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/SIN3A.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/SIN3A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/SIN3A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46277_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/SIN3A_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/SIN3A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45639_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/SIN3A_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46277_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_SIN3A_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45639_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_SIN3A_GM12878.bed

echo SP1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46300_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46135_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/SP1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/SP1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/SP1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/SP1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/SP1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/SP1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46300_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/SP1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/SP1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46135_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/SP1_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46300_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_SP1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46135_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_SP1_GM12878.bed

echo MXI1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45879_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45633_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/MXI1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/MXI1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/MXI1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/MXI1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/MXI1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/MXI1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45879_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/MXI1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/MXI1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45633_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/MXI1_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45879_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_MXI1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45633_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_MXI1_GM12878.bed

echo RFX5 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45893_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45649_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/RFX5.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/RFX5.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/RFX5.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/RFX5.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/RFX5.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/RFX5.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45893_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/RFX5_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/RFX5.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45649_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/RFX5_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45893_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_RFX5_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45649_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_RFX5_GM12878.bed

echo ZBTB33 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64738_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/48317_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/ZBTB33.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/ZBTB33.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/ZBTB33.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/ZBTB33.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/ZBTB33.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/ZBTB33.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64738_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/ZBTB33_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/ZBTB33.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/48317_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/ZBTB33_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64738_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_ZBTB33_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/48317_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_ZBTB33_GM12878.bed

echo ATF3 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64077_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46133_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/ATF3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/ATF3.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/ATF3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/ATF3.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/ATF3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/ATF3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64077_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/ATF3_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/ATF3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46133_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/ATF3_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64077_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_ATF3_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46133_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_ATF3_GM12878.bed

echo MAFK 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45954_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45621_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/MAFK.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/MAFK.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/MAFK.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/MAFK.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/MAFK.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/MAFK.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45954_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/MAFK_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/MAFK.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45621_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/MAFK_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45954_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_MAFK_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45621_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_MAFK_GM12878.bed

echo TBP 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45881_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45618_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/TBP.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/TBP.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/TBP.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/TBP.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/TBP.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/TBP.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45881_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/TBP_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/TBP.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45618_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/TBP_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45881_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_TBP_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45618_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_TBP_GM12878.bed

echo MAZ 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45924_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45661_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/MAZ.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/MAZ.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/MAZ.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/MAZ.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/MAZ.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/MAZ.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45924_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/MAZ_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/MAZ.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45661_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/MAZ_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45924_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_MAZ_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45661_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_MAZ_GM12878.bed

echo EGR1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46298_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46126_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/EGR1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/EGR1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/EGR1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/EGR1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/EGR1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/EGR1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46298_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/EGR1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/EGR1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46126_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/EGR1_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46298_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_EGR1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46126_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_EGR1_GM12878.bed

echo CEBPA 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8456_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/48277_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/CEBPA.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/CEBPA.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/CEBPA.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/CEBPA.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/CEBPA.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/CEBPA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8456_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/CEBPA_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/CEBPA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/48277_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/CEBPA_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8456_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_CEBPA_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/48277_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_CEBPA_GM12878.bed

echo FOS 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1990_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45606_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/FOS.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/FOS.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/FOS.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/FOS.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/FOS.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/FOS.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1990_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/FOS_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/FOS.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45606_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/FOS_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1990_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_FOS_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45606_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_FOS_GM12878.bed

echo BHLHE40 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45916_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45632_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/BHLHE40.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/BHLHE40.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/BHLHE40.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/BHLHE40.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/BHLHE40.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/BHLHE40.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45916_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/BHLHE40_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/BHLHE40.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45632_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/BHLHE40_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45916_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_BHLHE40_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45632_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_BHLHE40_GM12878.bed

echo NR2C2 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64395_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/669_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/NR2C2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/NR2C2.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/NR2C2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/NR2C2.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/NR2C2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/NR2C2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64395_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/NR2C2_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/NR2C2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/669_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/NR2C2_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64395_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_NR2C2_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/669_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_NR2C2_GM12878.bed

echo USF2 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63383_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45611_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/USF2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/USF2.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/USF2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/USF2.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/USF2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/USF2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63383_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/USF2_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/USF2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45611_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/USF2_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63383_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_USF2_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45611_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_USF2_GM12878.bed

echo CHD2 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45905_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45638_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/CHD2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/CHD2.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/CHD2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/CHD2.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/CHD2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/CHD2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45905_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/CHD2_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/CHD2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45638_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/CHD2_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45905_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_CHD2_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45638_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_CHD2_GM12878.bed

echo RCOR1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45903_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45652_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/RCOR1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/RCOR1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/RCOR1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/RCOR1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/RCOR1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/RCOR1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45903_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/RCOR1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/RCOR1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45652_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/RCOR1_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45903_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_RCOR1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45652_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_RCOR1_GM12878.bed

echo E2F4 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45951_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45617_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/E2F4.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/E2F4.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/E2F4.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/E2F4.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/E2F4.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/E2F4.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45951_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/E2F4_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/E2F4.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45617_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/E2F4_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45951_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_E2F4_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45617_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_E2F4_GM12878.bed

echo ETS1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46278_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46117_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/ETS1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/ETS1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/ETS1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/ETS1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/ETS1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/ETS1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46278_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/ETS1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/ETS1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46117_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/ETS1_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46278_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_ETS1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46117_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_ETS1_GM12878.bed

echo ZBTB40 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63068_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62937_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/ZBTB40.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/ZBTB40.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/ZBTB40.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/ZBTB40.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/ZBTB40.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/ZBTB40.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63068_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/ZBTB40_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/ZBTB40.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62937_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/ZBTB40_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63068_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_ZBTB40_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62937_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_ZBTB40_GM12878.bed

echo SIX5 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46302_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46142_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/SIX5.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/SIX5.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/SIX5.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/SIX5.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/SIX5.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/SIX5.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46302_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/SIX5_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/SIX5.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46142_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/SIX5_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46302_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_SIX5_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46142_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_SIX5_GM12878.bed

echo NFYA 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2094_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45626_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/NFYA.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/NFYA.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/NFYA.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/NFYA.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/NFYA.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/NFYA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2094_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/NFYA_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/NFYA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45626_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/NFYA_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2094_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_NFYA_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45626_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_NFYA_GM12878.bed

echo ELK1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63321_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45636_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/ELK1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/ELK1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/ELK1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/ELK1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/ELK1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/ELK1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63321_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/ELK1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/ELK1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45636_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/ELK1_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63321_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_ELK1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45636_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_ELK1_GM12878.bed

echo NFYB 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45936_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45625_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/NFYB.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/NFYB.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/NFYB.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/NFYB.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/NFYB.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/NFYB.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45936_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/NFYB_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/NFYB.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45625_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/NFYB_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45936_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_NFYB_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45625_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_NFYB_GM12878.bed

echo CUX1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45934_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45615_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/CUX1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/CUX1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/CUX1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/CUX1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/CUX1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/CUX1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45934_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/CUX1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/CUX1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45615_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/CUX1_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45934_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_CUX1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45615_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_CUX1_GM12878.bed

echo JUNB 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/7487_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64781_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/JUNB.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/JUNB.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/JUNB.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/JUNB.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/JUNB.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/JUNB.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/7487_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/JUNB_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/JUNB.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64781_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/JUNB_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/7487_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_JUNB_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64781_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_JUNB_GM12878.bed

echo MEF2A 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46283_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46123_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/MEF2A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/MEF2A.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/MEF2A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/MEF2A.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/MEF2A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/MEF2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46283_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/MEF2A_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/MEF2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46123_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/MEF2A_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46283_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_MEF2A_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46123_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_MEF2A_GM12878.bed

echo ZMYM3 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62709_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74652_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/ZMYM3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/ZMYM3.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/ZMYM3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/ZMYM3.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/ZMYM3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/ZMYM3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62709_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/ZMYM3_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/ZMYM3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74652_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/ZMYM3_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62709_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_ZMYM3_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74652_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_ZMYM3_GM12878.bed

echo ZNF384 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45922_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45634_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/ZNF384.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/ZNF384.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/ZNF384.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/ZNF384.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/ZNF384.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/ZNF384.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45922_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/ZNF384_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/ZNF384.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45634_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/ZNF384_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45922_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_ZNF384_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45634_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_ZNF384_GM12878.bed

echo PML 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46301_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46143_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/PML.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/PML.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/PML.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/PML.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/PML.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/PML.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46301_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/PML_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/PML.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46143_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/PML_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46301_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_PML_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46143_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_PML_GM12878.bed

echo STAT5A 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46280_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46146_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/STAT5A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/STAT5A.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/STAT5A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/STAT5A.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/STAT5A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/STAT5A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46280_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/STAT5A_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/STAT5A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46146_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/STAT5A_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46280_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_STAT5A_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46146_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_STAT5A_GM12878.bed

echo SUZ12 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63519_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64382_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/SUZ12.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/SUZ12.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/SUZ12.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/SUZ12.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/SUZ12.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/SUZ12.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63519_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/SUZ12_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/SUZ12.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64382_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/SUZ12_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63519_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_SUZ12_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64382_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_SUZ12_GM12878.bed

echo TBL1XR1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45876_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45660_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/TBL1XR1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/TBL1XR1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/TBL1XR1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/TBL1XR1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/TBL1XR1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/TBL1XR1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45876_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/TBL1XR1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/TBL1XR1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45660_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/TBL1XR1_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45876_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_TBL1XR1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45660_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_TBL1XR1_GM12878.bed

echo BCLAF1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46288_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/6958_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/BCLAF1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/BCLAF1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/BCLAF1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/BCLAF1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/BCLAF1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/BCLAF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46288_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/BCLAF1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/BCLAF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/6958_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/BCLAF1_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46288_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_BCLAF1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/6958_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_BCLAF1_GM12878.bed

echo IKZF1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63463_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64735_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/IKZF1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/IKZF1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/IKZF1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/IKZF1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/IKZF1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/IKZF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63463_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/IKZF1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/IKZF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64735_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/IKZF1_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63463_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_IKZF1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64735_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_IKZF1_GM12878.bed

echo MLLT1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62738_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63876_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/MLLT1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/MLLT1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/MLLT1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/MLLT1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/MLLT1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/MLLT1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62738_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/MLLT1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/MLLT1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63876_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/MLLT1_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62738_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_MLLT1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63876_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_MLLT1_GM12878.bed

echo MTA2 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62756_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63213_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/MTA2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/MTA2.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/MTA2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/MTA2.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/MTA2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/MTA2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62756_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/MTA2_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/MTA2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63213_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/MTA2_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62756_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_MTA2_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63213_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_MTA2_GM12878.bed

echo NFE2 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63878_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45608_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/NFE2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/NFE2.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/NFE2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_GM12878/NFE2.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_GM12878/NFE2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/NFE2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63878_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/NFE2_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_GM12878/NFE2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45608_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_GM12878/NFE2_merge_intersect_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63878_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_NFE2_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45608_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_GM12878/_NFE2_GM12878.bed


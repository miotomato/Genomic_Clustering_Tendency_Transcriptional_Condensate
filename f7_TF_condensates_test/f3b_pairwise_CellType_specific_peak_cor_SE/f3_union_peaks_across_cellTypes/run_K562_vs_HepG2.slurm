#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=8000
#SBATCH -t 4:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f3_union_peaks_across_cellTypes/slurm_K562_vs_HepG2.out

echo POLR2A 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/33859_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46256_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/POLR2A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/POLR2A.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/POLR2A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/POLR2A.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/POLR2A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/33859_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/POLR2A_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46256_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/POLR2A_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/33859_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_POLR2A_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46256_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_POLR2A_HepG2.bed

echo H3K4me3 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/72784_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45306_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K4me3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K4me3.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K4me3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K4me3.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K4me3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/72784_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K4me3_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45306_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K4me3_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/72784_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_H3K4me3_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45306_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_H3K4me3_HepG2.bed

echo H3K27ac 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/35400_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45299_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K27ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K27ac.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K27ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K27ac.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K27ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/35400_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K27ac_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45299_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K27ac_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/35400_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_H3K27ac_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45299_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_H3K27ac_HepG2.bed

echo CTCF 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45899_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1380_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/CTCF.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/CTCF.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/CTCF.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/CTCF.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/CTCF.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45899_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/CTCF_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1380_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/CTCF_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45899_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_CTCF_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1380_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_CTCF_HepG2.bed

echo H3K4me2 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45382_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45307_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K4me2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K4me2.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K4me2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K4me2.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K4me2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45382_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K4me2_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45307_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K4me2_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45382_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_H3K4me2_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45307_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_H3K4me2_HepG2.bed

echo EP300 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45963_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46225_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/EP300.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/EP300.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/EP300.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/EP300.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/EP300.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45963_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/EP300_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46225_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/EP300_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45963_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_EP300_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46225_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_EP300_HepG2.bed

echo H3K9ac 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45406_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45305_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K9ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K9ac.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K9ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K9ac.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K9ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45406_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K9ac_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45305_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K9ac_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45406_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_H3K9ac_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45305_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_H3K9ac_HepG2.bed

echo MYC 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45894_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45827_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/MYC.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/MYC.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/MYC.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/MYC.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/MYC.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45894_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/MYC_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45827_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/MYC_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45894_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_MYC_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45827_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_MYC_HepG2.bed

echo H3K4me1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/8462_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/51255_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K4me1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K4me1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K4me1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K4me1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K4me1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K4me1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/8462_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K4me1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K4me1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/51255_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K4me1_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/8462_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_H3K4me1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/51255_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_H3K4me1_HepG2.bed

echo REST 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46286_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46257_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/REST.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/REST.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/REST.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/REST.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/REST.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/REST.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46286_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/REST_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/REST.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46257_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/REST_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46286_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_REST_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46257_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_REST_HepG2.bed

echo GABPA 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46268_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56045_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/GABPA.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/GABPA.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/GABPA.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/GABPA.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/GABPA.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/GABPA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46268_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/GABPA_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/GABPA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56045_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/GABPA_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46268_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_GABPA_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56045_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_GABPA_HepG2.bed

echo YY1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46303_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46226_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/YY1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/YY1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/YY1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/YY1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/YY1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/YY1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46303_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/YY1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/YY1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46226_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/YY1_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46303_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_YY1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46226_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_YY1_HepG2.bed

echo MAX 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46292_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46244_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/MAX.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/MAX.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/MAX.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/MAX.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/MAX.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/MAX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46292_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/MAX_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/MAX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46244_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/MAX_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46292_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_MAX_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46244_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_MAX_HepG2.bed

echo RAD21 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46275_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46252_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/RAD21.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/RAD21.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/RAD21.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/RAD21.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/RAD21.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/RAD21.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46275_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/RAD21_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/RAD21.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46252_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/RAD21_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46275_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_RAD21_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46252_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_RAD21_HepG2.bed

echo TAF1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46282_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46258_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/TAF1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/TAF1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/TAF1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/TAF1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/TAF1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/TAF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46282_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/TAF1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/TAF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46258_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/TAF1_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46282_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_TAF1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46258_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_TAF1_HepG2.bed

echo CEBPB 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45923_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45814_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/CEBPB.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/CEBPB.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/CEBPB.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/CEBPB.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/CEBPB.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/CEBPB.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45923_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/CEBPB_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/CEBPB.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45814_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/CEBPB_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45923_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_CEBPB_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45814_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_CEBPB_HepG2.bed

echo SMC3 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45914_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45816_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/SMC3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/SMC3.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/SMC3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/SMC3.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/SMC3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/SMC3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45914_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/SMC3_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/SMC3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45816_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/SMC3_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45914_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_SMC3_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45816_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_SMC3_HepG2.bed

echo JUND 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45880_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45828_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/JUND.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/JUND.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/JUND.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/JUND.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/JUND.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/JUND.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45880_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/JUND_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/JUND.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45828_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/JUND_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45880_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_JUND_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45828_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_JUND_HepG2.bed

echo TEAD4 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46284_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46238_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/TEAD4.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/TEAD4.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/TEAD4.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/TEAD4.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/TEAD4.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/TEAD4.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46284_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/TEAD4_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/TEAD4.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46238_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/TEAD4_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46284_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_TEAD4_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46238_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_TEAD4_HepG2.bed

echo H3K9me3 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/56869_b_sort_peaks.broadPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45300_b_sort_peaks.broadPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K9me3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K9me3.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K9me3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K9me3.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K9me3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K9me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/56869_b_sort_peaks.broadPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K9me3_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K9me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45300_b_sort_peaks.broadPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/H3K9me3_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/56869_b_sort_peaks.broadPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_H3K9me3_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45300_b_sort_peaks.broadPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_H3K9me3_HepG2.bed

echo NRF1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74666_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64667_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/NRF1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/NRF1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/NRF1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/NRF1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/NRF1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/NRF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74666_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/NRF1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/NRF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64667_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/NRF1_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74666_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_NRF1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64667_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_NRF1_HepG2.bed

echo SRF 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46299_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46223_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/SRF.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/SRF.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/SRF.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/SRF.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/SRF.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/SRF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46299_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/SRF_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/SRF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46223_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/SRF_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46299_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_SRF_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46223_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_SRF_HepG2.bed

echo TCF7L2 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8435_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45825_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/TCF7L2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/TCF7L2.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/TCF7L2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/TCF7L2.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/TCF7L2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/TCF7L2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8435_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/TCF7L2_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/TCF7L2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45825_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/TCF7L2_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8435_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_TCF7L2_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45825_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_TCF7L2_HepG2.bed

echo ELF1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46266_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46254_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/ELF1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/ELF1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/ELF1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/ELF1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/ELF1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/ELF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46266_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/ELF1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/ELF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46254_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/ELF1_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46266_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_ELF1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46254_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_ELF1_HepG2.bed

echo CREB1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46274_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56047_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/CREB1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/CREB1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/CREB1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/CREB1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/CREB1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/CREB1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46274_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/CREB1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/CREB1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56047_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/CREB1_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46274_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_CREB1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56047_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_CREB1_HepG2.bed

echo JUN 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2078_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45830_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/JUN.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/JUN.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/JUN.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/JUN.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/JUN.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/JUN.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2078_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/JUN_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/JUN.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45830_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/JUN_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2078_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_JUN_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45830_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_JUN_HepG2.bed

echo MED1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74667_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/71877_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/MED1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/MED1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/MED1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/MED1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/MED1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/MED1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74667_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/MED1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/MED1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/71877_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/MED1_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74667_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_MED1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/71877_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_MED1_HepG2.bed

echo USF1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46279_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46248_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/USF1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/USF1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/USF1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/USF1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/USF1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/USF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46279_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/USF1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/USF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46248_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/USF1_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46279_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_USF1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46248_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_USF1_HepG2.bed

echo SIN3A 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46277_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46236_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/SIN3A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/SIN3A.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/SIN3A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/SIN3A.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/SIN3A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/SIN3A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46277_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/SIN3A_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/SIN3A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46236_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/SIN3A_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46277_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_SIN3A_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46236_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_SIN3A_HepG2.bed

echo SP1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46300_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46234_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/SP1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/SP1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/SP1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/SP1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/SP1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/SP1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46300_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/SP1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/SP1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46234_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/SP1_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46300_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_SP1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46234_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_SP1_HepG2.bed

echo MXI1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45879_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45810_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/MXI1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/MXI1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/MXI1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/MXI1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/MXI1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/MXI1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45879_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/MXI1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/MXI1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45810_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/MXI1_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45879_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_MXI1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45810_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_MXI1_HepG2.bed

echo RFX5 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45893_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45823_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/RFX5.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/RFX5.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/RFX5.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/RFX5.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/RFX5.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/RFX5.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45893_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/RFX5_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/RFX5.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45823_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/RFX5_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45893_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_RFX5_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45823_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_RFX5_HepG2.bed

echo ZBTB33 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64738_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46229_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/ZBTB33.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/ZBTB33.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/ZBTB33.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/ZBTB33.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/ZBTB33.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/ZBTB33.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64738_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/ZBTB33_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/ZBTB33.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46229_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/ZBTB33_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64738_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_ZBTB33_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46229_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_ZBTB33_HepG2.bed

echo HDAC2 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46304_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46235_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/HDAC2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/HDAC2.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/HDAC2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/HDAC2.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/HDAC2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/HDAC2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46304_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/HDAC2_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/HDAC2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46235_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/HDAC2_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46304_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_HDAC2_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46235_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_HDAC2_HepG2.bed

echo ATF3 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64077_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/7301_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/ATF3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/ATF3.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/ATF3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/ATF3.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/ATF3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/ATF3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64077_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/ATF3_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/ATF3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/7301_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/ATF3_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64077_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_ATF3_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/7301_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_ATF3_HepG2.bed

echo MAFK 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45954_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45815_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/MAFK.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/MAFK.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/MAFK.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/MAFK.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/MAFK.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/MAFK.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45954_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/MAFK_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/MAFK.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45815_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/MAFK_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45954_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_MAFK_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45815_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_MAFK_HepG2.bed

echo TBP 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45881_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45791_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/TBP.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/TBP.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/TBP.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/TBP.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/TBP.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/TBP.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45881_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/TBP_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/TBP.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45791_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/TBP_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45881_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_TBP_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45791_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_TBP_HepG2.bed

echo MAZ 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45924_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45788_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/MAZ.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/MAZ.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/MAZ.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/MAZ.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/MAZ.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/MAZ.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45924_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/MAZ_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/MAZ.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45788_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/MAZ_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45924_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_MAZ_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45788_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_MAZ_HepG2.bed

echo CEBPA 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8456_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1374_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/CEBPA.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/CEBPA.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/CEBPA.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/CEBPA.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/CEBPA.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/CEBPA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8456_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/CEBPA_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/CEBPA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1374_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/CEBPA_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8456_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_CEBPA_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1374_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_CEBPA_HepG2.bed

echo BHLHE40 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45916_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45807_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/BHLHE40.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/BHLHE40.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/BHLHE40.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/BHLHE40.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/BHLHE40.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/BHLHE40.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45916_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/BHLHE40_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/BHLHE40.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45807_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/BHLHE40_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45916_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_BHLHE40_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45807_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_BHLHE40_HepG2.bed

echo NR2C2 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64395_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45800_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/NR2C2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/NR2C2.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/NR2C2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/NR2C2.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/NR2C2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/NR2C2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64395_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/NR2C2_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/NR2C2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45800_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/NR2C2_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64395_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_NR2C2_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45800_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_NR2C2_HepG2.bed

echo USF2 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63383_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74655_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/USF2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/USF2.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/USF2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/USF2.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/USF2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/USF2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63383_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/USF2_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/USF2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74655_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/USF2_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63383_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_USF2_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74655_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_USF2_HepG2.bed

echo CHD2 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45905_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45802_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/CHD2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/CHD2.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/CHD2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/CHD2.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/CHD2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/CHD2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45905_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/CHD2_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/CHD2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45802_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/CHD2_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45905_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_CHD2_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45802_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_CHD2_HepG2.bed

echo RCOR1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45903_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45798_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/RCOR1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/RCOR1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/RCOR1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/RCOR1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/RCOR1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/RCOR1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45903_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/RCOR1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/RCOR1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45798_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/RCOR1_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45903_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_RCOR1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45798_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_RCOR1_HepG2.bed

echo ZBTB40 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63068_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63800_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/ZBTB40.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/ZBTB40.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/ZBTB40.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/ZBTB40.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/ZBTB40.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/ZBTB40.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63068_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/ZBTB40_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/ZBTB40.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63800_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/ZBTB40_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63068_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_ZBTB40_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63800_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_ZBTB40_HepG2.bed

echo SP2 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46287_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46246_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/SP2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/SP2.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/SP2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/SP2.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/SP2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/SP2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46287_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/SP2_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/SP2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46246_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/SP2_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46287_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_SP2_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46246_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_SP2_HepG2.bed

echo NR2F2 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46265_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46249_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/NR2F2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/NR2F2.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/NR2F2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/NR2F2.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/NR2F2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/NR2F2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46265_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/NR2F2_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/NR2F2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46249_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/NR2F2_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46265_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_NR2F2_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46249_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_NR2F2_HepG2.bed

echo CEBPD 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46294_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46230_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/CEBPD.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/CEBPD.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/CEBPD.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/CEBPD.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/CEBPD.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/CEBPD.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46294_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/CEBPD_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/CEBPD.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46230_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/CEBPD_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46294_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_CEBPD_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46230_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_CEBPD_HepG2.bed

echo HES1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62672_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74644_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/HES1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/HES1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/HES1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/HES1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/HES1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/HES1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62672_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/HES1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/HES1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74644_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/HES1_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62672_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_HES1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74644_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_HES1_HepG2.bed

echo ZMYM3 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62709_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64657_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/ZMYM3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/ZMYM3.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/ZMYM3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/ZMYM3.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/ZMYM3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/ZMYM3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62709_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/ZMYM3_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/ZMYM3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64657_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/ZMYM3_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62709_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_ZMYM3_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64657_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_ZMYM3_HepG2.bed

echo ZBTB7A 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46295_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46224_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/ZBTB7A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/ZBTB7A.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/ZBTB7A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/ZBTB7A.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/ZBTB7A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/ZBTB7A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46295_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/ZBTB7A_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/ZBTB7A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46224_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/ZBTB7A_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46295_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_ZBTB7A_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46224_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_ZBTB7A_HepG2.bed

echo FOXK2 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63743_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62913_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/FOXK2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/FOXK2.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/FOXK2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/FOXK2.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/FOXK2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/FOXK2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63743_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/FOXK2_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/FOXK2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62913_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/FOXK2_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63743_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_FOXK2_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62913_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_FOXK2_HepG2.bed

echo FOXO3 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74680_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74682_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/FOXO3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/FOXO3.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/FOXO3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/FOXO3.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/FOXO3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/FOXO3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74680_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/FOXO3_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/FOXO3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74682_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/FOXO3_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74680_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_FOXO3_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74682_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_FOXO3_HepG2.bed

echo TCF7 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64697_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63592_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/TCF7.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/TCF7.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/TCF7.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/TCF7.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/TCF7.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/TCF7.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64697_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/TCF7_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/TCF7.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63592_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/TCF7_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64697_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_TCF7_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63592_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_TCF7_HepG2.bed

echo CBX1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64916_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64582_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/CBX1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/CBX1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/CBX1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/CBX1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/CBX1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/CBX1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64916_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/CBX1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/CBX1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64582_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/CBX1_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64916_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_CBX1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64582_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_CBX1_HepG2.bed

echo ATF1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45871_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56053_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/ATF1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/ATF1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/ATF1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/ATF1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/ATF1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/ATF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45871_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/ATF1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/ATF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56053_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/ATF1_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45871_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_ATF1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56053_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_ATF1_HepG2.bed

echo ARID3A 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45897_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45806_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/ARID3A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/ARID3A.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/ARID3A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/ARID3A.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/ARID3A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/ARID3A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45897_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/ARID3A_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/ARID3A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45806_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/ARID3A_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45897_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_ARID3A_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45806_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_ARID3A_HepG2.bed

echo MAFF 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45950_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45801_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/MAFF.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/MAFF.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/MAFF.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/MAFF.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/MAFF.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/MAFF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45950_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/MAFF_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/MAFF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45801_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/MAFF_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45950_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_MAFF_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45801_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_MAFF_HepG2.bed

echo RAD51 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63793_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62652_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/RAD51.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/RAD51.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/RAD51.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/RAD51.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/RAD51.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/RAD51.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63793_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/RAD51_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/RAD51.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62652_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/RAD51_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63793_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_RAD51_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62652_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_RAD51_HepG2.bed

echo NFIA 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/77263_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/77267_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/NFIA.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/NFIA.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/NFIA.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/NFIA.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/NFIA.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/NFIA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/77263_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/NFIA_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/NFIA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/77267_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/NFIA_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/77263_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_NFIA_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/77267_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_NFIA_HepG2.bed

echo PBX2 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63133_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64661_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/PBX2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/PBX2.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/PBX2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/PBX2.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/PBX2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/PBX2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63133_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/PBX2_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/PBX2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64661_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/PBX2_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63133_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_PBX2_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64661_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_PBX2_HepG2.bed

echo NR2F6 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64300_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63772_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/NR2F6.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/NR2F6.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/NR2F6.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/NR2F6.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/NR2F6.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/NR2F6.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64300_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/NR2F6_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/NR2F6.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63772_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/NR2F6_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64300_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_NR2F6_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63772_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_NR2F6_HepG2.bed

echo HEY1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/7527_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/7341_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/HEY1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/HEY1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/HEY1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HepG2/HEY1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HepG2/HEY1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/HEY1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/7527_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/HEY1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HepG2/HEY1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/7341_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HepG2/HEY1_merge_intersect_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/7527_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_HEY1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/7341_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HepG2/_HEY1_HepG2.bed


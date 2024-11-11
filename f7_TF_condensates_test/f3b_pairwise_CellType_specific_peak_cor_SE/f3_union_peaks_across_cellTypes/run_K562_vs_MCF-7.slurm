#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=8000
#SBATCH -t 4:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f3_union_peaks_across_cellTypes/slurm_K562_vs_MCF-7.out

echo POLR2A 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/33859_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45985_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/POLR2A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/POLR2A.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/POLR2A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/POLR2A.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/POLR2A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/33859_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/POLR2A_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45985_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/POLR2A_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/33859_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_POLR2A_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45985_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_POLR2A_MCF-7.bed

echo H3K4me3 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/72784_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/47270_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K4me3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K4me3.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K4me3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K4me3.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K4me3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/72784_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K4me3_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/47270_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K4me3_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/72784_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_H3K4me3_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/47270_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_H3K4me3_MCF-7.bed

echo H3K27ac 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/35400_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/33141_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K27ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K27ac.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K27ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K27ac.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K27ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/35400_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K27ac_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/33141_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K27ac_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/35400_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_H3K27ac_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/33141_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_H3K27ac_MCF-7.bed

echo CTCF 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45899_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/83736_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/CTCF.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/CTCF.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/CTCF.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/CTCF.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/CTCF.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45899_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/CTCF_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/83736_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/CTCF_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45899_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_CTCF_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/83736_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_CTCF_MCF-7.bed

echo H3K4me2 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45382_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/75112_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K4me2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K4me2.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K4me2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K4me2.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K4me2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45382_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K4me2_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/75112_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K4me2_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45382_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_H3K4me2_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/75112_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_H3K4me2_MCF-7.bed

echo EP300 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45963_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/33155_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/EP300.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/EP300.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/EP300.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/EP300.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/EP300.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45963_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/EP300_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/33155_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/EP300_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45963_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_EP300_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/33155_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_EP300_MCF-7.bed

echo H3K9ac 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45406_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/86292_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K9ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K9ac.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K9ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K9ac.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K9ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45406_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K9ac_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/86292_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K9ac_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45406_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_H3K9ac_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/86292_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_H3K9ac_MCF-7.bed

echo MYC 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45894_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/38193_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/MYC.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/MYC.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/MYC.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/MYC.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/MYC.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45894_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/MYC_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/38193_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/MYC_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45894_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_MYC_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/38193_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_MYC_MCF-7.bed

echo H3K4me1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/8462_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/33151_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K4me1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K4me1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K4me1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K4me1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K4me1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K4me1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/8462_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K4me1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K4me1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/33151_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K4me1_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/8462_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_H3K4me1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/33151_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_H3K4me1_MCF-7.bed

echo REST 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46286_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46313_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/REST.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/REST.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/REST.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/REST.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/REST.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/REST.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46286_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/REST_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/REST.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46313_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/REST_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46286_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_REST_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46313_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_REST_MCF-7.bed

echo GABPA 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46268_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46314_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/GABPA.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/GABPA.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/GABPA.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/GABPA.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/GABPA.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/GABPA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46268_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/GABPA_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/GABPA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46314_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/GABPA_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46268_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_GABPA_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46314_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_GABPA_MCF-7.bed

echo MAX 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46292_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46324_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/MAX.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/MAX.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/MAX.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/MAX.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/MAX.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/MAX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46292_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/MAX_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/MAX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46324_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/MAX_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46292_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_MAX_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46324_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_MAX_MCF-7.bed

echo RAD21 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46275_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2387_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/RAD21.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/RAD21.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/RAD21.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/RAD21.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/RAD21.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/RAD21.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46275_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/RAD21_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/RAD21.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2387_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/RAD21_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46275_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_RAD21_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2387_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_RAD21_MCF-7.bed

echo TAF1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46282_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46320_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/TAF1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/TAF1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/TAF1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/TAF1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/TAF1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/TAF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46282_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/TAF1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/TAF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46320_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/TAF1_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46282_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_TAF1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46320_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_TAF1_MCF-7.bed

echo RUNX1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2121_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/70984_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/RUNX1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/RUNX1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/RUNX1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/RUNX1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/RUNX1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/RUNX1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2121_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/RUNX1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/RUNX1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/70984_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/RUNX1_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2121_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_RUNX1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/70984_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_RUNX1_MCF-7.bed

echo H2AZ 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45416_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/76117_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H2AZ.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H2AZ.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H2AZ.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H2AZ.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H2AZ.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H2AZ.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45416_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H2AZ_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H2AZ.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/76117_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H2AZ_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45416_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_H2AZ_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/76117_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_H2AZ_MCF-7.bed

echo CEBPB 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45923_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46316_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/CEBPB.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/CEBPB.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/CEBPB.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/CEBPB.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/CEBPB.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/CEBPB.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45923_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/CEBPB_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/CEBPB.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46316_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/CEBPB_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45923_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_CEBPB_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46316_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_CEBPB_MCF-7.bed

echo JUND 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45880_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46308_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/JUND.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/JUND.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/JUND.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/JUND.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/JUND.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/JUND.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45880_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/JUND_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/JUND.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46308_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/JUND_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45880_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_JUND_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46308_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_JUND_MCF-7.bed

echo ZNF143 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45946_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/71042_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ZNF143.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ZNF143.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ZNF143.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ZNF143.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ZNF143.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ZNF143.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45946_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ZNF143_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ZNF143.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/71042_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ZNF143_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45946_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_ZNF143_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/71042_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_ZNF143_MCF-7.bed

echo TEAD4 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46284_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46315_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/TEAD4.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/TEAD4.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/TEAD4.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/TEAD4.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/TEAD4.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/TEAD4.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46284_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/TEAD4_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/TEAD4.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46315_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/TEAD4_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46284_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_TEAD4_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46315_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_TEAD4_MCF-7.bed

echo H3K9me3 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/56869_b_sort_peaks.broadPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/49474_b_sort_peaks.broadPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K9me3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K9me3.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K9me3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K9me3.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K9me3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K9me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/56869_b_sort_peaks.broadPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K9me3_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K9me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/49474_b_sort_peaks.broadPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/H3K9me3_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/56869_b_sort_peaks.broadPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_H3K9me3_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/49474_b_sort_peaks.broadPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_H3K9me3_MCF-7.bed

echo NRF1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74666_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74638_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/NRF1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/NRF1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/NRF1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/NRF1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/NRF1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/NRF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74666_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/NRF1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/NRF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74638_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/NRF1_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74666_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_NRF1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74638_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_NRF1_MCF-7.bed

echo SRF 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46299_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46312_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/SRF.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/SRF.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/SRF.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/SRF.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/SRF.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/SRF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46299_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/SRF_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/SRF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46312_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/SRF_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46299_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_SRF_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46312_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_SRF_MCF-7.bed

echo TCF7L2 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8435_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45979_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/TCF7L2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/TCF7L2.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/TCF7L2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/TCF7L2.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/TCF7L2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/TCF7L2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8435_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/TCF7L2_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/TCF7L2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45979_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/TCF7L2_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8435_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_TCF7L2_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45979_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_TCF7L2_MCF-7.bed

echo ELF1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46266_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46323_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ELF1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ELF1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ELF1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ELF1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ELF1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ELF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46266_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ELF1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ELF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46323_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ELF1_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46266_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_ELF1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46323_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_ELF1_MCF-7.bed

echo MED1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74667_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/49652_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/MED1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/MED1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/MED1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/MED1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/MED1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/MED1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74667_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/MED1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/MED1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/49652_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/MED1_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74667_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_MED1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/49652_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_MED1_MCF-7.bed

echo ERG 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/37639_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/37633_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ERG.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ERG.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ERG.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ERG.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ERG.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ERG.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/37639_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ERG_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ERG.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/37633_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ERG_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/37639_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_ERG_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/37633_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_ERG_MCF-7.bed

echo SIN3A 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46277_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46327_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/SIN3A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/SIN3A.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/SIN3A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/SIN3A.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/SIN3A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/SIN3A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46277_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/SIN3A_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/SIN3A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46327_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/SIN3A_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46277_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_SIN3A_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46327_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_SIN3A_MCF-7.bed

echo HDAC2 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46304_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46326_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/HDAC2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/HDAC2.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/HDAC2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/HDAC2.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/HDAC2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/HDAC2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46304_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/HDAC2_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/HDAC2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46326_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/HDAC2_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46304_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_HDAC2_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46326_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_HDAC2_MCF-7.bed

echo EGR1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46298_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46319_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/EGR1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/EGR1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/EGR1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/EGR1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/EGR1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/EGR1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46298_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/EGR1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/EGR1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46319_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/EGR1_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46298_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_EGR1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46319_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_EGR1_MCF-7.bed

echo HCFC1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45890_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64266_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/HCFC1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/HCFC1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/HCFC1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/HCFC1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/HCFC1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/HCFC1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45890_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/HCFC1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/HCFC1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64266_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/HCFC1_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45890_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_HCFC1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64266_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_HCFC1_MCF-7.bed

echo FOS 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1990_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2318_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/FOS.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/FOS.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/FOS.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/FOS.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/FOS.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/FOS.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1990_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/FOS_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/FOS.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2318_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/FOS_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1990_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_FOS_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2318_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_FOS_MCF-7.bed

echo E2F4 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45951_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/36807_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/E2F4.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/E2F4.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/E2F4.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/E2F4.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/E2F4.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/E2F4.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45951_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/E2F4_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/E2F4.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/36807_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/E2F4_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45951_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_E2F4_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/36807_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_E2F4_MCF-7.bed

echo KDM5B 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8791_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/47267_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/KDM5B.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/KDM5B.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/KDM5B.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/KDM5B.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/KDM5B.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/KDM5B.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8791_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/KDM5B_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/KDM5B.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/47267_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/KDM5B_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8791_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_KDM5B_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/47267_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_KDM5B_MCF-7.bed

echo ZBTB40 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63068_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63278_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ZBTB40.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ZBTB40.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ZBTB40.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ZBTB40.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ZBTB40.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ZBTB40.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63068_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ZBTB40_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ZBTB40.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63278_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ZBTB40_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63068_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_ZBTB40_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63278_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_ZBTB40_MCF-7.bed

echo ELK1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63321_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63433_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ELK1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ELK1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ELK1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ELK1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ELK1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ELK1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63321_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ELK1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ELK1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63433_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ELK1_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63321_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_ELK1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63433_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_ELK1_MCF-7.bed

echo NR2F2 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46265_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46309_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/NR2F2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/NR2F2.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/NR2F2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/NR2F2.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/NR2F2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/NR2F2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46265_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/NR2F2_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/NR2F2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46309_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/NR2F2_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46265_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_NR2F2_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46309_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_NR2F2_MCF-7.bed

echo MBD2 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63028_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/48825_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/MBD2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/MBD2.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/MBD2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/MBD2.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/MBD2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/MBD2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63028_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/MBD2_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/MBD2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/48825_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/MBD2_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63028_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_MBD2_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/48825_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_MBD2_MCF-7.bed

echo CUX1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45934_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62548_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/CUX1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/CUX1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/CUX1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/CUX1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/CUX1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/CUX1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45934_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/CUX1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/CUX1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62548_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/CUX1_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45934_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_CUX1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62548_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_CUX1_MCF-7.bed

echo RFX1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62594_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62623_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/RFX1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/RFX1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/RFX1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/RFX1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/RFX1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/RFX1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62594_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/RFX1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/RFX1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62623_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/RFX1_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62594_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_RFX1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62623_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_RFX1_MCF-7.bed

echo ZFX 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64866_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/84928_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ZFX.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ZFX.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ZFX.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ZFX.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ZFX.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ZFX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64866_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ZFX_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ZFX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/84928_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ZFX_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64866_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_ZFX_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/84928_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_ZFX_MCF-7.bed

echo ZBTB11 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64291_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62870_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ZBTB11.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ZBTB11.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ZBTB11.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ZBTB11.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ZBTB11.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ZBTB11.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64291_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ZBTB11_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ZBTB11.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62870_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/ZBTB11_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64291_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_ZBTB11_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62870_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_ZBTB11_MCF-7.bed

echo TRIM24 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64819_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2405_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/TRIM24.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/TRIM24.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/TRIM24.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/TRIM24.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/TRIM24.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/TRIM24.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64819_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/TRIM24_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/TRIM24.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2405_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/TRIM24_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64819_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_TRIM24_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2405_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_TRIM24_MCF-7.bed

echo PML 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46301_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46307_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/PML.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/PML.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/PML.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/PML.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/PML.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/PML.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46301_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/PML_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/PML.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46307_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/PML_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46301_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_PML_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46307_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_PML_MCF-7.bed

echo DPF1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/77264_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/77255_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/DPF1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/DPF1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/DPF1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/DPF1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/DPF1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/DPF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/77264_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/DPF1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/DPF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/77255_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/DPF1_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/77264_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_DPF1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/77255_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_DPF1_MCF-7.bed

echo DPF2 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64319_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63061_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/DPF2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/DPF2.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/DPF2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/DPF2.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/DPF2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/DPF2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64319_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/DPF2_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/DPF2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63061_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/DPF2_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64319_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_DPF2_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63061_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_DPF2_MCF-7.bed

echo CEBPG 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64053_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62690_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/CEBPG.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/CEBPG.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/CEBPG.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/CEBPG.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/CEBPG.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/CEBPG.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64053_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/CEBPG_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/CEBPG.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62690_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/CEBPG_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64053_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_CEBPG_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62690_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_CEBPG_MCF-7.bed

echo CTCFL 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46293_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/55785_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/CTCFL.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/CTCFL.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/CTCFL.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_MCF-7/CTCFL.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/CTCFL.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/CTCFL.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46293_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/CTCFL_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_MCF-7/CTCFL.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/55785_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_MCF-7/CTCFL_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46293_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_CTCFL_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/55785_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_MCF-7/_CTCFL_MCF-7.bed


#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=8000
#SBATCH -t 4:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f3_union_peaks_across_cellTypes/slurm_K562_vs_HeLa.out

echo POLR2A 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/33859_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/73296_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/POLR2A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HeLa/POLR2A.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/POLR2A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HeLa/POLR2A.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/POLR2A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HeLa/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/33859_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HeLa/POLR2A_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HeLa/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/73296_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HeLa/POLR2A_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/33859_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HeLa/_POLR2A_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/73296_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HeLa/_POLR2A_HeLa.bed

echo H3K4me3 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/72784_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/53321_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/H3K4me3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HeLa/H3K4me3.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/H3K4me3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HeLa/H3K4me3.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/H3K4me3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HeLa/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/72784_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HeLa/H3K4me3_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HeLa/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/53321_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HeLa/H3K4me3_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/72784_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HeLa/_H3K4me3_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/53321_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HeLa/_H3K4me3_HeLa.bed

echo H3K27ac 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/35400_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62017_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/H3K27ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HeLa/H3K27ac.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/H3K27ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HeLa/H3K27ac.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/H3K27ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HeLa/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/35400_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HeLa/H3K27ac_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HeLa/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62017_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HeLa/H3K27ac_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/35400_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HeLa/_H3K27ac_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62017_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HeLa/_H3K27ac_HeLa.bed

echo CTCF 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45899_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8087_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/CTCF.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HeLa/CTCF.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/CTCF.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HeLa/CTCF.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/CTCF.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HeLa/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45899_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HeLa/CTCF_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HeLa/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8087_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HeLa/CTCF_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45899_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HeLa/_CTCF_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8087_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HeLa/_CTCF_HeLa.bed

echo H3K4me2 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45382_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1252_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/H3K4me2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HeLa/H3K4me2.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/H3K4me2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HeLa/H3K4me2.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/H3K4me2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HeLa/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45382_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HeLa/H3K4me2_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HeLa/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1252_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HeLa/H3K4me2_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45382_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HeLa/_H3K4me2_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1252_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HeLa/_H3K4me2_HeLa.bed

echo EP300 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45963_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56234_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/EP300.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HeLa/EP300.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/EP300.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HeLa/EP300.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/EP300.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HeLa/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45963_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HeLa/EP300_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HeLa/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56234_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HeLa/EP300_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45963_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HeLa/_EP300_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56234_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HeLa/_EP300_HeLa.bed

echo H3K9ac 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45406_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/67816_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/H3K9ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HeLa/H3K9ac.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/H3K9ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HeLa/H3K9ac.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/H3K9ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HeLa/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45406_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HeLa/H3K9ac_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HeLa/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/67816_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HeLa/H3K9ac_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45406_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HeLa/_H3K9ac_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/67816_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HeLa/_H3K9ac_HeLa.bed

echo MYC 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45894_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/38115_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/MYC.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HeLa/MYC.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/MYC.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HeLa/MYC.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/MYC.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HeLa/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45894_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HeLa/MYC_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HeLa/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/38115_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HeLa/MYC_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45894_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HeLa/_MYC_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/38115_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HeLa/_MYC_HeLa.bed

echo GABPA 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46268_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/9212_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/GABPA.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HeLa/GABPA.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/GABPA.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HeLa/GABPA.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/GABPA.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HeLa/GABPA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46268_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HeLa/GABPA_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HeLa/GABPA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/9212_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HeLa/GABPA_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46268_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HeLa/_GABPA_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/9212_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HeLa/_GABPA_HeLa.bed

echo YY1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46303_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74510_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/YY1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HeLa/YY1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/YY1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HeLa/YY1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/YY1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HeLa/YY1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46303_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HeLa/YY1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HeLa/YY1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74510_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HeLa/YY1_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46303_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HeLa/_YY1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74510_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HeLa/_YY1_HeLa.bed

echo MAX 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46292_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/38119_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/MAX.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HeLa/MAX.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/MAX.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HeLa/MAX.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/MAX.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HeLa/MAX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46292_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HeLa/MAX_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HeLa/MAX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/38119_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HeLa/MAX_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46292_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HeLa/_MAX_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/38119_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HeLa/_MAX_HeLa.bed

echo H2AZ 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45416_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1231_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/H2AZ.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HeLa/H2AZ.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/H2AZ.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HeLa/H2AZ.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/H2AZ.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HeLa/H2AZ.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45416_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HeLa/H2AZ_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HeLa/H2AZ.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1231_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HeLa/H2AZ_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45416_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HeLa/_H2AZ_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1231_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HeLa/_H2AZ_HeLa.bed

echo SMC3 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45914_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/87043_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/SMC3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HeLa/SMC3.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/SMC3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HeLa/SMC3.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/SMC3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HeLa/SMC3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45914_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HeLa/SMC3_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HeLa/SMC3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/87043_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HeLa/SMC3_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45914_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HeLa/_SMC3_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/87043_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HeLa/_SMC3_HeLa.bed

echo ZNF143 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45946_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/9210_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/ZNF143.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HeLa/ZNF143.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/ZNF143.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HeLa/ZNF143.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/ZNF143.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HeLa/ZNF143.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45946_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HeLa/ZNF143_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HeLa/ZNF143.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/9210_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HeLa/ZNF143_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45946_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HeLa/_ZNF143_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/9210_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HeLa/_ZNF143_HeLa.bed

echo H3K9me3 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/56869_b_sort_peaks.broadPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/71058_b_sort_peaks.broadPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/H3K9me3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HeLa/H3K9me3.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/H3K9me3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HeLa/H3K9me3.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/H3K9me3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HeLa/H3K9me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/56869_b_sort_peaks.broadPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HeLa/H3K9me3_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HeLa/H3K9me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/71058_b_sort_peaks.broadPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HeLa/H3K9me3_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/56869_b_sort_peaks.broadPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HeLa/_H3K9me3_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/71058_b_sort_peaks.broadPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HeLa/_H3K9me3_HeLa.bed

echo HCFC1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45890_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/9205_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/HCFC1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HeLa/HCFC1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/HCFC1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HeLa/HCFC1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/HCFC1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HeLa/HCFC1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45890_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HeLa/HCFC1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HeLa/HCFC1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/9205_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HeLa/HCFC1_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45890_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HeLa/_HCFC1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/9205_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HeLa/_HCFC1_HeLa.bed

echo NR2C2 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64395_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/50977_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/NR2C2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HeLa/NR2C2.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/NR2C2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HeLa/NR2C2.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/NR2C2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HeLa/NR2C2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64395_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HeLa/NR2C2_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HeLa/NR2C2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/50977_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HeLa/NR2C2_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64395_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HeLa/_NR2C2_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/50977_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HeLa/_NR2C2_HeLa.bed

echo PHF8 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8790_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1314_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/PHF8.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HeLa/PHF8.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/PHF8.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HeLa/PHF8.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/PHF8.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HeLa/PHF8.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8790_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HeLa/PHF8_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HeLa/PHF8.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1314_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HeLa/PHF8_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8790_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HeLa/_PHF8_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1314_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HeLa/_PHF8_HeLa.bed

echo MBD2 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63028_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/40334_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/MBD2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HeLa/MBD2.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/MBD2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HeLa/MBD2.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/MBD2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HeLa/MBD2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63028_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HeLa/MBD2_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HeLa/MBD2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/40334_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HeLa/MBD2_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63028_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HeLa/_MBD2_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/40334_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HeLa/_MBD2_HeLa.bed

echo GTF2B 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45932_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56961_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/GTF2B.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HeLa/GTF2B.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/GTF2B.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HeLa/GTF2B.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/GTF2B.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HeLa/GTF2B.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45932_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HeLa/GTF2B_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HeLa/GTF2B.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56961_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HeLa/GTF2B_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45932_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HeLa/_GTF2B_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56961_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HeLa/_GTF2B_HeLa.bed

echo SUPT5H 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64773_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/53325_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/SUPT5H.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_HeLa/SUPT5H.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/SUPT5H.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_HeLa/SUPT5H.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_HeLa/SUPT5H.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HeLa/SUPT5H.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64773_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HeLa/SUPT5H_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_HeLa/SUPT5H.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/53325_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_HeLa/SUPT5H_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64773_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HeLa/_SUPT5H_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/53325_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_HeLa/_SUPT5H_HeLa.bed


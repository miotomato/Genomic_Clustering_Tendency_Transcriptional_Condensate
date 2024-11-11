#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=8000
#SBATCH -t 4:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f3_union_peaks_across_cellTypes/slurm_GM12878_vs_HeLa.out

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46137_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/73296_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/POLR2A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/POLR2A.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/POLR2A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/POLR2A.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/POLR2A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46137_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/POLR2A_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/73296_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/POLR2A_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46137_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/_POLR2A_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/73296_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/_POLR2A_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42291_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/53321_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/H3K4me3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/H3K4me3.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/H3K4me3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/H3K4me3.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/H3K4me3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42291_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/H3K4me3_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/53321_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/H3K4me3_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42291_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/_H3K4me3_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/53321_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/_H3K4me3_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/48285_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62017_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/H3K27ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/H3K27ac.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/H3K27ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/H3K27ac.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/H3K27ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/48285_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/H3K27ac_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62017_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/H3K27ac_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/48285_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/_H3K27ac_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62017_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/_H3K27ac_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45631_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8087_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/CTCF.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/CTCF.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/CTCF.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/CTCF.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/CTCF.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45631_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/CTCF_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8087_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/CTCF_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45631_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/_CTCF_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8087_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/_CTCF_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45214_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1252_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/H3K4me2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/H3K4me2.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/H3K4me2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/H3K4me2.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/H3K4me2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45214_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/H3K4me2_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1252_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/H3K4me2_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45214_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/_H3K4me2_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1252_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/_H3K4me2_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/48272_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56234_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/EP300.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/EP300.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/EP300.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/EP300.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/EP300.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/48272_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/EP300_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56234_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/EP300_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/48272_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/_EP300_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56234_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/_EP300_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/36652_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/67816_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/H3K9ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/H3K9ac.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/H3K9ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/H3K9ac.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/H3K9ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/36652_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/H3K9ac_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/67816_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/H3K9ac_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/36652_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/_H3K9ac_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/67816_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/_H3K9ac_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/666_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/38115_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/MYC.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/MYC.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/MYC.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/MYC.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/MYC.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/666_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/MYC_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/38115_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/MYC_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/666_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/_MYC_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/38115_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/_MYC_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46136_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/9212_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/GABPA.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/GABPA.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/GABPA.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/GABPA.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/GABPA.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/GABPA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46136_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/GABPA_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/GABPA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/9212_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/GABPA_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46136_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/_GABPA_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/9212_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/_GABPA_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/7060_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74510_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/YY1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/YY1.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/YY1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/YY1.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/YY1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/YY1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/7060_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/YY1_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/YY1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74510_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/YY1_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/7060_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/_YY1_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74510_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/_YY1_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/661_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/38119_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/MAX.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/MAX.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/MAX.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/MAX.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/MAX.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/MAX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/661_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/MAX_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/MAX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/38119_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/MAX_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/661_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/_MAX_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/38119_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/_MAX_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45221_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1231_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/H2AZ.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/H2AZ.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/H2AZ.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/H2AZ.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/H2AZ.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/H2AZ.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45221_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/H2AZ_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/H2AZ.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1231_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/H2AZ_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45221_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/_H2AZ_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1231_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/_H2AZ_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45640_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/87043_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/SMC3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/SMC3.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/SMC3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/SMC3.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/SMC3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/SMC3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45640_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/SMC3_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/SMC3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/87043_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/SMC3_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45640_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/_SMC3_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/87043_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/_SMC3_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45629_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/9210_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/ZNF143.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/ZNF143.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/ZNF143.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/ZNF143.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/ZNF143.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/ZNF143.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45629_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/ZNF143_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/ZNF143.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/9210_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/ZNF143_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45629_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/_ZNF143_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/9210_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/_ZNF143_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45641_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1336_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/RELA.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/RELA.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/RELA.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/RELA.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/RELA.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/RELA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45641_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/RELA_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/RELA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1336_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/RELA_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45641_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/_RELA_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1336_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/_RELA_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/669_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/50977_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/NR2C2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/NR2C2.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/NR2C2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/NR2C2.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/NR2C2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/NR2C2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/669_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/NR2C2_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/NR2C2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/50977_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/NR2C2_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/669_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/_NR2C2_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/50977_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_HeLa/_NR2C2_HeLa.bed


#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=8000
#SBATCH -t 4:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f3_union_peaks_across_cellTypes/slurm_HCT-116_vs_Jurkat.out

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46201_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/5199_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/POLR2A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/POLR2A.cat.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/POLR2A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/POLR2A.sort.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/POLR2A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46201_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/POLR2A_merge_intersect_HCT-116.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/5199_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/POLR2A_merge_intersect_Jurkat.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46201_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/_POLR2A_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/5199_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/_POLR2A_Jurkat.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57094_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/49527_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/H3K4me3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/H3K4me3.cat.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/H3K4me3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/H3K4me3.sort.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/H3K4me3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57094_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/H3K4me3_merge_intersect_HCT-116.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/49527_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/H3K4me3_merge_intersect_Jurkat.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57094_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/_H3K4me3_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/49527_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/_H3K4me3_Jurkat.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62110_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/54237_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/H3K27ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/H3K27ac.cat.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/H3K27ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/H3K27ac.sort.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/H3K27ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62110_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/H3K27ac_merge_intersect_HCT-116.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/54237_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/H3K27ac_merge_intersect_Jurkat.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62110_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/_H3K27ac_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/54237_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/_H3K27ac_Jurkat.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/42149_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/54151_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/CTCF.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/CTCF.cat.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/CTCF.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/CTCF.sort.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/CTCF.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/42149_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/CTCF_merge_intersect_HCT-116.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/54151_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/CTCF_merge_intersect_Jurkat.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/42149_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/_CTCF_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/54151_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/_CTCF_Jurkat.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42910_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/74134_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/H3K4me2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/H3K4me2.cat.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/H3K4me2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/H3K4me2.sort.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/H3K4me2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42910_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/H3K4me2_merge_intersect_HCT-116.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/74134_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/H3K4me2_merge_intersect_Jurkat.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42910_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/_H3K4me2_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/74134_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/_H3K4me2_Jurkat.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/70806_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/69616_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/MYC.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/MYC.cat.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/MYC.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/MYC.sort.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/MYC.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/70806_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/MYC_merge_intersect_HCT-116.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/69616_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/MYC_merge_intersect_Jurkat.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/70806_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/_MYC_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/69616_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/_MYC_Jurkat.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46203_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/47863_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/REST.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/REST.cat.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/REST.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/REST.sort.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/REST.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/REST.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46203_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/REST_merge_intersect_HCT-116.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/REST.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/47863_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/REST_merge_intersect_Jurkat.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46203_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/_REST_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/47863_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/_REST_Jurkat.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/68783_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/69617_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/BRD4.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/BRD4.cat.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/BRD4.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/BRD4.sort.bed > f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/BRD4.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/BRD4.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/68783_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/BRD4_merge_intersect_HCT-116.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/BRD4.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/69617_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/BRD4_merge_intersect_Jurkat.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/68783_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/_BRD4_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/69617_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HCT-116_vs_Jurkat/_BRD4_Jurkat.bed


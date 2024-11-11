#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=8000
#SBATCH -t 4:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f3_union_peaks_across_cellTypes/slurm_HepG2_vs_Jurkat.out

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46256_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/5199_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/POLR2A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/POLR2A.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/POLR2A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/POLR2A.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/POLR2A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46256_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/POLR2A_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/5199_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/POLR2A_merge_intersect_Jurkat.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46256_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/_POLR2A_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/5199_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/_POLR2A_Jurkat.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45306_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/49527_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/H3K4me3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/H3K4me3.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/H3K4me3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/H3K4me3.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/H3K4me3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45306_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/H3K4me3_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/49527_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/H3K4me3_merge_intersect_Jurkat.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45306_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/_H3K4me3_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/49527_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/_H3K4me3_Jurkat.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45299_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/54237_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/H3K27ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/H3K27ac.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/H3K27ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/H3K27ac.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/H3K27ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45299_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/H3K27ac_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/54237_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/H3K27ac_merge_intersect_Jurkat.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45299_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/_H3K27ac_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/54237_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/_H3K27ac_Jurkat.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1380_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/54151_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/CTCF.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/CTCF.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/CTCF.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/CTCF.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/CTCF.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1380_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/CTCF_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/54151_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/CTCF_merge_intersect_Jurkat.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1380_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/_CTCF_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/54151_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/_CTCF_Jurkat.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45307_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/74134_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/H3K4me2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/H3K4me2.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/H3K4me2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/H3K4me2.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/H3K4me2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45307_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/H3K4me2_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/74134_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/H3K4me2_merge_intersect_Jurkat.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45307_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/_H3K4me2_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/74134_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/_H3K4me2_Jurkat.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45827_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/69616_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/MYC.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/MYC.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/MYC.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/MYC.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/MYC.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45827_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/MYC_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/69616_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/MYC_merge_intersect_Jurkat.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45827_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/_MYC_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/69616_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/_MYC_Jurkat.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46257_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/47863_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/REST.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/REST.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/REST.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/REST.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/REST.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/REST.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46257_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/REST_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/REST.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/47863_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/REST_merge_intersect_Jurkat.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46257_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/_REST_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/47863_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/_REST_Jurkat.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56045_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/47413_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/GABPA.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/GABPA.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/GABPA.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/GABPA.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/GABPA.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/GABPA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56045_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/GABPA_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/GABPA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/47413_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/GABPA_merge_intersect_Jurkat.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56045_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/_GABPA_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/47413_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_Jurkat/_GABPA_Jurkat.bed


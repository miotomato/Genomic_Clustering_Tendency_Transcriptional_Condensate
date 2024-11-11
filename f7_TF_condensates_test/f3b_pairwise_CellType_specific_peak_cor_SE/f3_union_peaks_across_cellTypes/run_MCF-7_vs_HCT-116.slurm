#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=8000
#SBATCH -t 4:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f3_union_peaks_across_cellTypes/slurm_MCF-7_vs_HCT-116.out

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45985_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46201_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/POLR2A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/POLR2A.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/POLR2A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/POLR2A.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/POLR2A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45985_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/POLR2A_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46201_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/POLR2A_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45985_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_POLR2A_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46201_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_POLR2A_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/47270_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57094_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H3K4me3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H3K4me3.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H3K4me3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H3K4me3.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H3K4me3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/47270_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H3K4me3_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57094_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H3K4me3_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/47270_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_H3K4me3_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57094_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_H3K4me3_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/33141_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62110_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H3K27ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H3K27ac.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H3K27ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H3K27ac.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H3K27ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/33141_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H3K27ac_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62110_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H3K27ac_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/33141_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_H3K27ac_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62110_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_H3K27ac_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/83736_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/42149_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/CTCF.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/CTCF.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/CTCF.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/CTCF.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/CTCF.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/83736_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/CTCF_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/42149_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/CTCF_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/83736_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_CTCF_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/42149_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_CTCF_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/75112_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42910_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H3K4me2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H3K4me2.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H3K4me2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H3K4me2.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H3K4me2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/75112_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H3K4me2_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42910_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H3K4me2_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/75112_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_H3K4me2_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42910_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_H3K4me2_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/33155_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/42907_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/EP300.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/EP300.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/EP300.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/EP300.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/EP300.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/33155_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/EP300_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/42907_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/EP300_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/33155_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_EP300_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/42907_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_EP300_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/86292_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/87280_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H3K9ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H3K9ac.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H3K9ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H3K9ac.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H3K9ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/86292_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H3K9ac_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/87280_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H3K9ac_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/86292_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_H3K9ac_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/87280_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_H3K9ac_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/38193_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/70806_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/MYC.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/MYC.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/MYC.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/MYC.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/MYC.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/38193_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/MYC_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/70806_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/MYC_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/38193_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_MYC_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/70806_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_MYC_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/33151_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57095_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H3K4me1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H3K4me1.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H3K4me1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H3K4me1.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H3K4me1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H3K4me1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/33151_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H3K4me1_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H3K4me1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57095_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H3K4me1_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/33151_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_H3K4me1_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57095_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_H3K4me1_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46313_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46203_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/REST.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/REST.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/REST.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/REST.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/REST.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/REST.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46313_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/REST_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/REST.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46203_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/REST_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46313_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_REST_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46203_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_REST_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46324_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46216_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/MAX.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/MAX.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/MAX.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/MAX.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/MAX.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/MAX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46324_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/MAX_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/MAX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46216_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/MAX_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46324_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_MAX_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46216_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_MAX_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2387_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46207_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/RAD21.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/RAD21.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/RAD21.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/RAD21.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/RAD21.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/RAD21.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2387_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/RAD21_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/RAD21.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46207_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/RAD21_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2387_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_RAD21_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46207_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_RAD21_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/76117_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/83065_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H2AZ.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H2AZ.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H2AZ.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H2AZ.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H2AZ.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H2AZ.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/76117_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H2AZ_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H2AZ.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/83065_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/H2AZ_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/76117_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_H2AZ_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/83065_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_H2AZ_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46316_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46206_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/CEBPB.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/CEBPB.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/CEBPB.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/CEBPB.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/CEBPB.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/CEBPB.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46316_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/CEBPB_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/CEBPB.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46206_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/CEBPB_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46316_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_CEBPB_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46206_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_CEBPB_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46308_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46213_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/JUND.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/JUND.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/JUND.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/JUND.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/JUND.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/JUND.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46308_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/JUND_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/JUND.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46213_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/JUND_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46308_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_JUND_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46213_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_JUND_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46315_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46210_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/TEAD4.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/TEAD4.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/TEAD4.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/TEAD4.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/TEAD4.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/TEAD4.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46315_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/TEAD4_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/TEAD4.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46210_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/TEAD4_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46315_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_TEAD4_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46210_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_TEAD4_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46312_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46211_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/SRF.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/SRF.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/SRF.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/SRF.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/SRF.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/SRF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46312_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/SRF_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/SRF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46211_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/SRF_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46312_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_SRF_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46211_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_SRF_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45979_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45714_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/TCF7L2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/TCF7L2.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/TCF7L2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/TCF7L2.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/TCF7L2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/TCF7L2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45979_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/TCF7L2_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/TCF7L2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45714_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/TCF7L2_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45979_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_TCF7L2_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45714_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_TCF7L2_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46323_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46202_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/ELF1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/ELF1.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/ELF1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/ELF1.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/ELF1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/ELF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46323_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/ELF1_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/ELF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46202_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/ELF1_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46323_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_ELF1_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46202_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_ELF1_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46327_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46208_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/SIN3A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/SIN3A.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/SIN3A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/SIN3A.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/SIN3A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/SIN3A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46327_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/SIN3A_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/SIN3A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46208_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/SIN3A_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46327_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_SIN3A_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46208_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_SIN3A_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/71860_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/85289_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/SMC1A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/SMC1A.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/SMC1A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/SMC1A.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/SMC1A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/SMC1A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/71860_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/SMC1A_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/SMC1A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/85289_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/SMC1A_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/71860_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_SMC1A_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/85289_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_SMC1A_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46319_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46214_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/EGR1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/EGR1.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/EGR1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/EGR1.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/EGR1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/EGR1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46319_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/EGR1_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/EGR1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46214_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/EGR1_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46319_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_EGR1_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46214_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_EGR1_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/71862_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/84953_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/NIPBL.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/NIPBL.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/NIPBL.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/NIPBL.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/NIPBL.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/NIPBL.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/71862_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/NIPBL_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/NIPBL.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/84953_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/NIPBL_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/71862_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_NIPBL_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/84953_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_NIPBL_HCT-116.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/84928_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/81160_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/ZFX.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/ZFX.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/ZFX.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/ZFX.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/ZFX.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/ZFX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/84928_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/ZFX_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/ZFX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/81160_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/ZFX_merge_intersect_HCT-116.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/84928_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_ZFX_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/81160_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HCT-116/_ZFX_HCT-116.bed


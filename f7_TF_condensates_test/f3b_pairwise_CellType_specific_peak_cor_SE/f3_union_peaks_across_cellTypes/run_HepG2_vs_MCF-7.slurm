#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=8000
#SBATCH -t 4:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f3_union_peaks_across_cellTypes/slurm_HepG2_vs_MCF-7.out

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46256_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45985_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/POLR2A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/POLR2A.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/POLR2A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/POLR2A.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/POLR2A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46256_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/POLR2A_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45985_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/POLR2A_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46256_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_POLR2A_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45985_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_POLR2A_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45306_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/47270_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K4me3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K4me3.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K4me3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K4me3.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K4me3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45306_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K4me3_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/47270_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K4me3_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45306_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_H3K4me3_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/47270_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_H3K4me3_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45299_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/33141_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K27ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K27ac.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K27ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K27ac.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K27ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45299_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K27ac_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/33141_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K27ac_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45299_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_H3K27ac_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/33141_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_H3K27ac_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1380_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/83736_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/CTCF.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/CTCF.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/CTCF.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/CTCF.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/CTCF.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1380_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/CTCF_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/83736_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/CTCF_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1380_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_CTCF_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/83736_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_CTCF_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45307_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/75112_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K4me2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K4me2.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K4me2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K4me2.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K4me2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45307_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K4me2_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/75112_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K4me2_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45307_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_H3K4me2_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/75112_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_H3K4me2_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46225_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/33155_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/EP300.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/EP300.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/EP300.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/EP300.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/EP300.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46225_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/EP300_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/33155_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/EP300_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46225_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_EP300_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/33155_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_EP300_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45305_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/86292_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K9ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K9ac.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K9ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K9ac.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K9ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45305_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K9ac_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/86292_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K9ac_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45305_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_H3K9ac_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/86292_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_H3K9ac_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45827_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/38193_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/MYC.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/MYC.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/MYC.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/MYC.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/MYC.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45827_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/MYC_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/38193_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/MYC_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45827_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_MYC_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/38193_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_MYC_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/51255_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/33151_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K4me1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K4me1.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K4me1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K4me1.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K4me1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K4me1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/51255_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K4me1_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K4me1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/33151_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K4me1_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/51255_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_H3K4me1_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/33151_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_H3K4me1_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46257_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46313_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/REST.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/REST.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/REST.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/REST.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/REST.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/REST.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46257_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/REST_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/REST.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46313_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/REST_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46257_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_REST_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46313_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_REST_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56045_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46314_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/GABPA.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/GABPA.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/GABPA.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/GABPA.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/GABPA.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/GABPA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56045_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/GABPA_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/GABPA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46314_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/GABPA_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56045_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_GABPA_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46314_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_GABPA_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46244_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46324_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/MAX.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/MAX.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/MAX.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/MAX.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/MAX.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/MAX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46244_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/MAX_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/MAX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46324_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/MAX_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46244_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_MAX_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46324_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_MAX_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46252_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2387_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/RAD21.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/RAD21.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/RAD21.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/RAD21.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/RAD21.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/RAD21.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46252_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/RAD21_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/RAD21.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2387_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/RAD21_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46252_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_RAD21_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2387_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_RAD21_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46258_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46320_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/TAF1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/TAF1.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/TAF1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/TAF1.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/TAF1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/TAF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46258_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/TAF1_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/TAF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46320_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/TAF1_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46258_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_TAF1_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46320_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_TAF1_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45814_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46316_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/CEBPB.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/CEBPB.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/CEBPB.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/CEBPB.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/CEBPB.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/CEBPB.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45814_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/CEBPB_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/CEBPB.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46316_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/CEBPB_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45814_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_CEBPB_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46316_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_CEBPB_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45828_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46308_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/JUND.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/JUND.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/JUND.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/JUND.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/JUND.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/JUND.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45828_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/JUND_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/JUND.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46308_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/JUND_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45828_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_JUND_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46308_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_JUND_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/73725_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74153_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/NR3C1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/NR3C1.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/NR3C1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/NR3C1.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/NR3C1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/NR3C1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/73725_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/NR3C1_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/NR3C1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74153_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/NR3C1_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/73725_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_NR3C1_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74153_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_NR3C1_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46238_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46315_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/TEAD4.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/TEAD4.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/TEAD4.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/TEAD4.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/TEAD4.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/TEAD4.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46238_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/TEAD4_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/TEAD4.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46315_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/TEAD4_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46238_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_TEAD4_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46315_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_TEAD4_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45300_b_sort_peaks.broadPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/49474_b_sort_peaks.broadPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K9me3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K9me3.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K9me3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K9me3.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K9me3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K9me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45300_b_sort_peaks.broadPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K9me3_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K9me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/49474_b_sort_peaks.broadPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/H3K9me3_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45300_b_sort_peaks.broadPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_H3K9me3_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/49474_b_sort_peaks.broadPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_H3K9me3_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46228_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46321_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/TCF12.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/TCF12.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/TCF12.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/TCF12.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/TCF12.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/TCF12.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46228_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/TCF12_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/TCF12.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46321_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/TCF12_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46228_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_TCF12_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46321_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_TCF12_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/7323_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/33143_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/FOXA1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/FOXA1.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/FOXA1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/FOXA1.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/FOXA1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/FOXA1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/7323_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/FOXA1_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/FOXA1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/33143_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/FOXA1_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/7323_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_FOXA1_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/33143_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_FOXA1_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64667_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74638_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/NRF1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/NRF1.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/NRF1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/NRF1.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/NRF1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/NRF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64667_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/NRF1_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/NRF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74638_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/NRF1_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64667_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_NRF1_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74638_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_NRF1_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46223_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46312_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/SRF.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/SRF.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/SRF.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/SRF.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/SRF.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/SRF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46223_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/SRF_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/SRF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46312_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/SRF_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46223_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_SRF_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46312_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_SRF_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45825_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45979_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/TCF7L2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/TCF7L2.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/TCF7L2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/TCF7L2.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/TCF7L2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/TCF7L2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45825_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/TCF7L2_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/TCF7L2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45979_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/TCF7L2_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45825_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_TCF7L2_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45979_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_TCF7L2_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46254_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46323_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/ELF1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/ELF1.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/ELF1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/ELF1.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/ELF1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/ELF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46254_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/ELF1_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/ELF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46323_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/ELF1_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46254_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_ELF1_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46323_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_ELF1_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/71877_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/49652_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/MED1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/MED1.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/MED1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/MED1.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/MED1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/MED1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/71877_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/MED1_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/MED1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/49652_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/MED1_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/71877_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_MED1_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/49652_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_MED1_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46236_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46327_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/SIN3A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/SIN3A.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/SIN3A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/SIN3A.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/SIN3A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/SIN3A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46236_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/SIN3A_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/SIN3A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46327_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/SIN3A_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46236_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_SIN3A_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46327_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_SIN3A_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/71876_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/71860_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/SMC1A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/SMC1A.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/SMC1A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/SMC1A.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/SMC1A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/SMC1A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/71876_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/SMC1A_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/SMC1A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/71860_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/SMC1A_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/71876_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_SMC1A_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/71860_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_SMC1A_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46235_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46326_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/HDAC2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/HDAC2.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/HDAC2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/HDAC2.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/HDAC2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/HDAC2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46235_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/HDAC2_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/HDAC2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46326_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/HDAC2_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46235_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_HDAC2_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46326_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_HDAC2_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1492_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2393_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/STAG1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/STAG1.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/STAG1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/STAG1.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/STAG1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/STAG1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1492_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/STAG1_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/STAG1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2393_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/STAG1_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1492_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_STAG1_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2393_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_STAG1_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46243_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46317_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/FOSL2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/FOSL2.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/FOSL2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/FOSL2.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/FOSL2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/FOSL2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46243_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/FOSL2_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/FOSL2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46317_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/FOSL2_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46243_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_FOSL2_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46317_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_FOSL2_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63800_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63278_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/ZBTB40.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/ZBTB40.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/ZBTB40.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/ZBTB40.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/ZBTB40.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/ZBTB40.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63800_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/ZBTB40_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/ZBTB40.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63278_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/ZBTB40_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63800_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_ZBTB40_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63278_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_ZBTB40_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46249_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46309_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/NR2F2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/NR2F2.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/NR2F2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/NR2F2.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/NR2F2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/NR2F2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46249_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/NR2F2_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/NR2F2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46309_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/NR2F2_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46249_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_NR2F2_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46309_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_NR2F2_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74686_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74639_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/KLF10.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/KLF10.cat.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/KLF10.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/KLF10.sort.bed > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/KLF10.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/KLF10.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74686_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/KLF10_merge_intersect_HepG2.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/KLF10.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74639_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/KLF10_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74686_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_KLF10_HepG2.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74639_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/HepG2_vs_MCF-7/_KLF10_MCF-7.bed


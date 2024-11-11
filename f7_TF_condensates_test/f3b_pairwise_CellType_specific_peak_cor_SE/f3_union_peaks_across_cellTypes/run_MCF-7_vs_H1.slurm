#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=8000
#SBATCH -t 4:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f3_union_peaks_across_cellTypes/slurm_MCF-7_vs_H1.out

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45985_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46199_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/POLR2A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/POLR2A.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/POLR2A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/POLR2A.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/POLR2A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45985_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/POLR2A_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46199_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/POLR2A_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45985_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_POLR2A_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46199_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_POLR2A_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/47270_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/908_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K4me3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K4me3.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K4me3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K4me3.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K4me3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/47270_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K4me3_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/908_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K4me3_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/47270_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_H3K4me3_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/908_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_H3K4me3_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/33141_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/6531_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K27ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K27ac.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K27ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K27ac.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K27ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/33141_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K27ac_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/6531_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K27ac_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/33141_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_H3K27ac_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/6531_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_H3K27ac_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/83736_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/50987_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/CTCF.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/CTCF.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/CTCF.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/CTCF.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/CTCF.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/83736_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/CTCF_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/50987_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/CTCF_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/83736_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_CTCF_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/50987_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_CTCF_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/75112_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/897_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K4me2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K4me2.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K4me2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K4me2.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K4me2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/75112_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K4me2_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/897_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K4me2_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/75112_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_H3K4me2_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/897_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_H3K4me2_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/33155_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46165_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/EP300.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/EP300.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/EP300.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/EP300.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/EP300.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/33155_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/EP300_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46165_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/EP300_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/33155_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_EP300_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46165_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_EP300_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/86292_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45227_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K9ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K9ac.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K9ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K9ac.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K9ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/86292_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K9ac_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45227_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K9ac_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/86292_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_H3K9ac_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45227_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_H3K9ac_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/38193_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45691_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/MYC.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/MYC.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/MYC.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/MYC.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/MYC.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/38193_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/MYC_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45691_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/MYC_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/38193_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_MYC_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45691_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_MYC_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/33151_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/894_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K4me1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K4me1.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K4me1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K4me1.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K4me1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K4me1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/33151_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K4me1_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K4me1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/894_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K4me1_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/33151_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_H3K4me1_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/894_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_H3K4me1_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46313_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46172_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/REST.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/REST.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/REST.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/REST.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/REST.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/REST.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46313_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/REST_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/REST.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46172_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/REST_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46313_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_REST_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46172_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_REST_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46314_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46185_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/GABPA.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/GABPA.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/GABPA.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/GABPA.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/GABPA.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/GABPA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46314_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/GABPA_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/GABPA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46185_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/GABPA_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46314_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_GABPA_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46185_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_GABPA_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46324_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46180_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/MAX.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/MAX.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/MAX.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/MAX.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/MAX.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/MAX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46324_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/MAX_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/MAX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46180_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/MAX_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46324_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_MAX_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46180_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_MAX_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2387_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46164_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/RAD21.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/RAD21.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/RAD21.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/RAD21.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/RAD21.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/RAD21.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2387_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/RAD21_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/RAD21.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46164_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/RAD21_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2387_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_RAD21_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46164_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_RAD21_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46320_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46196_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/TAF1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/TAF1.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/TAF1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/TAF1.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/TAF1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/TAF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46320_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/TAF1_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/TAF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46196_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/TAF1_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46320_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_TAF1_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46196_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_TAF1_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46316_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45681_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/CEBPB.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/CEBPB.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/CEBPB.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/CEBPB.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/CEBPB.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/CEBPB.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46316_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/CEBPB_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/CEBPB.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45681_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/CEBPB_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46316_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_CEBPB_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45681_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_CEBPB_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46308_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45696_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/JUND.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/JUND.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/JUND.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/JUND.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/JUND.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/JUND.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46308_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/JUND_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/JUND.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45696_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/JUND_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46308_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_JUND_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45696_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_JUND_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/71042_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45688_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/ZNF143.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/ZNF143.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/ZNF143.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/ZNF143.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/ZNF143.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/ZNF143.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/71042_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/ZNF143_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/ZNF143.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45688_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/ZNF143_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/71042_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_ZNF143_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45688_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_ZNF143_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46315_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46189_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/TEAD4.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/TEAD4.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/TEAD4.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/TEAD4.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/TEAD4.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/TEAD4.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46315_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/TEAD4_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/TEAD4.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46189_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/TEAD4_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46315_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_TEAD4_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46189_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_TEAD4_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/49474_b_sort_peaks.broadPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/48924_b_sort_peaks.broadPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K9me3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K9me3.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K9me3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K9me3.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K9me3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K9me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/49474_b_sort_peaks.broadPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K9me3_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K9me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/48924_b_sort_peaks.broadPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K9me3_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/49474_b_sort_peaks.broadPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_H3K9me3_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/48924_b_sort_peaks.broadPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_H3K9me3_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46321_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46166_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/TCF12.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/TCF12.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/TCF12.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/TCF12.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/TCF12.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/TCF12.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46321_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/TCF12_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/TCF12.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46166_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/TCF12_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46321_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_TCF12_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46166_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_TCF12_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74638_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45697_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/NRF1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/NRF1.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/NRF1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/NRF1.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/NRF1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/NRF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74638_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/NRF1_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/NRF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45697_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/NRF1_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74638_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_NRF1_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45697_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_NRF1_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46312_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46167_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/SRF.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/SRF.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/SRF.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/SRF.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/SRF.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/SRF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46312_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/SRF_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/SRF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46167_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/SRF_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46312_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_SRF_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46167_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_SRF_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46327_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45699_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/SIN3A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/SIN3A.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/SIN3A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/SIN3A.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/SIN3A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/SIN3A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46327_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/SIN3A_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/SIN3A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45699_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/SIN3A_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46327_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_SIN3A_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45699_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_SIN3A_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46326_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46181_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/HDAC2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/HDAC2.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/HDAC2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/HDAC2.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/HDAC2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/HDAC2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46326_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/HDAC2_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/HDAC2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46181_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/HDAC2_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46326_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_HDAC2_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46181_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_HDAC2_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46319_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46169_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/EGR1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/EGR1.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/EGR1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/EGR1.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/EGR1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/EGR1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46319_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/EGR1_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/EGR1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46169_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/EGR1_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46319_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_EGR1_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46169_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_EGR1_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/68129_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/38503_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K4ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K4ac.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K4ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K4ac.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K4ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K4ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/68129_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K4ac_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K4ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/38503_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H3K4ac_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/68129_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_H3K4ac_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/38503_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_H3K4ac_H1.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/48131_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/36797_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H2A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H2A.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H2A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H2A.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H2A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/48131_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H2A_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/36797_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_H1/H2A_merge_intersect_H1.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/48131_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_H2A_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/36797_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_H1/_H2A_H1.bed


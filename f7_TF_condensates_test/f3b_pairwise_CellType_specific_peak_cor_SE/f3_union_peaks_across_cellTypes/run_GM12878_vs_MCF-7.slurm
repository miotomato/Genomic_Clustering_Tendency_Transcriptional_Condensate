#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=8000
#SBATCH -t 4:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f3_union_peaks_across_cellTypes/slurm_GM12878_vs_MCF-7.out

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46137_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45985_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/POLR2A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/POLR2A.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/POLR2A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/POLR2A.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/POLR2A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46137_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/POLR2A_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45985_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/POLR2A_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46137_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_POLR2A_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45985_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_POLR2A_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42291_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/47270_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H3K4me3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H3K4me3.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H3K4me3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H3K4me3.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H3K4me3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42291_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H3K4me3_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/47270_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H3K4me3_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/42291_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_H3K4me3_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/47270_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_H3K4me3_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/48285_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/33141_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H3K27ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H3K27ac.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H3K27ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H3K27ac.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H3K27ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/48285_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H3K27ac_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/33141_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H3K27ac_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/48285_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_H3K27ac_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/33141_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_H3K27ac_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45631_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/83736_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/CTCF.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/CTCF.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/CTCF.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/CTCF.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/CTCF.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45631_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/CTCF_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/83736_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/CTCF_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45631_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_CTCF_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/83736_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_CTCF_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45214_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/75112_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H3K4me2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H3K4me2.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H3K4me2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H3K4me2.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H3K4me2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45214_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H3K4me2_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/75112_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H3K4me2_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45214_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_H3K4me2_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/75112_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_H3K4me2_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/48272_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/33155_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/EP300.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/EP300.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/EP300.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/EP300.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/EP300.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/48272_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/EP300_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/33155_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/EP300_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/48272_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_EP300_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/33155_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_EP300_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/36652_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/86292_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H3K9ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H3K9ac.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H3K9ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H3K9ac.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H3K9ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/36652_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H3K9ac_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/86292_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H3K9ac_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/36652_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_H3K9ac_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/86292_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_H3K9ac_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/666_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/38193_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/MYC.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/MYC.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/MYC.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/MYC.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/MYC.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/666_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/MYC_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/38193_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/MYC_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/666_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_MYC_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/38193_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_MYC_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/49297_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/33151_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H3K4me1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H3K4me1.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H3K4me1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H3K4me1.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H3K4me1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H3K4me1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/49297_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H3K4me1_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H3K4me1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/33151_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H3K4me1_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/49297_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_H3K4me1_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/33151_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_H3K4me1_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46113_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46313_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/REST.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/REST.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/REST.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/REST.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/REST.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/REST.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46113_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/REST_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/REST.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46313_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/REST_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46113_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_REST_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46313_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_REST_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46136_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46314_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/GABPA.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/GABPA.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/GABPA.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/GABPA.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/GABPA.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/GABPA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46136_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/GABPA_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/GABPA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46314_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/GABPA_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46136_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_GABPA_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46314_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_GABPA_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/661_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46324_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/MAX.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/MAX.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/MAX.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/MAX.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/MAX.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/MAX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/661_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/MAX_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/MAX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46324_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/MAX_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/661_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_MAX_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46324_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_MAX_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46150_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2387_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/RAD21.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/RAD21.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/RAD21.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/RAD21.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/RAD21.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/RAD21.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46150_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/RAD21_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/RAD21.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2387_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/RAD21_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46150_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_RAD21_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2387_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_RAD21_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46108_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46320_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/TAF1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/TAF1.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/TAF1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/TAF1.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/TAF1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/TAF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46108_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/TAF1_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/TAF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46320_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/TAF1_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46108_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_TAF1_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46320_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_TAF1_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45221_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/76117_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H2AZ.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H2AZ.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H2AZ.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H2AZ.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H2AZ.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H2AZ.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45221_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H2AZ_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H2AZ.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/76117_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/H2AZ_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45221_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_H2AZ_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/76117_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_H2AZ_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/660_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46308_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/JUND.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/JUND.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/JUND.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/JUND.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/JUND.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/JUND.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/660_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/JUND_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/JUND.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46308_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/JUND_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/660_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_JUND_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46308_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_JUND_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45629_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/71042_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/ZNF143.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/ZNF143.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/ZNF143.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/ZNF143.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/ZNF143.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/ZNF143.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45629_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/ZNF143_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/ZNF143.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/71042_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/ZNF143_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45629_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_ZNF143_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/71042_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_ZNF143_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45641_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74392_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/RELA.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/RELA.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/RELA.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/RELA.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/RELA.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/RELA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45641_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/RELA_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/RELA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74392_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/RELA_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45641_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_RELA_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74392_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_RELA_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46120_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46321_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/TCF12.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/TCF12.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/TCF12.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/TCF12.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/TCF12.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/TCF12.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46120_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/TCF12_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/TCF12.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46321_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/TCF12_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46120_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_TCF12_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46321_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_TCF12_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45624_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74638_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/NRF1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/NRF1.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/NRF1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/NRF1.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/NRF1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/NRF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45624_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/NRF1_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/NRF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74638_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/NRF1_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45624_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_NRF1_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74638_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_NRF1_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46148_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46312_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/SRF.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/SRF.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/SRF.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/SRF.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/SRF.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/SRF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46148_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/SRF_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/SRF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46312_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/SRF_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46148_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_SRF_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46312_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_SRF_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46127_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46323_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/ELF1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/ELF1.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/ELF1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/ELF1.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/ELF1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/ELF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46127_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/ELF1_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/ELF1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46323_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/ELF1_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46127_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_ELF1_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46323_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_ELF1_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45639_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46327_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/SIN3A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/SIN3A.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/SIN3A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/SIN3A.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/SIN3A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/SIN3A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45639_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/SIN3A_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/SIN3A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46327_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/SIN3A_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45639_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_SIN3A_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46327_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_SIN3A_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46134_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/39443_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/FOXM1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/FOXM1.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/FOXM1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/FOXM1.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/FOXM1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/FOXM1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46134_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/FOXM1_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/FOXM1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/39443_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/FOXM1_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46134_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_FOXM1_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/39443_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_FOXM1_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46126_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46319_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/EGR1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/EGR1.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/EGR1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/EGR1.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/EGR1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/EGR1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46126_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/EGR1_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/EGR1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46319_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/EGR1_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46126_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_EGR1_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46319_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_EGR1_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/42295_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2393_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/STAG1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/STAG1.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/STAG1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/STAG1.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/STAG1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/STAG1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/42295_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/STAG1_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/STAG1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2393_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/STAG1_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/42295_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_STAG1_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2393_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_STAG1_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45606_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2318_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/FOS.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/FOS.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/FOS.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/FOS.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/FOS.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/FOS.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45606_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/FOS_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/FOS.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2318_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/FOS_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45606_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_FOS_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2318_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_FOS_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45617_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/36807_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/E2F4.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/E2F4.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/E2F4.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/E2F4.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/E2F4.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/E2F4.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45617_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/E2F4_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/E2F4.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/36807_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/E2F4_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45617_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_E2F4_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/36807_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_E2F4_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62937_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63278_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/ZBTB40.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/ZBTB40.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/ZBTB40.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/ZBTB40.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/ZBTB40.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/ZBTB40.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62937_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/ZBTB40_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/ZBTB40.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63278_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/ZBTB40_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62937_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_ZBTB40_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63278_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_ZBTB40_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45636_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63433_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/ELK1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/ELK1.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/ELK1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/ELK1.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/ELK1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/ELK1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45636_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/ELK1_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/ELK1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63433_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/ELK1_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45636_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_ELK1_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/63433_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_ELK1_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45615_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62548_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/CUX1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/CUX1.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/CUX1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/CUX1.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/CUX1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/CUX1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45615_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/CUX1_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/CUX1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62548_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/CUX1_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45615_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_CUX1_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/62548_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_CUX1_MCF-7.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46143_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46307_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/PML.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/PML.cat.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/PML.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/PML.sort.bed > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/PML.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/PML.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46143_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/PML_merge_intersect_GM12878.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/PML.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46307_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/PML_merge_intersect_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46143_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_PML_GM12878.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46307_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/GM12878_vs_MCF-7/_PML_MCF-7.bed


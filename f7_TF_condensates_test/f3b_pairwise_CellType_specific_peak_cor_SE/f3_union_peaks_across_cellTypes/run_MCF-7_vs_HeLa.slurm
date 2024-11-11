#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=8000
#SBATCH -t 4:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f3_union_peaks_across_cellTypes/slurm_MCF-7_vs_HeLa.out

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45985_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/73296_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/POLR2A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/POLR2A.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/POLR2A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/POLR2A.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/POLR2A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45985_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/POLR2A_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/73296_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/POLR2A_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45985_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/_POLR2A_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/73296_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/_POLR2A_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/47270_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/53321_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H3K4me3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H3K4me3.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H3K4me3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H3K4me3.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H3K4me3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/47270_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H3K4me3_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/53321_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H3K4me3_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/47270_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/_H3K4me3_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/53321_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/_H3K4me3_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/33141_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62017_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H3K27ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H3K27ac.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H3K27ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H3K27ac.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H3K27ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/33141_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H3K27ac_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62017_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H3K27ac_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/33141_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/_H3K27ac_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/62017_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/_H3K27ac_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/83736_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8087_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/CTCF.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/CTCF.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/CTCF.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/CTCF.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/CTCF.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/83736_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/CTCF_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8087_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/CTCF_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/83736_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/_CTCF_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/8087_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/_CTCF_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/75112_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1252_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H3K4me2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H3K4me2.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H3K4me2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H3K4me2.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H3K4me2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/75112_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H3K4me2_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1252_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H3K4me2_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/75112_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/_H3K4me2_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/1252_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/_H3K4me2_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/33155_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56234_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/EP300.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/EP300.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/EP300.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/EP300.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/EP300.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/33155_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/EP300_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/EP300.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56234_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/EP300_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/33155_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/_EP300_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/56234_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/_EP300_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/86292_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/67816_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H3K9ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H3K9ac.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H3K9ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H3K9ac.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H3K9ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/86292_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H3K9ac_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H3K9ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/67816_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H3K9ac_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/86292_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/_H3K9ac_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/67816_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/_H3K9ac_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/38193_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/38115_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/MYC.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/MYC.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/MYC.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/MYC.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/MYC.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/38193_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/MYC_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/38115_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/MYC_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/38193_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/_MYC_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/38115_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/_MYC_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46314_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/9212_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/GABPA.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/GABPA.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/GABPA.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/GABPA.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/GABPA.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/GABPA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46314_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/GABPA_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/GABPA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/9212_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/GABPA_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46314_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/_GABPA_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/9212_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/_GABPA_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46324_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/38119_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/MAX.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/MAX.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/MAX.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/MAX.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/MAX.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/MAX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46324_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/MAX_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/MAX.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/38119_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/MAX_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46324_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/_MAX_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/38119_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/_MAX_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/76117_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1231_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H2AZ.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H2AZ.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H2AZ.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H2AZ.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H2AZ.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H2AZ.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/76117_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H2AZ_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H2AZ.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1231_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H2AZ_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/76117_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/_H2AZ_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1231_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/_H2AZ_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74153_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1300_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/NR3C1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/NR3C1.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/NR3C1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/NR3C1.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/NR3C1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/NR3C1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74153_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/NR3C1_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/NR3C1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1300_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/NR3C1_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74153_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/_NR3C1_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1300_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/_NR3C1_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/71042_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/9210_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/ZNF143.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/ZNF143.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/ZNF143.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/ZNF143.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/ZNF143.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/ZNF143.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/71042_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/ZNF143_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/ZNF143.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/9210_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/ZNF143_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/71042_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/_ZNF143_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/9210_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/_ZNF143_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74392_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1336_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/RELA.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/RELA.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/RELA.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/RELA.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/RELA.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/RELA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74392_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/RELA_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/RELA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1336_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/RELA_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/74392_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/_RELA_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1336_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/_RELA_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/49474_b_sort_peaks.broadPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/71058_b_sort_peaks.broadPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H3K9me3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H3K9me3.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H3K9me3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H3K9me3.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H3K9me3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H3K9me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/49474_b_sort_peaks.broadPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H3K9me3_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H3K9me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/71058_b_sort_peaks.broadPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H3K9me3_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/49474_b_sort_peaks.broadPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/_H3K9me3_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/71058_b_sort_peaks.broadPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/_H3K9me3_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2278_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1212_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/E2F1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/E2F1.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/E2F1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/E2F1.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/E2F1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/E2F1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2278_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/E2F1_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/E2F1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1212_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/E2F1_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2278_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/_E2F1_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1212_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/_E2F1_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64266_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/9205_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/HCFC1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/HCFC1.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/HCFC1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/HCFC1.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/HCFC1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/HCFC1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64266_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/HCFC1_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/HCFC1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/9205_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/HCFC1_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/64266_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/_HCFC1_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/9205_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/_HCFC1_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/69816_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/53324_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/NELFE.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/NELFE.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/NELFE.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/NELFE.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/NELFE.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/NELFE.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/69816_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/NELFE_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/NELFE.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/53324_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/NELFE_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/69816_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/_NELFE_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/53324_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/_NELFE_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/48825_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/40334_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/MBD2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/MBD2.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/MBD2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/MBD2.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/MBD2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/MBD2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/48825_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/MBD2_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/MBD2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/40334_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/MBD2_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/48825_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/_MBD2_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/40334_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/_MBD2_HeLa.bed

cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/48131_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/43235_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H2A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H2A.cat.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H2A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H2A.sort.bed > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H2A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/48131_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H2A_merge_intersect_MCF-7.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/43235_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/H2A_merge_intersect_HeLa.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/48131_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/_H2A_MCF-7.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/43235_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/MCF-7_vs_HeLa/_H2A_HeLa.bed


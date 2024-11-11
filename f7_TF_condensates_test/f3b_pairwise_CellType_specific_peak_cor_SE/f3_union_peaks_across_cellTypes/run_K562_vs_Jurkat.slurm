#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=8000
#SBATCH -t 4:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o f3_union_peaks_across_cellTypes/slurm_K562_vs_Jurkat.out

echo POLR2A 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/33859_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/5199_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/POLR2A.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_Jurkat/POLR2A.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/POLR2A.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_Jurkat/POLR2A.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/POLR2A.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_Jurkat/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/33859_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/POLR2A_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_Jurkat/POLR2A.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/5199_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/POLR2A_merge_intersect_Jurkat.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/33859_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_Jurkat/_POLR2A_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/5199_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_Jurkat/_POLR2A_Jurkat.bed

echo H3K4me3 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/72784_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/49527_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/H3K4me3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_Jurkat/H3K4me3.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/H3K4me3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_Jurkat/H3K4me3.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/H3K4me3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_Jurkat/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/72784_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/H3K4me3_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_Jurkat/H3K4me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/49527_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/H3K4me3_merge_intersect_Jurkat.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/72784_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_Jurkat/_H3K4me3_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/49527_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_Jurkat/_H3K4me3_Jurkat.bed

echo H3K27ac 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/35400_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/54237_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/H3K27ac.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_Jurkat/H3K27ac.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/H3K27ac.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_Jurkat/H3K27ac.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/H3K27ac.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_Jurkat/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/35400_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/H3K27ac_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_Jurkat/H3K27ac.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/54237_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/H3K27ac_merge_intersect_Jurkat.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/35400_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_Jurkat/_H3K27ac_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/54237_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_Jurkat/_H3K27ac_Jurkat.bed

echo CTCF 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45899_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/54151_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/CTCF.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_Jurkat/CTCF.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/CTCF.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_Jurkat/CTCF.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/CTCF.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_Jurkat/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45899_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/CTCF_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_Jurkat/CTCF.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/54151_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/CTCF_merge_intersect_Jurkat.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45899_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_Jurkat/_CTCF_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/54151_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_Jurkat/_CTCF_Jurkat.bed

echo H3K4me2 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45382_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/74134_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/H3K4me2.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_Jurkat/H3K4me2.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/H3K4me2.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_Jurkat/H3K4me2.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/H3K4me2.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_Jurkat/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45382_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/H3K4me2_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_Jurkat/H3K4me2.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/74134_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/H3K4me2_merge_intersect_Jurkat.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/45382_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_Jurkat/_H3K4me2_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/74134_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_Jurkat/_H3K4me2_Jurkat.bed

echo MYC 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45894_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/69616_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/MYC.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_Jurkat/MYC.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/MYC.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_Jurkat/MYC.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/MYC.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_Jurkat/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45894_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/MYC_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_Jurkat/MYC.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/69616_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/MYC_merge_intersect_Jurkat.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45894_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_Jurkat/_MYC_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/69616_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_Jurkat/_MYC_Jurkat.bed

echo REST 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46286_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/47863_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/REST.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_Jurkat/REST.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/REST.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_Jurkat/REST.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/REST.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_Jurkat/REST.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46286_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/REST_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_Jurkat/REST.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/47863_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/REST_merge_intersect_Jurkat.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46286_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_Jurkat/_REST_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/47863_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_Jurkat/_REST_Jurkat.bed

echo GABPA 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46268_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/47413_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/GABPA.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_Jurkat/GABPA.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/GABPA.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_Jurkat/GABPA.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/GABPA.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_Jurkat/GABPA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46268_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/GABPA_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_Jurkat/GABPA.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/47413_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/GABPA_merge_intersect_Jurkat.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/46268_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_Jurkat/_GABPA_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/47413_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_Jurkat/_GABPA_Jurkat.bed

echo RUNX1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2121_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/54234_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/RUNX1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_Jurkat/RUNX1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/RUNX1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_Jurkat/RUNX1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/RUNX1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_Jurkat/RUNX1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2121_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/RUNX1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_Jurkat/RUNX1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/54234_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/RUNX1_merge_intersect_Jurkat.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/2121_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_Jurkat/_RUNX1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/54234_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_Jurkat/_RUNX1_Jurkat.bed

echo H3K27me3 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57482_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/74127_b_sort_peaks.broadPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/H3K27me3.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_Jurkat/H3K27me3.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/H3K27me3.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_Jurkat/H3K27me3.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/H3K27me3.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_Jurkat/H3K27me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57482_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/H3K27me3_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_Jurkat/H3K27me3.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/74127_b_sort_peaks.broadPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/H3K27me3_merge_intersect_Jurkat.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/57482_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_Jurkat/_H3K27me3_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_hm/74127_b_sort_peaks.broadPeak.bed f3_union_peaks_across_cellTypes/K562_vs_Jurkat/_H3K27me3_Jurkat.bed

echo ERG 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/37639_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/47416_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/ERG.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_Jurkat/ERG.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/ERG.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_Jurkat/ERG.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/ERG.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_Jurkat/ERG.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/37639_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/ERG_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_Jurkat/ERG.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/47416_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/ERG_merge_intersect_Jurkat.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/37639_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_Jurkat/_ERG_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/47416_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_Jurkat/_ERG_Jurkat.bed

echo TAL1 
cat /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45872_sort_peaks.narrowPeak.bed \
/nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1957_sort_peaks.narrowPeak.bed > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/TAL1.cat.bed
bedtools sort -i f3_union_peaks_across_cellTypes/K562_vs_Jurkat/TAL1.cat.bed > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/TAL1.sort.bed
bedtools merge -i f3_union_peaks_across_cellTypes/K562_vs_Jurkat/TAL1.sort.bed > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/TAL1.merge.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_Jurkat/TAL1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45872_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/TAL1_merge_intersect_K562.bed
bedtools intersect -a f3_union_peaks_across_cellTypes/K562_vs_Jurkat/TAL1.merge.bed -b /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1957_sort_peaks.narrowPeak.bed -c > f3_union_peaks_across_cellTypes/K562_vs_Jurkat/TAL1_merge_intersect_Jurkat.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/45872_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_Jurkat/_TAL1_K562.bed
cp /nv/vol190/zanglab/zw5j/data/cistrome/cistrome_db_2019/human_factor/1957_sort_peaks.narrowPeak.bed f3_union_peaks_across_cellTypes/K562_vs_Jurkat/_TAL1_Jurkat.bed


#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=200000
#SBATCH -t 12:00:00
#SBATCH -p parallel
#SBATCH -A cphg_cz3d
#SBATCH -o out_py1.log3

# Run program

module load gcc/9.2.0
module load bedtools/2.29.2

cat cistrome_data/H3K27me3.cat.bed |cut -f1-3 > cistrome_data/H3K27me3.cat.cut.bed
bedtools sort -i cistrome_data/H3K27me3.cat.cut.bed > cistrome_data/H3K27me3.sort.bed
bedtools merge -i cistrome_data/H3K27me3.sort.bed > cistrome_data/H3K27me3.merge.bed

cat cistrome_data/H3K36me3.cat.bed |cut -f1-3 > cistrome_data/H3K36me3.cat.cut.bed
bedtools sort -i cistrome_data/H3K36me3.cat.cut.bed > cistrome_data/H3K36me3.sort.bed
bedtools merge -i cistrome_data/H3K36me3.sort.bed > cistrome_data/H3K36me3.merge.bed

# bedtools sort -i cistrome_data/H3K4me1.cat.bed > cistrome_data/H3K4me1.sort.bed
# bedtools merge -i cistrome_data/H3K4me1.sort.bed > cistrome_data/H3K4me1.merge.bed

cat cistrome_data/H3K9me3.cat.bed |cut -f1-3 > cistrome_data/H3K9me3.cat.cut.bed
bedtools sort -i cistrome_data/H3K9me3.cat.cut.bed > cistrome_data/H3K9me3.sort.bed
bedtools merge -i cistrome_data/H3K9me3.sort.bed > cistrome_data/H3K9me3.merge.bed



wc -l cistrome_data/H3K27me3.merge.bed
wc -l cistrome_data/H3K36me3.merge.bed
wc -l cistrome_data/H3K4me1.merge.bed
wc -l cistrome_data/H3K9me3.merge.bed






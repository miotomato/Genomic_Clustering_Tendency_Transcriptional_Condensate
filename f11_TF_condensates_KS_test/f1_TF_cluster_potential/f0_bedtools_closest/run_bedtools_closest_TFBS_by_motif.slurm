#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 48:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o out_TFBS_by_motif.log

module load gcc/7.1.0
module load bedtools/2.26.0
	 
peak_dir=../../data/cistrome/

time for condition in all overlap_motif NOT_overlap_motif
do
	file_dir=${peak_dir}/peaks_${condition}/
	outdir=data_TFBS_${condition}
	mkdir $outdir
	for ii in $file_dir/*.bed
	do
	  prename=$(basename $ii .bed)
	  echo $condition $prename
	  bedtools sort -i $ii > $outdir/${prename}.sorted.bed
	  bedtools closest \
	  -a $outdir/${prename}.sorted.bed \
	  -b $outdir/${prename}.sorted.bed -D ref -fd -io -t first > $outdir/${prename}.tsv
	  rm $outdir/${prename}.sorted.bed
	done
done
# bedtools closest -a $notch_motif -b $notch_motif -d -io -t first > $outdir/notch_motif_closest.tsv


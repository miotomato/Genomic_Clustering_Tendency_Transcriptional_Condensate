#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 48:00:00
#SBATCH -p standard
#SBATCH -A cphg_cz3d
#SBATCH -o out_run1.log

file_dir=/nv/vol190/zanglab/shared/Motif/sites/hg38_seqpos_cistrome/cutoff10k/
outdir=f1_bedtools_closest
mkdir $outdir

for ii in $file_dir/*.bed
do
  prename=$(basename $ii .txt.bed)
  echo $ii $prename
  bedtools closest -a $ii -b $ii -d -io -t first > $outdir/${prename}.tsv
done

# bedtools closest -a $notch_motif -b $notch_motif -d -io -t first > $outdir/notch_motif_closest.tsv


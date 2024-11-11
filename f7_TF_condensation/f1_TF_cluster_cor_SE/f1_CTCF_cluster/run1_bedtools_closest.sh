
outdir=f1_bedtools_closest
mkdir $outdir

ctcf_union_peak=data/union_binding_occupancy_score_GT3.bed
ctcf_motif=data/CTCF_MA0139.1.bed
myc_motif=data/MYC_MA0147.3.bed
notch_motif=data/RBPJ_MA1116.1.bed

bedtools closest -a $ctcf_union_peak -b $ctcf_union_peak -d -io -t first > $outdir/ctcf_union_peak_closest.tsv
bedtools closest -a $ctcf_motif -b $ctcf_motif -d -io -t first > $outdir/ctcf_motif_closest.tsv
bedtools closest -a $myc_motif -b $myc_motif -d -io -t first > $outdir/myc_motif_closest.tsv
bedtools closest -a $notch_motif -b $notch_motif -d -io -t first > $outdir/notch_motif_closest.tsv


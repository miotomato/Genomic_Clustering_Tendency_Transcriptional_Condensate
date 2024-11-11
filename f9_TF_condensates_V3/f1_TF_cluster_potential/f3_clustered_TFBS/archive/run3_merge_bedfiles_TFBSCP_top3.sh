outdir=f2_bedfiles_merge
mkdir $outdir
rm $outdir/*

# merge file per factor per cell type
for prename in HCT-116_SRF HCT-116_JUND HCT-116_CEBPB MCF-7_ERG MCF-7_E2F1 MCF-7_MYC
do
	echo
	ls f1_clustered_TFBS_percentile-5/bedfile/${prename}_*.merge.bed 
	echo $outdir/$prename.cat.bed
    cat f1_clustered_TFBS_percentile-5/bedfile/${prename}_*.merge.bed > $outdir/$prename.cat.bed
    bedtools sort -i $outdir/$prename.cat.bed > $outdir/$prename.sort.bed
    bedtools merge -i $outdir/$prename.sort.bed > $outdir/$prename.merge.bed
done


# merge file per cell type
for prename in HCT-116 MCF-7
do 
	echo
	ls $outdir/${prename}_*.merge.bed
	echo $outdir/$prename.cat.bed
    cat $outdir/${prename}_*.merge.bed > $outdir/$prename.cat.bed
    bedtools sort -i $outdir/$prename.cat.bed > $outdir/$prename.sort.bed
    bedtools merge -i $outdir/$prename.sort.bed > $outdir/$prename.merge.bed
done


# bedtools intersect for each cell type
ct_name=MCF-7
for prename in MCF-7_ERG MCF-7_E2F1 MCF-7_MYC
do 
    afile=$outdir/$ct_name.merge.bed
    bfile=$outdir/$prename.merge.bed
    bedtools intersect -a $afile -b $bfile -wa -c > $outdir/${ct_name}_merge_overlap_${prename}.bed
done

ct_name=HCT-116
for prename in HCT-116_SRF HCT-116_JUND HCT-116_CEBPB
do 
    afile=$outdir/$ct_name.merge.bed
    bfile=$outdir/$prename.merge.bed
    bedtools intersect -a $afile -b $bfile -wa -c > $outdir/${ct_name}_merge_overlap_${prename}.bed
done

# bedtools intersect with promoter/exon/intron 
data_dir=/nv/vol190/zanglab/zw5j/data/geneID_annotation/hg38/
for ct_name in HCT-116 MCF-7
do 
  for prename in hg38_exons hg38_introns hg38_4k_promoter_geneID 
  do 
    afile=$outdir/$ct_name.merge.bed
    bfile=$data_dir/$prename.bed
    bedtools intersect -a $afile -b $bfile -wa -c > $outdir/${ct_name}_merge_overlap_${prename}.bed
  done
done


# bedtools intersect with SE
data_dir=../../data/SE_hg38/
for ct_name in HCT-116 MCF-7
do 
    afile=$outdir/$ct_name.merge.bed
    bfile=$data_dir/$ct_name.bed
    bedtools intersect -a $afile -b $bfile -wa -c > $outdir/${ct_name}_merge_overlap_${ct_name}_SE.bed
done




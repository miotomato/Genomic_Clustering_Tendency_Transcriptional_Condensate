
mkdir SE_hg38
mapchain=/nv/vol190/zanglab/zw5j/data/liftover/hg19ToHg38.over.chain

for ii in SE_hg19/*bed
do
  outname=$(basename $ii)
  echo $outname
  liftOver $ii $mapchain SE_hg38/${outname} unMapped
done


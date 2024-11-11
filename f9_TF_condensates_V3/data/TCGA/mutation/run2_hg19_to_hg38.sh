
mkdir icgc_mutation_bed_hg38

for cancerType in BRCA CESC COAD LIHC PRAD
# for cancerType in PRAD
do

liftOver icgc_mutation_bed_hg19/${cancerType}_mutation.bed /nv/vol190/zanglab/zw5j/data/liftover/hg19ToHg38.over.chain icgc_mutation_bed_hg38/${cancerType}_mutation.bed unmapped

done


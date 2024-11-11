

for cancerType in BRCA CESC COAD LIHC PRAD

do 
  wget https://dcc.icgc.org/api/v1/download?fn=/release_28/Projects/${cancerType}-US/simple_somatic_mutation.open.${cancerType}-US.tsv.gz -O simple_somatic_mutation.open.${cancerType}-US.tsv.gz
  gzip -d simple_somatic_mutation.open.${cancerType}-US.tsv.gz

done




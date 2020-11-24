sample=$1 #The STEX sample to be tested

echo "$sample	ANN_SPECIAL" > tmp.$sample.groups.txt
cat Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.admixr.notex.groups.txt >> tmp.$sample.groups.txt

cat ../Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.DA.txt | perl ../DA2Dstatfreqs.pl tmp.$sample.groups.txt ANN_CEN.ANN_SPECIAL.DEB.groups.txt | gzip > Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.$sample.dstat.txt.gz

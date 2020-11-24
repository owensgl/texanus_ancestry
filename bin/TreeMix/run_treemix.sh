/home/owens/bin/treemix-1.13/src/treemix -i Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.0p2ld.treemix.gz -o Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.0p2ld.treemix -root OUT
/home/owens/bin/treemix-1.13/src/treemix -i Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.treemix.gz -o Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.treemix -k 1000 -root OUT

for i in `seq 10`
do
/home/owens/bin/treemix-1.13/src/treemix -i Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.0p2ld.treemix.gz -o Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.0p2ld.treemix.$i -m $i -root OUT
/home/owens/bin/treemix-1.13/src/treemix -i Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.treemix.gz -o Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.treemix.$i -k 1000 -m $i -root OUT
done

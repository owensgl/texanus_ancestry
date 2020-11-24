rep=$1

zcat Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.0p2ld.annargdeb.vcf.gz | perl /home/owens/bin/reformat/vcf2subsample.pl 0.05  | perl /home/owens/bin/reformat/vcf2structure.pl Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.admixr.annargdeb.groups.txt > Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.0p2ld.annargdeb.0p05sub.rep$rep.str

nloci=$(head -n 1 Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.0p2ld.annargdeb.0p05sub.rep$rep.str | cut -f 7- | tr '\t' '\n' | wc -l)

/home/owens/bin/console/structure -K 3 -L $nloci -N 234 -i Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.0p2ld.annargdeb.0p05sub.rep$rep.str -o Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.0p2ld.annargdeb.0p05.rep$rep.k3.str.out
/home/owens/bin/console/structure -K 5 -L $nloci -N 234 -i Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.0p2ld.annargdeb.0p05sub.rep$rep.str -o Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.0p2ld.annargdeb.0p05.rep$rep.k5.str.out 


i=$1
/home/owens/bin/admixture_linux-1.3.0/admixture Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.0p2ld.annargdeb.numchr.bed -B --cv $i | tee Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.0p2ld.annargdeb.numchr.$i.log

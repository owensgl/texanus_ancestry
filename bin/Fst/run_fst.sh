i=$1 #eg ANNCEN_ANNTEX

zcat ../Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.vcf.gz | perl /home/owens/bin/pop_gen/vcf2fst.pl Texanus.tranche90.snp.gwas.90.bi.groups.txt $i.pops.txt | gzip > Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.$i.fst.txt.gz;
zcat ../Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.vcf.gz | perl /home/owens/bin/pop_gen/vcf2fst.pl Texanus.tranche90.snp.gwas.90.bi.equalsize.groups.txt $i.pops.txt | gzip > Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.equalsize.$i.fst.txt.gz;

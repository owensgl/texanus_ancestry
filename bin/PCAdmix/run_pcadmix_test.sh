chr=$1
##/home/owens/bin/PCAdmix3_linux -anc Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.beagle.ARG.$chr.haplotypes.txt Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.beagle.ANN_CEN.$chr.haplotypes.txt  -adm Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.beagle.ANN_STEX.$chr.haplotypes.txt -map Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.beagle.$chr.physmap.txt -rho Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.beagle.$chr.genmap.txt -o Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.beagle.ANN_STEX.ARG-ANN.$chr -lab ARG ANN TEX -r2 0.8 -w 100 -bed $chr
#/home/owens/bin/PCAdmix3_linux -anc Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.beagle.DEB.$chr.haplotypes.txt Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.beagle.ANN_CEN.$chr.haplotypes.txt  -adm Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.beagle.ANN_STEX.$chr.haplotypes.txt -map Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.beagle.$chr.physmap.txt -rho Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.beagle.$chr.genmap.txt -o Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.beagle.ANN_STEX.DEB-ANN.$chr -lab DEB ANN TEX -r2 0.8 -w 100 -bed $chr
/home/owens/bin/PCAdmix3_linux -anc Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.beagle.ARG.$chr.haplotypes.txt Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.beagle.DEB.$chr.haplotypes.txt Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.beagle.ANN_CEN.$chr.haplotypes.txt  -adm Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.beagle.ANN_STEX.$chr.haplotypes.txt -map Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.beagle.$chr.physmap.txt -rho Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.beagle.$chr.genmap.txt -o Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.beagle.ANN_STEX.ARG-DEB-ANN.$chr -lab ARG DEB ANN TEX -r2 0.8 -w 100 -bed $chr
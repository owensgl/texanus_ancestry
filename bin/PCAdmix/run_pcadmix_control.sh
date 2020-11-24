chr=$1
/home/owens/bin/PCAdmix3_linux -anc tmp.ARG.$chr.haplotypes.txt tmp.DEB.$chr.haplotypes.txt tmp.ANN_CEN.$chr.haplotypes.txt  -adm tmp.CONTROL.$chr.haplotypes.txt -map tmp.$chr.physmap.txt -rho tmp.$chr.genmap.txt -o tmp.CONTROL.ARG-DEB-ANN.$chr -lab ARG DEB ANN CONTROL -r2 0.8 -w 100 -bed $chr

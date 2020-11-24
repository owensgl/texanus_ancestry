mig_file="tmp.arg.mig"
for i in `seq 0.01 0.01 0.4`
do
echo "ARG ANN_STEX $i" > $mig_file

/home/owens/bin/treemix-1.13/src/treemix -i Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.0p2ld.treemix.gz -o arg_introgression_test/Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.0p2ld.treemix.$i -root OUT  -cor_mig $mig_file
/home/owens/bin/treemix-1.13/src/treemix -i Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.treemix.gz -o arg_introgression_test/Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.treemix.$i -k 1000 -root OUT  -cor_mig $mig_file
done

for i in `seq 0.001 0.001 0.009`
do
echo "ARG ANN_STEX $i" > $mig_file

/home/owens/bin/treemix-1.13/src/treemix -i Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.0p2ld.treemix.gz -o arg_introgression_test/Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.0p2ld.treemix.$i -root OUT  -cor_mig $mig_file
/home/owens/bin/treemix-1.13/src/treemix -i Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.treemix.gz -o arg_introgression_test/Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.treemix.$i -k 1000 -root OUT  -cor_mig $mig_file
done


##perl Fst.hap.pl 0.1 spring/chr1_partA.snp.vcf.gz  winter/chr1_partA.snp.vcf.gz chr1_partA.gz ;perl combine.pl chr1_partA.gz chr1_partA.txt ###Other chromosome is same.
less chr.len |while read i st end;do mkdir $i; echo "perl Fst.hap.pl 0.1 spring/$i.snp.vcf.gz winter/$i.snp.vcf.gz $i/$i.gz;perl combine.pl  $i/$i.gz $i/$i.txt">>ww.sh;done
less chr.len |while read chr st end;do echo "perl windows_value_cal.Fst.pl $end 10000 2000 0.1 ECR/split/$chr.gz $chr/$chr.gz $chr/winadow_Fst.$chr.xls" >>ww2.sh;done

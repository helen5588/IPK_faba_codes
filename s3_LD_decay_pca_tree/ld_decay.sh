#https://github.com/BGI-shenzhen/PopLDdecay.git
cat chr.len |while read i st end;do mkdir $i;echo "PopLDdecay -InVCF $i.snp.vcf.gz -OutType  1 -MaxDist 10000  -OutStat $i.stat ">>ld.sh;done

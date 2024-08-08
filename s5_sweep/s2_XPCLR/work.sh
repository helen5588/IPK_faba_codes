module load xpclr
module load vcftools

for k in {1..7};do echo "vcftools --gzvcf $k.gz --plink --out $k.MT;awk 'BEGIN{OFS=" "} {print 1,".",$4/2000000,$4}' chr{k}.MT.map > chr{k}.MT.map.distance" >>ww2.sh;done

for k in {1..7}; do echo "/opt/Bio/xpclr/1.1.2/bin/xpclr --out ./$k.out --format vcf  --input $k.gz  --samplesA spring.lst   --samplesB winter.lst --map $k.MT.map.distance  --chr ${k} --gdistkey None --phased  --size 50000  --step 10000;">>xpclr.sh;done
#module load xpclr;/opt/Bio/xpclr/1.1.2/bin/xpclr --out 1/1.out --format vcf --input 1.gz  --samplesA spring.lst --samplesB winter.lst --map 1.MT.map.distance  --chr 1 --gdistkey None --phased  --size 50000  --step 10000;
for k in  {1..7}
do awk  '{print $2,$3,$4,$12,$13}'   ${k} > ${k}.chart.xpclr.50kb.windows;
 cat ./*.chart.xpclr.50kb.windows > all.xpclr.50kb.windows;
done


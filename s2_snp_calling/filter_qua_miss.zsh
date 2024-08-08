#!/bin/zsh

vcf=$1
shift
vcfout=$1
shift
output=$1
shift
maf=$1


perl filter_qua_miss_hete_nopoly.pl 0 0.1 0.1 $vcf $output/$vcfout.gz 2>>$output/filter.err
/opt/Bio/vcftools/0.1.14/bin/vcftools --gzvcf $output/$vcfout.gz --maf $maf --recode --stdout | /opt/Bio/bcftools/1.15.1/bin/bgzip >$output/$vcfout.$maf.vcf.gz 2>>$output/maf.err


echo "This work is completed"

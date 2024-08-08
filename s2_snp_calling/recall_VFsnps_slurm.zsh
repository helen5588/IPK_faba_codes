#!/bin/zsh

bgzip='/opt/Bio/bcftools/1.10.2/bin/bgzip'
tabix='/opt/Bio/bcftools/1.10.2/bin/tabix'
bcftools='/opt/Bio/bcftools/1.10.2/bin/bcftools'

bed=$1
shift
vcf=$1
shift
out=$1


{
 date 1>&2
 {
  pigz -cd $vcf | ./recall_all_vcf.awk -v header=1 | $bgzip
  $bcftools view $vcf | ./recall_all_vcf.awk \
    -v header=0 \
    -v bed=$bed\
    -v dphom=2 \
    -v dphet=4 \
    -v minqual=40 \
    -v mindp=0 \
    -v dphommax=50 \
    -v minhomn=0 \
    -v minhomp=0 \
    -v tol=0.2 \
    -v minmaf=0 \
    -v minpresent=0.01 | $bgzip
 } > $out.vcf.gz && \
  /opt/Bio/bcftools/1.5/bin/tabix -Cp vcf $out.vcf.gz 
 date 1>&2 
} 2> ${out:r:r}.err




#!/bin/zsh

#SBATCH -p ipk
#SBATCH --output=bcf_snpcall.txt
#SBATCH --error=bcf_snpcall.err

bcftools='/opt/Bio/bcftools/1.8/bin/bcftools'
bgzip='/opt/Bio/bcftools/1.8/bin/bgzip'
tabix='/opt/Bio/bcftools/1.8/bin/tabix'

ref=$1
shift
bamlist=$1
shift
prefix=$1
shift
threads=$1
shift
part=$1

base=$prefix/${prefix}_${part:t:r}

ulimit -HSn 102400

{
 $bcftools mpileup --threads $threads -d 10000000 -q 20 -Q 20 --excl-flags 3332 -a DP,DV,AD -b $bamlist -R $part -f $ref \
  | $bcftools call --threads $threads  -mv -f GQ - | $bgzip -@ 2 > $base.vcf.gz

 echo $pipestatus | tee ${base}_pipestatus.txt  | tr ' ' '\n' | grep -q '^[^0$]' || $tabix -C -p vcf $base.vcf.gz 
} 2> $base.err

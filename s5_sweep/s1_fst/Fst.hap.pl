=pod

1. the perl is to calculate the Fst per SNPs between two pompared populations (Bhatia G, Patterson N J, Sankararaman S, et al. Estimating and interpreting FST: the impact of rare variants. Genome research 2013)
2. the diploid SNPs data will be convert to haploid SNPs data in this perl, by coding all heterozygous genotype calls as missing data.
3. the converted haploid SNPs data will be filtered by an missing rate
4. the perl is only suitable for : (i) VCFv4.1 or other similar format (GT:AD:DP:GQ:PGT:PID:PL or GT:AD:DP:GQ:PL) (ii) bi allele variant
5. Usage: perl ./vcf_Fst.hap.pl maximum_missing_rate pop1.vcf.gz pop2.vcf.gz output 

maximum_missing_rate FLOAT the maximum proportion of the missing genotype calls after conversion

pop1.vcf.gz FILE the compressed input VCF file of population 1

pop2.vcf.gz FILE the compressed input VCF file of population 2

output FILE the value of nucleotide diversity per SNPs position

note: the output files contant three columns: chromosome, SNPs position, Fst

'miss' means 
'nan' means fm==0

=cut

#!/usr/bin/perl -w
use strict;
use warnings;

die '@ARGV is required' if @ARGV != 4;

my $max_miss_rate=shift;
my $vcf1=shift;
my $vcf2=shift;
my $output=shift;

open IN1,"gzip -dc $vcf1|grep -v \"^#\"|" or die $!;
open IN2,"gzip -dc $vcf2|grep -v \"^#\"|" or die $!;
open OUT,"| gzip >$output" or die $!;
print OUT "chr\tpos\tvalue\n";
while (my $ln1=<IN1>)
{
    my $ln2=<IN2>;
    chomp ($ln1,$ln2);

    my @ln1=split /\s+/,$ln1;
    my @ln2=split /\s+/,$ln2;
    die $! if $ln1[1] ne $ln2[1]; 
   
    my $miss1=0;
    my $geno1_0=0;
    for (my $i=9;$i<@ln1;$i++)
    {
        if ($ln1[$i]=~/^([^:])\/([^:]):/)
        {
            if ($1 eq '0' && $2 eq '0'){$geno1_0++}
            elsif ($1 eq '1' && $2 eq '1'){next}
            else{$miss1++}
        }
    }
    if ($miss1/(@ln1-9) > $max_miss_rate)
    {
        print OUT "$ln1[0]\t$ln1[1]\tmiss\n";
        next;
    }   
    my $sample_size1=@ln1-9-$miss1;
    my $fre1=$geno1_0/$sample_size1;

    my $miss2=0;
    my $geno2_0=0;
    for (my $i=9;$i<@ln2;$i++)
    {
        if ($ln2[$i]=~/^([^:])\/([^:]):/)
        {   
            if ($1 eq '0' && $2 eq '0'){$geno2_0++}
            elsif ($1 eq '1' && $2 eq '1'){next}
            else{$miss2++}
        }
    }
    if ($miss2/(@ln2-9) > $max_miss_rate)
    {
        print OUT "$ln1[0]\t$ln1[1]\tmiss\n";
        next;
    }
    my $sample_size2=@ln2-9-$miss2;
    my $fre2=$geno2_0/$sample_size2;    
    
    my $Fst=&Fst($fre1,$fre2,$sample_size1,$sample_size2);
    if ($Fst eq 'nan'){print OUT "$ln1[0]\t$ln1[1]\tnan\n"}
    else {print OUT "$ln1[0]\t$ln1[1]\t$Fst\n"}
}
close IN1;
close IN2;
close OUT;

sub Fst
{
    my $p1=$_[0];
    my $p2=$_[1];
    my $n1=$_[2];
    my $n2=$_[3];

    my $fz1=($p1-$p2)*($p1-$p2);
    my $fz2_1=$p1*(1-$p1)/($n1-1);
    my $fz2_2=$p2*(1-$p2)/($n2-1);
    my $fm=($p1*(1-$p2))+($p2*(1-$p1));

    my $value;
    if ($fm==0){$value='nan'}
    else {$value=($fz1-$fz2_1-$fz2_2)/$fm}
    return $value;
}    

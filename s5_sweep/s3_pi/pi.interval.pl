=pod

1. the perl is to calculate the nucleotide diversity (pi) for each SNPs site (Tajima F. Evolutionary relationship of DNA sequences in finite populations. Genetics 1983)
2. non-polymorphic variants will be pass 
3. Usage: perl ./vcf_pi.hap.pl maximum_missing_rate input output

maximum_missing_rate FLOAT the maximum proportion of the missing genotype calls after conversion

input FILE     the compressed input VCF file   eg. in.vcf.gz

output FILE    the value of nucleotide diversity per SNPs position

note: the output files contant three columns: chromosome, SNPs position, pi

=cut

#!/usr/bin/perl -w
use strict;
use warnings;

die '@ARGV is required' if @ARGV != 4;

my $max_miss_rate=shift;
my $interval=shift; # 1000,2000
my $vcf=shift;
my $output=shift;

my @region=split /,/,$interval;

open IN,"gzip -dc $vcf | grep -v \"^#\"|" or die $!;
open OUT,"| gzip >$output" or die $!;
print OUT "#chr\tpos\tvalue\n";
while (<IN>)
{
    chomp;
    my @a=split;
    next if $a[1]<$region[0];
    last if $a[1]>$region[1];
    my (@geno,%geno);
    my $miss=0;
    for (my $i=9;$i<@a;$i++)
    {        
        if ($a[$i]=~/^([^:])\/([^:]):/)
        {
            if ($1 eq '0' && $2 eq '0'){push @geno,'a','a';$geno{"0"}=1}
            elsif ($1 eq '1' && $2 eq '1'){push @geno,'b','b';$geno{"1"}=1}
            elsif ($1 eq '0' && $2 eq '1'){push @geno,'a','b';$geno{"0"}=1;$geno{"1"}=1}
            else {$miss++}
        }
        else{die $!}
    }
    next if $miss/(@a-9)>$max_miss_rate;

    my $pi;
    if (not exists $geno{"0"} or not exists $geno{"1"}){$pi=0}
    else{$pi=&pi(@geno)}
    print OUT "$a[0]\t$a[1]\t$pi\n";
}
close IN;
close OUT;

sub pi
{
    my @a=@_;
    my $diff=0;
    my $time=0;

    for (my $i=0;$i<@a;$i++)
    {
        for (my $m=0;$m<@a;$m++)
        {
            if ($m>$i)
            {
                $diff++ if $a[$i] ne $a[$m];
                $time++;
            }
        }
    }
    my $value=$diff/$time;
    return $value;
}    

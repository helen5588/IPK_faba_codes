# only for bi-allel

#!/usr/bin/perl -w
use strict;
use warnings;

die '@ARGV is required' if @ARGV != 5;

my $min_qua=shift;
my $max_missing=shift;
my $max_hete_pro=shift;
my $vcf=shift;
my $output=shift;

open IN,"gzip -dc $vcf |" or die $!;
open OUT,"| /opt/Bio/bcftools/1.15.1/bin/bgzip > $output";
while (<IN>)
{
    chomp;
    if ($_=~/^#/){print OUT "$_\n"}
    else
    {
        my @a=split;
        next if $a[4]=~/,/; # only bi-allele 
        next if $a[5]<$min_qua; # filter by quality
        if ($a[4] eq '.'){print "$_\n"} 
        else
        {
            my %x;
            my $miss=0;
            my $hete=0;
            my $hemo=0;    
            for (my $i=9;$i<@a;$i++)
            {
                if ($a[$i]=~/^([^:])[\/\|]([^:]):/)
                {
                    if ($1 eq '.'){$miss++}
                    elsif ($1 eq '0' && $2 eq '1'){$hete++}
                    else{$hemo++}        
                    my $key=$1.$2;
                    $x{$key}=1;
                }
	        else{die $!}
            }
            if (exists $x{"00"} && exists $x{"11"}) # filter no poly
            {
                next if $miss/(@a-9)>$max_missing; # filter missing
                next if $hete/($hete+$hemo)>$max_hete_pro; # filter_hete
                print OUT "$_\n"
            }
        }
    }
}
close IN;
close OUT;

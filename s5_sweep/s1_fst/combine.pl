#!/usr/bin/perl -w
use strict;
use warnings;

die '@ARGV is required' if @ARGV != 2;

my $fst=shift;
my $output=shift;

my ($sum,$site);
open IN,"gzip -dc $fst|" or die $!;
<IN>;
while (<IN>)
{
    chomp;
    my @a=split;
    next if $a[2]=~/miss/;
    $a[2]=0 if $a[2]=~/nan/;
    $a[2]=0 if $a[2]<0;
    $sum+=$a[2];
    $site++;
}
close IN;

open OUT,">$output" or die $!;
print OUT "$sum\t$site\n";
close OUT;

=pod

1. the perl is to calculate the windows value of pi, theta-w or He based on their single SNP value.
2. the perl is only fit for a single chromosome input file 
3. Usage: perl ./windows_value_cal.pi_theta-w_He.pl chr_length window_size step_size minimum_effective_covered_proportion effective_covered_file single_value output

chr_length     INT   the length of the chromosome 
window_size    INT   the windows size  eg. 10000
step_size      INT   the step size  eg. 2000
minimum_effective_covered_proportion  FLOAT  the minimum effective covered region proportion in a window eg. 0.2
effective_covered_file FILE  the compressed file of effective covered region, contain three columns: chromosome start end
single_value   FILE  the single value file of pi, theta-w or He with a header, contain three columns: chromosome position value
output         FILE  the output file

=cut

#!/usr/bin/perl -w
use strict;
use warnings;

die '@ARGV is required' if @ARGV != 7;

my $chr_len=shift;
my $window=shift;  # 10
my $step=shift; # 2
my $min_ER_pro=shift; # 0.5 
my $ER_region=shift; # effective_covered_file
my $site_file=shift;
my $output=shift;

my $step_num;
if ($step !=0){$step_num=$window/$step} # 5
else{$step_num=1}

my %ER;
open IN,"gzip -dc $ER_region|" or die $!;
while (<IN>)
{
    chomp;
    my @a=split;
    for my $pos($a[1]..$a[2])
    {
        for my $n(0..$step_num-1)
        {
            next if $pos-1-$step*$n<0;
            my $num=int (($pos-1-$step*$n)/$window); #(25-1-2*1)/10
            my $st=$num*$window+1+$step*$n; #print "$a[1]\t$n\t$num\t$st\n";
            $ER{$st}++;
        }
    }   
}
close IN;

my (%win,%snp_num,$chr);
if ($site_file=~/gz$/){open IN,"gzip -dc $site_file |" or die $!}
else{open IN,$site_file or die $!}
<IN>;
while (<IN>)
{
    chomp;
    my @a=split;
    $chr=$a[0];
    next unless $a[-1]=~/^[0-9]/;

    for my $n(0..$step_num-1) # (0..4)
    {
        next if $a[1]-1-$step*$n<0; 
        my $num=int (($a[1]-1-$step*$n)/$window); #(25-1-2*1)/10
        my $st=$num*$window+1+$step*$n; #print "$a[1]\t$n\t$num\t$st\n";
        $win{$st}+=$a[-1];
        $snp_num{$st}++;
    }               
}
close IN;

open OUT,">$output";
print OUT "#chr\tstart\tend\teff_cov_reg\tsnp_num\tvalue\n";
my $window_num=int($chr_len/$window);
for (0..$window_num-2)
{
    for my $n(0..$step_num-1)
    {
        my $st=$_*$window+1+$step*$n;
        my $end=$st+$window-1;
        if (exists $ER{$st})
        {
            next if $ER{$st}/$window<$min_ER_pro;
            my ($win_value,$snp_num);
            if (exists $win{$st}){$win_value=$win{$st}/$ER{$st};$snp_num=$snp_num{$st}}
            else{$win_value=0;$snp_num=0}
            print OUT "$chr\t$st\t$end\t$ER{$st}\t$snp_num\t$win_value\n";
        } 
    }   
}
close OUT;

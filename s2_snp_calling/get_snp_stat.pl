#!/usr/bin/perl -w
use strict;
unless(@ARGV)
{
	die "<input file> <window:10000>\n";
}
my $ltr=shift;
my $window=shift;
if($ltr=~/\.gz$/){
   open IN,"pigz -dc $ltr|" || die $!;
  }else{
     open IN,$ltr || die $!;
    }
my %save;
my %zz;
while(<IN>)
{
	next if(/^#/);
	chomp;
	my @line=split(/\s+/,$_);
	#$line[0]=~s/chr(0)*/hs/;
	$line[0]="$line[0]";

	#my @z=split //,$line[0];
	#$line[0]="$z[1]$z[2]$z[0]";
	my $pos1=$line[1];
	my $pos2=$line[1];
	my $len=$pos2-$pos1+1;
#	for(my $i=0;$i<$len;$i++){
#		my $pos=$pos1+$i;
	$save{$line[0]}->{int($pos1/$window)}+=1;##250000

#	if(not defined $save{$line[0]}->{int($pos1/$window)}){$save{$line[0]}->{int($pos1/$window)}=0;}
#	elsif($save{$line[0]}->{int($pos1/$window)} < $pos2){$save{$line[0]}->{int($pos1/$window)}=$pos2;}
	$zz{$line[0]}->{int($pos1/$window)}++;
#	}
}
print "#chr\tStart\tEnd\tGene\n";
foreach my $i(sort{$a cmp $b}keys %save){
	foreach my $j(sort{$a <=> $b}keys %{$save{$i}})	{
		my $str=$j*$window;
		my $end=($j+1)*$window;

#		if($save{$i}->{$j}>100){$save{$i}->{$j}=1;}##50
#		else{$save{$i}->{$j}=$save{$i}->{$j}/100;}
#		print  "$i $str $end $save{$i}->{$j}\n";

		$save{$i}->{$j}=$save{$i}->{$j}/$window;
#		$save{$i}->{$j}=$save{$i}->{$j};
		print  "$i\t$str\t$end\t$save{$i}->{$j}\n";
	}
}


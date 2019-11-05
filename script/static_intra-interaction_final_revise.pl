#!/usr/bin/perl
#static intra-interaction at a distance from validPairs files#
use strict;
use warnings;
my %cal = ();
my $count;
my $fc;
my $len;
my $dis;
my $key;
my $in;
my $out;
print "input allvalidPairs file: $ARGV[0]\ninput defined distance: $ARGV[1]\ninput defined FC of distance:$ARGV[2]\n";
open (ORA, "<$ARGV[0]") or die "error reading $ARGV[0] for reading\n";
$len = $ARGV[1];
$fc = $ARGV[2];
$count = 0;	
$in = 0;
$out = 0;
while (<ORA>)
{
	chop($_) while ($_ =~/[\r\n]$/);
	my @arr = split(/\t/,$_);
	if ($arr[1] eq $arr[4])
	{
		$count = $count+1;
		$dis = abs($arr[5]-$arr[2]);
		for (my $i=1;$i<=$fc;$i++)
		{
			$key = $i*$len;
			if ($dis<=$key)
			{
				$cal{$key} =$cal{$key} +1;
			}		
		}
		if ($dis<=2000000)
		{
			$in = $in+1;
		}
		else
		{
			$out = $out+1;
		}
	}
}
for my $base (keys(%cal))
{
	my $ratio1 = $cal{$base}/$count;
	my $ratio2 = $in/$count;
	my $ratio3 = $out/$count;
	print "$base\t$cal{$base}\t$in\t$out\t$count\t$ratio1\t$ratio2\t$ratio3\n";
}
close ORA;

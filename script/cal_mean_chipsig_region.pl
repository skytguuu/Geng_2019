#!/usr/bin/perl
#unique boundary identification 
$my_region = @ARGV[0];#输入region文件
$my_chip = @ARGV[1];
$my_output = @ARGV[2];#输出output结果
open (REG, "$my_region");
open (CHIP, "$my_chip");
open (R,">$my_output");
%hash = ();
while (my $line = <CHIP>)
{
	chop($line) while ($line =~/[\r\n]$/);
	@arr = split(/\t/,$line);
	$key = "$arr[0]\t$arr[1]\t$arr[2]";
	if ($arr[3]<0)
	{
		$hash{$key}=0;
	}
	else
	{
		$hash{$key} = $arr[3];
	}
}
close CHIP;
while (my $line = <REG>)
{
	chop($line) while ($line =~/[\r\n]$/);
	@pos = split(/\t/,$line);
	$len = ($pos[2]-$pos[1]+1)/1000;
	$sig = 0;
	foreach $base (keys(%hash))
	{
		@loc = split(/\t/,$base);
		if ($loc[0] eq $pos[0])
		{
			if (($pos[1]<=$loc[1]) and ($pos[2]>=$loc[2]))
			{
				$sig = $sig + $hash{$base};
			}
		}
	}
	$mean = $sig/$len;
	print R "$pos[0]\t$pos[1]\t$pos[2]\t$mean\n";
}
close REG;
close R;
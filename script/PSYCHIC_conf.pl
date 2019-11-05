#!/usr/bin/perl
$dir = "/home/dell/mnt1/HiC_reference/development/HiTC/20kb/PSY/Input";
opendir (DIR, $dir ) || die "Error in opening dir $dir\n";
@dir_list = ();
while($filename = readdir(DIR))
{
	push @dir_list,$filename;
}
foreach $base (@dir_list)
{
	$dirname = "$dir/$base";
	opendir (D, $dirname ) || die "Error in opening dir $dirname\n";
	while($file = readdir(D))
	{
		if ($file =~ /chr/){
		@id = split(/_/,$file);
		$name = uc($base.$id[0]);
		$out = $id[0]."_".$id[1].".conf";
		open ($name,">$out") || die "Error in creat $name\t$out\n";
		print $name "[DefaultSection]
res: 20000
win: 2000000
chrname: $id[0]
chrsize: /home/dell/softwares/tools/Hi-C/PSYCHIC-master/examples/mm9.chrom.sizes
output_prefix: $out
output_dir: /home/dell/mnt1/HiC_reference/development/HiTC/20kb/PSY/Output/$base/
input_matrix: $dirname/$file
genes_file: /home/dell/softwares/tools/Hi-C/PSYCHIC-master/examples/mm9.genes2.bed\n";}
	close $name;
	}
}

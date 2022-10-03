#!/usr/bin/perl

#please run this job in a new clear fold
die "usage: perl run_mark_duplicates.pl <Data_Index> <fastq-dir> <START> <END> \n" if @ARGV!= 4;

$index_table = $ARGV[0];
$root = $ARGV[1];


open INDEX, $index_table;
$i = 0;
$id = 0;
while($line = <INDEX>){
	chomp($line);
	my @temp=split('\t',$line);
	$sample[$i][0] = $temp[1];
	$sample[$i][1] = $temp[2];
	$i++;
}
close INDEX;

#`mkdir Dir_tmp_1`;
#chdir "./Dir_tmp_1";
#
$start = $ARGV[2];
$end = $ARGV[3];

for($j=$start;$j<=$end;$j++){          #run the first 100 samples as a try
	`sh ./do_mark_duplicates.sh $sample[$j][1]`;
	$completed_No = $j;
        print "MarkingDup completed!!! No.$completed_No : $sample[$j][0]\t$sample[$j][1]\n\n";
}
#chdir "../";



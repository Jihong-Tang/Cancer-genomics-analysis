#!/usr/bin/perl

#please run this job in a new clear fold
die "usage: perl maping.pl <Data_Index> <START> <END> \n" if @ARGV!= 3;

$index_table = $ARGV[0];


open INDEX, $index_table;
$i = 0;
$id = 0;
while($line = <INDEX>){
	chomp($line);
	my @temp=split('\t',$line);
	$sample[$i][0] = $temp[1];
	$i++;
}
close INDEX;

#`mkdir Dir_tmp_1`;
#chdir "./Dir_tmp_1";
#
$start = $ARGV[1];
$end = $ARGV[2];

for($j=$start;$j<=$end;$j++){
	`sh do_savi_WES_2sample_GRCh38.p7.RefSeq.sh  $sample[$j][0] ../all_mapping_hg38/$sample[$j][0]_B_WES_hg38.sorted.MD.bam ../all_mapping_hg38/$sample[$j][0]_T_WES_hg38.sorted.MD.bam`;
	$completed_No = $j;
        print "bwa-MEM completed!!! No.$completed_No : $sample[$j][0]\n";
}
#chdir "../";



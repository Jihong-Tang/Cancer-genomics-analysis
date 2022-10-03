#!/usr/bin/bash
A='WangLab'
Sample_ID=$1
FAQ1=$2
FAQ2=$3
result=$4


echo [Directory] `pwd`
echo [Machine] `uname -n`
echo [Start] `date`
echo [args] $*
time1=$( date "+%s" )


#line_q1=`zcat $FAQ1 | wc -l`
#line_q2=`zcat $FAQ2 | wc -l`

bwa mem -t 20 -T 0 -R "@RG\tCN:$A\tID:$Sample_ID\tSM:$Sample_ID" /scratch/PI/jgwang/jtangbd/reference/TCGA_GRCh38_bwa_ref/GRCh38.d1.vd1.fa $FAQ1 $FAQ2 2> ./tmp/tmp.$result | samtools view -Sbh1 -@ 8 - -o $result\.bam

#### mapping duplicates
###counting the num of bam reads

samtools sort -@ 8 -O bam -o $result.sorted.bam -T the_temp_$result $result\.bam
samtools index $result\.sorted.bam 
rm $result\.bam
#rm $FAQ1
#rm $FAQ2
#samtools idxstats result_test.sorted.bam |head 
#samtools view -H /home/qmu/projects/AsianpGBM/TCGAGBM/C484.TCGA-06-6699-10A-01D-1845-08.2_gdc_realn.bam |grep "^@PG"
#samtools flagstat result_test.sorted.bam
#samtools view -f 4 result_test.sorted.bam chr7 | wc -l
#samtools mpileup ../raw_data/results_B_DNA/IGCT_S002_B.sorted.bam -r chr10:63207680 -f /home/dsongad/software/my_lib/STAR_index_hg38_RefSeq_hg38/hg38.fa |head -n20

time2=$( date "+%s" )
echo [deltat] $(( $time2 - $time1 ))
echo [End]

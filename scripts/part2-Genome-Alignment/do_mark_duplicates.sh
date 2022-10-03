#!/usr/bin/bash

#set -euxo pipefail
#export LD_LIBRARY_PATH=/home/qmu/tools/lib:/opt/intel/advisor_xe_2016.1.40.463413/lib64:/home/qmu/tools/xz-5.2.3/tmp/lib 
#picard=/home/share/jgwang/softwares/picard/build/libs/picard.jar
picard=/scratch/PI/jgwang/dsongad/software/Picard/picard.jar

NAME=$1

mkdir md_tmp/tmp_$NAME

java -Djava.io.tmpdir=md_tmp/tmp_$NAME \
	-Xmx120g -XX:+UseParallelGC -XX:ParallelGCThreads=4 \
	 -jar ${picard} MarkDuplicates \
	INPUT=${NAME}.sorted.bam \
	OUTPUT=${NAME}.sorted.MD.bam \
	METRICS_FILE=./md_metrix/${NAME}.metrics.txt \
	ASSUME_SORTED=true \
	REMOVE_DUPLICATES=true \
	VALIDATION_STRINGENCY=LENIENT

rm ${NAME}.sorted.bam
rm ${NAME}.sorted.bam.bai
samtools index ${NAME}.sorted.MD.bam
#rm ${NAME}.sorted.MD.bam.txt

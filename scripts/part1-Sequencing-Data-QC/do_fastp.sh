#!/usr/bin/bash

rawName=$1
cleanName=$2_fp

#mkdir -p qc_res

fastp -w 2 -i ./$rawName\.R1.fq.gz -o ./$cleanName\.R1.fq.gz \
 -I ./$rawName\.R2.fq.gz -O ./$cleanName\.R2.fq.gz \
 -q 15 -u 40 --length_required 45 \
 --detect_adapter_for_pe \
 -h ./qc_res/fastp_$cleanName\.html -j ./qc_res/fastp_$cleanName\.json -R ./qc_res/report_$cleanName\.txt

rm ./$rawName\.R1.fq.gz ./$rawName\.R2.fq.gz

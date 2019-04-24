#!/bin/bash

st=$1
ed=$2

for i in $(seq $st $ed) :
do
	target=`sed -n $i"p" /home/haeun/DATA/splicing/src/QoRT_bam_list.txt`
	out=`echo $target | cut -f 1 -d "." | cut -f 10 -d "/"`
	echo "Filtering "$out
	samtools view -b -F 4 -F 8 $target > "filtered_"$out".bam"	# filter out unmapped read
	echo "Complete filtering "$out
	mkdir /home/haeun/DATA/splicing/Analysis/Quant/QoRT/$out
	java -jar /home/haeun/DATA/tools/QoRTs.jar QC "filtered_"$out".bam" ~/DATA/splicing/data/GENCODE/gencode.v29.annotation.gtf.gz /home/haeun/DATA/splicing/Analysis/Quant/QoRT/$out
	rm "filtered_"$out".bam"
done


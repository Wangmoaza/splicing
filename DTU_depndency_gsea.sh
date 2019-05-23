#!/bin/bash

DIR=/home/haeun/DATA/splicing/Analysis/transcript_usage/CCLE-OV/dependency
doit(){
	DIR=/home/haeun/DATA/splicing/Analysis/transcript_usage/CCLE-OV/dependency
	while read line
	do
		geneset=$(echo $line | cut -d, -f1)
		gmt=$(echo $line | cut -d, -f2)
		echo $gmt
		mkdir -p $DIR/gsea/${2}_${geneset}/
		gseapy prerank -r $1 -g $gmt -o "$DIR/gsea/${2}_${geneset}"
	done < $DIR/genset_files2.txt
			
}

export -f doit

#for i in $(find $DIR -type f | egrep "rnk$")
#do
#	echo $i
#	doit $i
#done
#doit /home/haeun/DATA/splicing/Analysis/transcript_usage/CCLE-OV/dependency/sig3_abs_d2_corr.rnk sig3_abs
find $DIR -type f | egrep "rnk$" | parallel -j2 doit {} {/.}

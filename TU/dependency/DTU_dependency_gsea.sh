#!/bin/bash

DIR=/home/haeun/DATA/splicing/Analysis/transcript_usage/PDC/dependency
doit(){
	DIR=/home/haeun/DATA/splicing/Analysis/transcript_usage/PDC/dependency
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
find $DIR -type f | egrep ".tsv$" | parallel -j4 'tail -n+2 {} | cut -f1,2 > {.}.rnk'
find $DIR -type f | egrep "rnk$" | grep -v "not" | parallel -j2 doit {} {/.}
#find $DIR -type f | grep "minor_d2_corr" | egrep ".rnk$" | parallel -j2 doit {} {/.}

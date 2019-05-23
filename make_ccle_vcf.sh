#!/bin/bash

# Haeun Hwangbo 2019-05-03

doawk(){
	INPUT=$1
	gawk -F "\t" 'BEGIN {OFS="\t"} match($4, /[0-9]+/) && ($9 == "SNP") {if($12~"rs") {print $4, $5, $12, $10, $11, ".", ".", "."} else {print $4, $5, ".", $10, $11, ".", ".", "."}}' "${INPUT}.txt" | sort -k1,2 -n > "${INPUT}_snp.tmp"
	cat /home/haeun/DATA/splicing/data/Cell_lines/CCLE/variant/CCLE-OV/vcf_header ${INPUT}_snp.tmp > ${INPUT}_snp.vcf
	rm ${INPUT}_snp.tmp
}

export -f doawk
parallel -j3 doawk {.} ::: /home/haeun/DATA/splicing/data/Cell_lines/CCLE/variant/CCLE-OV/*.txt

#!/bin/bash

TOOL_DIR="/home/haeun/Tools/gffcompare-0.11.2.Linux_x86_64"
REF="/home/haeun/DATA/splicing/data/GENCODE/gencode.v29.annotation.gtf"
OUT_DIR="/home/haeun/DATA/splicing/Analysis/Quant/gffcompare"

$TOOL_DIR/gffcompare -i gffcompare_file.txt \
	-r $REF -o $OUT_DIR/stringtie.BRCA_751.gencode.v29.gffcompare

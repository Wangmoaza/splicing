#!/bin/bash

TOOL_DIR="/home/haeun/DATA/tools/stringtie-1.3.5.Linux_x86_64"
INPUT_DIR="/home/haeun/DATA1/splicing/TCGA-OV/bamfiles"
OUTPUT_DIR="/home/haeun/DATA/splicing/Analysis/Quant/TCGA-OV/assembled"
REFERENCE="/home/haeun/DATA/splicing/data/GENCODE/gencode.v29.annotation.gtf"

while read i
do
        echo $i
	$TOOL_DIR/stringtie $INPUT_DIR/$i -p 4 -G $REFERENCE -o $OUTPUT_DIR/${i}.gtf
done < /home/haeun/DATA/splicing/data/TCGA/TCGA-OV/bamfile_list_killed.txt

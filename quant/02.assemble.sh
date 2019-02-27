#!/bin/bash

TOOL_DIR="/home/haeun/2TB_disk/tools/stringtie-1.3.5.Linux_x86_64"
INPUT_DIR="/home/haeun/2TB_disk/splicing/data/TCGA/TCGA-BRCA/bamfiles"
OUTPUT_DIR="/home/haeun/2TB_disk/splicing/data/TCGA/TCGA-BRCA/assembled"
REFERENCE="/home/haeun/2TB_disk/splicing/data/GENCODE/gencode.v29.annotation.gtf"

while read i
do
        echo $i
	$TOOL_DIR/stringtie $INPUT_DIR/$i -p 4 -G $REFERENCE -o $OUTPUT_DIR/${i}.gtf
done < /home/haeun/2TB_disk/splicing/data/TCGA/TCGA-BRCA/bamfile_list.txt

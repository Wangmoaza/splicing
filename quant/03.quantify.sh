#!/bin/bash


IN_DIR="/home/haeun/2TB_disk/splicing/data/TCGA/TCGA-BRCA/bamfiles"
OUT_DIR="/home/haeun/2TB_disk/splicing/Analysis/Quant/ballgown"
REF="/home/omics/DATA1/06_haeun/splicing/Analysis/Quant/stringtie.BRCA_751.merged.gtf"

for i in $(ls $IN_DIR | grep "\.bam" | sed 's/.bam//')
do
/home/haeun/2TB_disk/tools/stringtie-1.3.5.Linux_x86_64/stringtie -p 4 -f 0.05 -a 5 -e -B \
    -G /home/haeun/2TB_disk/splicing/Analysis/Quant/merged/stringtie.BRCA_751.merged.tier_123.gtf \
    -A $OUT_DIR/${i}/${i}_tier123_gene.tab \
    -o $OUT_DIR/${i}/${i}_tier123_quant.gtf \
    $IN_DIR/${i}.bam
done

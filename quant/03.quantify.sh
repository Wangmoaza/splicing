#!/bin/bash


IN_DIR="/home/haeun/DATA/splicing/data/TCGA/TCGA-BRCA/bamfiles"
OUT_DIR="/home/haeun/DATA/splicing/Analysis/Quant/ballgown-all"
#REF="/home/omics/DATA1/06_haeun/splicing/Analysis/Quant/stringtie.BRCA_751.merged.gtf"

#for i in $(ls $IN_DIR | grep "\.bam" | sed 's/.bam//')
#while read line
#do
#	input=$(echo $line | sed 's/.bam//')
#	/home/haeun/DATA/tools/stringtie-1.3.4d.Linux_x86_64/stringtie -p 4 -e -B \
#    	-G /home/haeun/DATA/splicing/Analysis/Quant/merged/stringtie.BRCA_751.merged.gtf \
#    	-A $OUT_DIR/${input}/${input}_gene.tab \
#    	-o $OUT_DIR/${input}/${input}_quant.gtf \
#    	$IN_DIR/${input}.bam
#done < to_do_list.txt

ls $IN_DIR | egrep "bam$" | sed s/.bam// | \
	parallel --eta -k -j4 --load 80% --noswap '/home/haeun/DATA/tools/stringtie-1.3.4d.Linux_x86_64/stringtie -p 1 -e -B -G /home/haeun/DATA/splicing/Analysis/Quant/TCGA-BRCA/merged/stringtie.TCGA-BRCA.merged_strict.gtf -A /home/haeun/DATA/splicing/Analysis/Quant/TCGA-BRCA/ballgown-strict/{}/{}_gene.tab -o /home/haeun/DATA/splicing/Analysis/Quant/TCGA-BRCA/ballgown-strict/{}/{}_quant.gtf /home/haeun/DATA/splicing/data/TCGA/TCGA-BRCA/bamfiles/{}.bam'


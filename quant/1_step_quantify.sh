#!/bin/bash


cat /home/haeun/DATA/splicing/data/TCGA/TCGA-BRCA/All.HRD_sig3_median_groups_bams.txt | sed s/.bam// | \
	        parallel --eta -k -j4 --load 80% --noswap '/home/haeun/DATA/tools/stringtie-1.3.5.Linux_x86_64/stringtie -p 1 -e -B -G /home/haeun/DATA/splicing/data/GENCODE/gencode.v29.annotation.protein_coding.gtf -A /home/haeun/DATA1/splicing/1step_quant/{}/{}_gene.tab -o /home/haeun/DATA1/splicing/1step_quant/{}/{}_quant.gtf /home/haeun/DATA/splicing/data/TCGA/TCGA-BRCA/bamfiles/{}.bam'

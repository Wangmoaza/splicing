#!/bin/bash

bedtools window -a FDR_005_sig_pos.bed -b ../../data/GENCODE/gencode.v29lift37.annotation.protein_coding.gene.bed -w 1000000 > 1mb_window_genes.bed 

cut -f16 1mb_window_genes.bed | awk -F "; " '{ if ($3 ~ "gene_name") {print $3} else {print $4}}' | awk -F " " '{print $2}' | sed 's/\"//g' | sort -u > 1mb_window_genes.list.txt

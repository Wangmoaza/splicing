#!/usr/bin/env bash

cut -f2 /home/haeun/DATA/splicing/data/KEGG/KEGG_plus_curated_genes_list.ensembl.txt > tmp.txt

find /home/haeun/DATA/splicing/Analysis/Quant-all/htseq -type f | \
    parallel -j4 'grep -f tmp.txt {} > {//}_tier123/{/.}_tier123.txt'

rm tmp.txt

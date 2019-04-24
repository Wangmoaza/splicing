#!/bin/bash
cut -f2 $TIER > tmp.txt

find /home/haeun/DATA/splicing/Analysis/Quant-all/htseq -type f | \
    parallel -j4 'grep -f tmp.txt {} > {//}_tier123/{/.}_tier123.txt'

rm tmp.txt

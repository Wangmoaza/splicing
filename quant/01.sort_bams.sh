#!/bin/bash

DIR="/home/haeun/2TB_disk/splicing/data/TCGA/TCGA-BRCA/bamfiles"

find $DIR -type f | grep "bam" | parallel -j 2 'samtools sort {} -o {}.sorted.bam'

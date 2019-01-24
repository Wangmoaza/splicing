#!/usr/bin/env python2
import sys

INFILE = sys.argv[1]
OUTFILE = sys.argv[2]

with open(OUTFILE, 'w') as fw:
    fw.write("probeset_id\tCall\n")
    with open(INFILE, 'r') as fr:
        fr.readline() # skip header
        for line in fr.readlines():
            tokens = line.strip().split(',')
            snp, alt, genotype = tokens[4], tokens[7], tokens[13]
            fw.write("{0}\t{1}\n".format(snp, genotype.count(alt)))

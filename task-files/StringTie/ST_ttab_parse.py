#!/usr/bin/env python
# Script to parse Stringtie transcript coverage table into two column count matrix
# First Write April 27th, 2016
# William Poehlman <wpoehlm@clemson.edu>

import sys

#Load input file
FileName = sys.argv[1]
RawData = [line for line in open(FileName)]

#Print header
print 'TranscriptID' + '\t' + 'FPKM'

#Extract transcript ID's and FPKM values  
for line in RawData:
    line = line.strip()

    if not line or line.startswith('t_id'):
        continue
        
    try:
        t_id, chr, strand, start, end, t_name,\
        num_exons, length, gene_id, gene_name, cov, FPKM = line.split('\t')
        
    except ValueError:
        continue
    
    print t_name + '\t' + FPKM

#!/usr/bin/env python
# Script to parse Stringtie Or Cufflinks GTF Output files into two column matrix
# First Write April 27th, 2016
# William Poehlman <wpoehlm@clemson.edu>
import sys

#Load input file

FileName = sys.argv[1]
RawData = [line for line in open(FileName)]

#Extract transcript ID and FPKM values without redundancy

print 'TranscriptID' + '\t' + 'FPKM'
for line in RawData:
    line = line.strip()
    if not line or line.startswith('#'):
        continue 
        
    try:
        chromosome, software, feature, start,\
        stop, val, strand, point, expression = line.split('\t')
    
    except ValueError:
        continue
    
    if software =='StringTie' and feature == 'transcript':
        trans = expression.split(';')[2]
        Num = expression.split(';')[5]
        #print transcriptID + '\t' + FPKM 
        
        Tran = trans.split(' ')[2]
        TranscriptID = Tran.strip('"')
        
        
        Nums = Num.split(' ')[2]
        FPKM = Nums.strip('"')
        
         
        print TranscriptID + '\t' + FPKM
        
    if software == 'Cufflinks' and feature == 'transcript':
        trans = expression.split(';')[1]
        Num = expression.split(';')[2]
        
        Tran = trans.split(' ')[2]
        TranscriptID = Tran.strip('"')
        
        
        Nums = Num.split(' ')[2]
        FPKM = Nums.strip('"')
        
        
        print TranscriptID + '\t' + FPKM




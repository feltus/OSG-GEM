## William Poehlman ##
## LogParse.py v1.1 ##
## wpoehlm@clemson.edu ##
## First Write February 1st, 2017 ##
## This script merges the standard output log files of the OSG-GEM workflow for each sample, and prints a quality control report ##

import glob
import os
import sys

def ParseLogs(trimlist, maplist, software, sample):
    """
    Merges statistics from lists of split trimmomatic and hisat/tophat log files for a given sample
    
    
    arguments:
    trimlist --> a list of trimmomatic log files
    maplist --> a list of hisat2 or tophat2 log files
    software --> read mapping software that generated log files.  'hisat' or 'tophat'
    sample --> A sample number that is present at the beginning of each filename
    
    output:
    A dictionary, where the sample ID number is the key, and the following statistics are stored in the list value:
    Number of raw input read pairs
    Number of read pairs following trimming
    Number of concordant, single alignments
    Overall Alignment rate
    
    The overall Alignment rate for hisat2 is calculated as follows:
    
    (concordant single alignments + concordant multiple alignments + non-paired single alignments)/(cleaned input pairs)
    
    The overall Alignment rate for tophat2 is calculated as follows:
    
    ((left mapped reads + right mapped reads)/2)/(cleaned input pairs)
    
    Logfile names must be formatted as follows:
    
    <sample number>-<unique file identifier>-align_summary.txt (only for tophat)
    <sample number>-<unique file identifier>-fastq-out.txt (for hisat)
    <sample number>-unique file identifier>-fastq-trimmomatic.txt
    
    for example:
    2-000000-000000.fastq-align_summary.txt
    1-000000-000000.fastq-out.txt
    1-000000-000000.fastq-trimmomatic.txt    
    
    NOTE: do not mix log files for tophat and hisat
    
    """
    triminput = []
    trimoutput = []
    mapped = []
    concordant = []
    sample = str(sample)
    software = str(software)
    report = {}
    if software == 'hisat2':
        for file in trimlist:
            basename = os.path.splitext(file)[0]
            dataset = basename.split('-')[0]
            if dataset == sample:
                for line in open(file):
                    if 'Input Read Pairs' in line:
                        inp = int(line.split(' ')[3])
                        triminput.append(inp)
                        out = int(line.split(' ')[6])
                        trimoutput.append(out)
        for file in maplist:
            basename = os.path.splitext(file)[0]
            dataset = basename.split('-')[0]
            if dataset == sample:
                for line in open(file):
                    if 'aligned concordantly exactly 1 time' in line:
                        conc = int(line.split('(')[0])
                        concordant.append(conc)
                        mapped.append(conc)
                    elif 'aligned concordantly >1 times' in line:
                        mult = int(line.split('(')[0])
                        mapped.append(mult)
                    elif 'aligned exactly 1 time' in line:
                        singlemap = int(line.split('(')[0])
                        mapped.append(singlemap)
    
    elif software == 'tophat2':
        for file in trimlist:
            basename = os.path.splitext(file)[0]
            dataset = basename.split('-')[0]
            if dataset == sample:
                for line in open(file):
                    if 'Input Read Pairs' in line:
                        inp = int(line.split(' ')[3])
                        triminput.append(inp)
                        out = int(line.split(' ')[6])
                        trimoutput.append(out)
        for file in maplist:
            basename = os.path.splitext(file)[0]
            dataset = basename.split('-')[0]
            if dataset == sample:
                totalpair = []
                multipair = []
                discordpair = []
                leftmap = []
                rightmap = []
            
                count = 0
                for line in open(file):
                    count +=1 
                    if count == 11:
                        totalpairs = int(line.split(':')[1].strip())
                        totalpair.append(totalpairs)
                    elif count == 12:
                        multipairs = int(line.split(':')[1].split('(')[0].strip())
                        multipair.append(multipairs)
                    elif count == 13:
                        discordpairs = int(line.split('(')[0].strip())
                        discordpair.append(discordpairs)
                    elif count == 3:
                        leftmapped = int(line.split(':')[1].split('(')[0].strip())
                        leftmap.append(leftmapped)
                    elif count == 7:
                        rightmapped = int(line.split(':')[1].split('(')[0].strip())
                        rightmap.append(rightmapped)
                
                totalpair = sum(totalpair)
                multipair = sum(multipair)
                discordpair = sum(discordpair)
                leftmap = sum(leftmap)
                rightmap = sum(rightmap)
                mappedsingles = leftmap + rightmap
                concordantpairs = totalpair - multipair - discordpair   
                concordant.append(concordantpairs)
                mapped.append(mappedsingles)
    
    if software == 'tophat2':
        
        inputraw = sum(triminput)
        inputclean = sum(trimoutput)
        concordant = sum(concordant)
        mapped = sum(mapped)
        float1 = float(mapped)
        float2 = float(inputclean)
        float3 = float1/2
        rate = float3/float2
    
        report[sample] = [inputraw, inputclean, concordant, rate*100]
        return report
    
    elif software == 'hisat2':
        inputraw = sum(triminput)
        inputclean = sum(trimoutput)
        concordant = sum(concordant)
        mapped = sum(mapped)
        float1 = float(mapped)
        float2 = float(inputclean)
        rate = float1/float2
    
        report[sample] = [inputraw, inputclean, concordant, rate*100]
        return report
                         
                    

def main():

    software = sys.argv[1]

    #Print header of QC report
    print 'Sample' + '\t' + 'raw_input_pairs' + '\t' + 'cleaned_input_pairs' + '\t' + 'concordant_single_alignments' + '\t' + 'Overall_Alignment_Rate'

    #locate input log files and append to lists

    #tophat alignment summary
    maploglist = [file for file in glob.glob(('*summary.txt'))]

    if len(maploglist) == 0:
        #If there are no tophat alignment summaries, look for hisat standard output instead
        maploglist = [file for file in glob.glob('*fastq-out.txt')]

    trimloglist = [file for file in glob.glob('*-trimmomatic.txt')]
    if len(maploglist) != len(trimloglist):
        print "WARNING: the number of trimmomatic log files is different than the number of hisat log files. Output may have frameshift."
    inputlist = []
    for file in maploglist:
        basename = os.path.splitext(file)[0]
        dataset = basename.split('-')[0]
        if dataset not in inputlist:
            inputlist.append(dataset)
        
    #for each unique sample, call the ParseLogs function and print results in tab-delimited format        
    for dataset in inputlist:
        stats = ParseLogs(trimloglist, maploglist, software, dataset)
        print dataset + '\t' + str(stats[dataset][0]) + '\t' + str(stats[dataset][1]) + '\t' + str(stats[dataset][2]) + '\t' + str(stats[dataset][3]) + '%'


if __name__ == '__main__':
    main()








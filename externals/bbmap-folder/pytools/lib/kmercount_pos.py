#!/usr/bin/env python

from __future__ import division
import sys
import os
import argparse
import string
from collections import Counter
import random
sys.path.append(os.path.dirname(__file__));
import readSeq
import multiprocessing

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

KMER = 16
SUBLEN = 1000000


def getArgs():
    parser = argparse.ArgumentParser(description="Count occurance of database kmers in reads", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-k', default=KMER, dest='kmer', metavar='<int>', type=int, help="kmer length")
    parser.add_argument('-l', dest='sublen', metavar='<int>', type=int, help="perform analysis on first <int> bases [RDLEN - K + 1]")
    parser.add_argument('-c', default=2, dest='cutoff', metavar='<int>', type=int, help="minimum allowed coverage")
    parser.add_argument('-t', default=30, dest='targetCov', metavar='<int>', type=int, help="target coverage")
    parser.add_argument('-p', '--plot', dest='plot', type=str, default=None, metavar="<file>", help='plot data and save as png to <file>')
    parser.add_argument('fastaFile', type=str, help='Input FASTA file(s). Text or gzip')
    parser.add_argument('fastqFile', nargs='+', help='Input FASTQ file(s). Text or gzip')
    args = parser.parse_args()
    return args

def getMers(seq, merLen):
    for i in xrange(min(len(seq) - merLen + 1,SUBLEN)):
        yield seq[i:i+merLen]

complements = string.maketrans('acgtACGT', 'tgcaTGCA')
def revComp(seq):
    revCompSeq = seq.translate(complements)[::-1]
    return revCompSeq

def tallyPositions(seq,counts,sublen=None):
    readPos = 0
    for mer in getMers(seq, KMER):
        if mer in merCnts:
            if readPos > len(counts):
                counts.append(1)
            else:
                counts[readPos] += 1
        readPos += 1
        if sublen:
            if readPos == sublen:
                break
#end tallyPositions
            

def main():

    """
    kmer based normization
    call rand once
    004 
        if avecov1|2 < mincov: 
            trash 
        elif avecov1|2 < target: 
            keep
        elif random < target/avecov1|2: 
            keep
        else:
            trash
    """

    args = getArgs()

    global KMER
    KMER = args.kmer
    
    out=sys.stdout
    
    #check to make sure input files exist
    for fq in args.fastqFile:
        if not os.path.exists(fq):
            sys.stderr.write("ERROR: Input file '%s' does not exist. Exiting.\n" % fq)
            sys.exit()
     
    #count mers / create database
    sys.stderr.write("Making k-mer database\n")
    global merCnts
    merCnts = Counter()
    for record in readSeq.readSeq(args.fastaFile,fileType='fasta'):
        for mer in getMers(record.seq, args.kmer):
            sortMers = sorted([mer,  revComp(mer)])
            merCnts[mer] += 1
            merCnts[revComp(mer)] += 1
    #normalize reads
    
    sys.stderr.write("Tallying occurrences of database kmers in reads\n")

    seqIt = readSeq.readSeq(fq,paired=True)
    record1,record2 = seqIt.next()
    readLen = len(record1.seq)
    sys.stderr.write("Read length = %d\n" % readLen)
    
    tallyLen = readLen - args.kmer + 1
    if args.sublen:
        if args.sublen > tallyLen:
            sys.stderr.write("sublen (-l) must be less than readlen - k + 1 : (found reads of length %d\n" % readLen)
            sys.exit(-1)
        tallyLen = args.sublen

    counts1 = [0 for x in range(tallyLen)]
    counts2 = [0 for x in range(tallyLen)]

    tallyPositions(record1.seq,counts1,tallyLen)
    tallyPositions(record2.seq,counts2,tallyLen)
    total_reads = 1
    
    fqName = ""
    
    for fq in args.fastqFile:
        for record1, record2 in seqIt:
            tallyPositions(record1.seq,counts1,tallyLen)
            tallyPositions(record2.seq,counts2,tallyLen)
            total_reads += 1
        fqName += fq+" "
    

    counts1_perc = [ 100 * float(x)/total_reads for x in counts1 ]
    counts2_perc = [ 100 * float(x)/total_reads for x in counts2 ]
    out.write("#pos\tread1_count\tread1_perc\tread2_count\tread2_perc\n")
    for i in range(tallyLen):
        out.write("%i\t%i\t%0.2f\t%i\t%0.2f\n" % (i+1,counts1[i],counts1_perc[i],counts2[i],counts2_perc[i]))

    if args.plot:
        sys.stderr.write("Plotting data. Saving to %s\n" % args.plot)
        plt.ioff()
        xcoord = range(1,tallyLen+1)
        plt.plot(xcoord,counts1_perc,color="red",linewidth=1.0,linestyle="-",label="Read 1")
        plt.plot(xcoord,counts2_perc,color="green",linewidth=1.0,linestyle="-",label="Read 2")
        leg_loc="upper right"
        max_cnt = 0
        max_pos = 0
        for i in range(len(counts1)):
            if counts1[i] > max_cnt or counts2[i] > max_cnt:
                max_cnt = counts1[i] if counts1[i] > counts2[i] else counts2[i]
                max_pos = i
        if max_pos > 0.5*len(counts1) :
            leg_loc="upper left"
        plt.legend(loc=leg_loc,prop={'size':10})
        plt.xlabel("Read Position")
        plt.ylabel("Percent Reads with Database k-mer")
        plt.title("Occurrence of reference k-mers (k = %i) in \n%s (# reads = %d)" % (args.kmer,fqName,total_reads))
        plt.savefig(args.plot)

 
if __name__ == '__main__':

    try:
        main()
    except KeyboardInterrupt:
        pass
        

#!/usr/bin/env python

import sys
import math
import random
import fileinput
import re
#import os
import copy

class readSeq:

    def __init__(self, fileName, paired = False, fileType = 'fastq', compression = 'auto', sampleRate = None):
        self.fileName = fileName
        self.fileType = fileType
        self.fileType = self.__detectFileFormat()
        self.compression = compression
        self.__readObj = self.readFile()
        self.paired = paired
        self.sampleRate = sampleRate
        
    def __iter__(self):
        return self

    def next(self):
        record = self.__readObj.next()
        self.header = record['header']
        self.seq = record['seq']
        self.seqNum = record['seqNum']
        if self.fileType == 'fastq':
            self.qual = record['qual']

        if self.paired is True:
            self2 = copy.copy(self)
            record = self.__readObj.next()
            self2.header = record['header']
            self2.seq = record['seq']
            self2.seqNum = record['seqNum']
            if self.fileType == 'fastq':
                self2.qual = record['qual']
            return self, self2
        else:
            return self

    def __detectFileFormat(self):
        if self.fileType != 'auto':
            return self.fileType
         
        if ( re.search('fasta$', self.fileName) or 
             re.search('fa$', self.fileName)    or
             re.search('fna$', self.fileName) ):
            self.fileType = 'fasta'

        elif ( re.search('fastq$', self.fileName) or
               re.search('fq$', self.fileName) ) :
            self.fileType = 'fastq'

        return self.fileType

    def readFile(self):

        record = dict()       
 
        # allow user to specify if compression should be auto detected using file names
        # specified by fileinput documentation
        if self.compression == 'auto':
            inFH = fileinput.FileInput(self.fileName, openhook=fileinput.hook_compressed)    
        else:
            inFH = fileinput.FileInput(self.fileName)
  
        #TODO RAISE exception if file doesn't exist
  
        # read fastq
        if self.fileType == 'fastq':
            self.seqNum = 0
            for record['header'] in inFH:

                record['seqNum'] = int(math.ceil(inFH.lineno()/8.0))
                record['seq'] = inFH.readline().strip()
                record['r1Plus'] = inFH.readline()
                record['qual'] = inFH.readline().strip()
                     
                if not re.search('^@', record['header'], flags=0):
                    pass
                    #TODO add exception for not being a fastq
                          
                record['header'] = re.sub('^@', '', record['header'].strip())
                yield record

        elif self.fileType == 'fasta':
            record['header'] = None
            record['seq'] = None
            record['seqNum'] = 0
            
            for line in inFH:
 
                line = line.strip()
                if re.search('^>', line):
                    line = re.sub('^>', '', line)
                    if record['header'] is None:
                        record['header'] = line
                        if not re.search('>', record['header']):
                            pass
                            # TODO add exception for not being fasta   
 
                    if record['seq'] is not None:
                        record['seqNum'] = record['seqNum'] + 1
                        
                        if self.sampleRate:
                        #subsampling of reads is desired
                            if random.random() < self.sampleRate:
                                yield record
                        else:
                            yield record
                        record['seq'] = None

                    record['header'] = line                
                else:
                    if record['seq'] is not None:
                        record['seq'] += line
                    else:
                        record['seq'] = line
                        
            record['seqNum'] = record['seqNum'] + 1
            
            if self.sampleRate:
                #subsampling of reads is desired
                if random.random() < self.sampleRate:
                    yield record
            else:
                yield record

        try:
            inFH = fileinput.close()
        except:
            pass

    #except:
        sys.exc_info()[0]

    def printRecord(self):
        if self.fileType == 'fastq':
            print("@%s" % self.header)
            print(self.seq)
            print("+") 
            print(self.qual)
        elif self.fileType == 'fasta':
            print(">%s" % self.header)
            print(self.seq)
    
if __name__ == '__main__':

    for inputfile in sys.argv[1:]:
        #get a pair of reads
        for r1, r2 in readSeq(inputfile, paired=True):
            print r1.header
            print r2.header

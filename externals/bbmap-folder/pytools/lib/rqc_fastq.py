#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import re
import gzip
import argparse

pipe_root = os.path.abspath(os.path.dirname(os.path.abspath(__file__)) + '/' + os.pardir)
sys.path.append(pipe_root + '/lib')   # common

from common import run_command
# from db_access import jgi_connect_db


"""
These are functions for fastq file analysis


"""

def _err_return(msg, log=None):
    if log:
        log.error(msg)
    print msg
    return None


'''

For a given fastq file (zipped or unzipped), check the format by comparing the string lengths of seq and quality

'''
def check_fastq_format(fastq, log=None):
    if not fastq:
        return _err_return('Function read_length_from_file() requires a fastq file.', log)

    if not os.path.isfile(fastq):
        return _err_return('The given fastq file [%s] does not exist' % fastq, log)

    if fastq.endswith(".gz"):
        fh = gzip.open(fastq, 'r')
    else:
        fh = open(fastq, 'r')

    lCnt = 0    # line counter
    seq = None
    lstr = None
    missing = False

    while 1:
        lstr = fh.readline()
        if not lstr:
            break

        lCnt += 1
        seqLen = 0
        qualLen = 0

        idx = lCnt % 4  # line index in the 4-line group

        if idx == 1:
            continue
        elif idx == 2:
            seq = lstr.strip()
        elif idx == 3:
            continue
        elif idx == 0:
            seqLen = len(seq)
            qualLen = len(lstr.strip())

            if seqLen != qualLen:
                missing = True
                break

    fh.close()

    if missing:
        log.error("Incorrect fastq file: missing quality score character or base character.\n seq=%s\n qual=%s" % (seq, lstr))
        return -1

    if lCnt % 4 != 0:
        log.error("Incorrect fastq file: missing fastq record item. Number of lines in the output fastq file = %d" % (lCnt))
        return -2


    return lCnt



'''
    Ref : JGI_Utility::read_length_from_file
    The first read count (read 1 and read 2) meets the sszie will stop the record reading.
'''
def read_length_from_file(fastq, log=None, ssize=10000):
    ''' For a given fastq file (zipped or unzipped), reads the first *ssize* reads to compute average length
        and return (avgLen1, avgLen2, isPE)
        The file scanning will stop when the ssize is met by either read (1 or 2).
    '''

    is_pe = False   # def to none paired end
    readi = 0   # counter for pe read 1
    readj = 0   # counter for pe read 2
    leni = 0    # total length of read 1
    lenj = 0    # total length of read 2

    if not fastq:
        return _err_return('Function read_length_from_file() requires a fastq file.', log)

    if not os.path.isfile(fastq):
        return _err_return('The given fastq file [%s] does not exist' % fastq, log)

    if fastq.endswith(".gz"):
        fh = gzip.open(fastq, 'r')
    else:
        fh = open(fastq, 'r')

    lCnt = 0    # line counter
    done = False

    while not done:
        lstr = fh.readline()
        lstr = lstr.rstrip()

        lCnt += 1
        idx = lCnt % 4  # line index in the 4-line group
        if idx == 1:
            header = lstr
        elif idx == 2:
            seq = lstr
        elif idx == 3:
            plus = lstr
        elif idx == 4:
            quality = lstr
        else:
            if not header or not seq:   # end of file
                done = True

            if header.find('enter') != -1 and header.find('exit') != -1:    # copied from perl's logic
                continue

            match = re.match(r'^@([^/]+)([ |/]1)', header)  # for pe read 1
            aLen = len(seq.strip())
            if match:
                # to let read2 match up
                readi += 1
                if readi > ssize:   #read 1 meet the max count; done regardless situations of read2 count.
                    readi -= 1
                    done = True
                else:
                    leni += aLen
            else:
                match = re.match(r'^@([^/]+)([ |/]2)', header)  # for pe read 2
                if match:
                    readj += 1
                    if readj > ssize:   #read2 meet the max count; done regardless situation of read1 count
                        readj -= 1
                        done = True
                    else:
                        lenj += aLen
                    if leni > 0:  ### Only set is_pe to true if leni > 0, which means the 1st read in the pair was found
                        is_pe = True
    fh.close()

    # debug to be sure the max are met properly for both reads
    #print('read1len=%d; read1cnt=%d; read2len=%d; read2cnt=%d' %(leni, readi, lenj, readj))
    if leni > 0 and readi > 0:
        leni = leni / readi
    if lenj > 0 and readj > 0:
        lenj = lenj / readj

    return (leni, lenj, is_pe)


def get_working_read_length(fastq, log):
    read_length = 0
    read_length_1 = 0
    read_length_2 = 0

    (read_length_1, read_length_2, is_pe) = read_length_from_file(fastq, log)

    if not is_pe:
        log.info("It is NOT pair-ended. read_length_1 %s" % (read_length_1))
        read_length = read_length_1
    else:
        log.info("It is pair-ended. read_length_1 %s read_length_2 %s" % (read_length_1, read_length_2))
        if read_length_1 != read_length_2:
            log.warning("The lengths of read 1 (" + str(read_length_1) + ") and read 2 (" + str(read_length_2) + ") are not equal")

        if read_length_1 < read_length_2:
            read_length = read_length_1
        else:
            read_length = read_length_2

        #if read_length < 10:
        #    log.error("File name: " + fastq + ". Read length is less than 10 bps. Is this paired end? " + str(is_pe) + ". Read one length: " + str(read_length_1) + "; Read two length: " + str(read_length_2) )
        #    return (0, read_length_1, read_length_2, is_pe)

    return (read_length, read_length_1, read_length_2, is_pe)



def read_count(fpath):
    'return the raw read count in a fastq file, assuming each record occupy 4 lines in file.'
    if os.path.isfile(fpath):
        cdir = os.path.dirname(__file__)
        cmd = os.path.join(cdir, '../../testformat2.sh')
        cmd = '%s in=%s | grep "^Reads"' % (cmd, fpath)
        stdOut, stdErr, exitCode = run_command(cmd, True)

        if exitCode == 0:
            toks = stdOut.strip().split()
            if len(toks) == 2:
                return int(toks[1])
            else:
                print('error in %s:read_count (%s): wrong run output [%s]' % (__file__, cmd, stdOut.strip()))
        else:
            print('error in %s:read_count (%s): %s' % (__file__, cmd, stdErr.strip()))

    return None

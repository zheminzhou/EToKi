#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
stand alone filter script based on jgi-rqc-pipeline/filter/fastq_filter2.py v3.2.8

Command
    $ filter.py -f FASTQ.GZ -o OUT_DIR --skip-blast -html

Outputs
  - normal filter outputs + index.html

Created: March 15, 2018

Shijie Yao (syao@lbl.gov)

Revision:

"""

import os
import sys
import argparse
import datetime
import shutil
import glob
import re
# import numpy as np
# import matplotlib
# matplotlib.use("Agg") ## This needs to skip the DISPLAY env var checking
# import matplotlib.pyplot as plt
# import mpld3
from pprint import pprint

## append the pipeline lib and tools relative path:
SRC_ROOT = os.path.abspath(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(SRC_ROOT + "/lib")   # common

import rqc_fastq as fastqUtil
from common import get_logger, get_status, checkpoint_step, set_colors, get_subsample_rate, append_rqc_file, append_rqc_stats
from os_utility import run_sh_command, make_dir_p

from rqc_utility import get_dict_obj, pipeline_val
from readqc import do_html_body as do_readqc_html_body
from html_utility import html_tag


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## global vars
VERSION = "1.0.0"
SCRIPT_NAME = __file__
logLevel = "DEBUG"
DEBUG = False
color = {}
color = set_colors(color, True)

READQC_DIR = 'filtered-readqc'

PYDIR = os.path.abspath(os.path.dirname(__file__))
BBDIR = os.path.join(PYDIR, os.path.pardir)

FILTER_READ_COUNT = "read count"
FILTER_READ_SAMPLED_COUNT = "read sampled count"

STATUS_LOG_FNAME = 'status.log'
BB_STATS_LIST_FILE_NAME = 'filterStats.txt'
STATS_LIST_FILE_NAME = 'filter-stats.txt'

PIPE_START        = 'start'
RQCFILTER_START   = 'start rqcfilter'
RQCFILTER_END     = 'end rqcfilter'
POST_START        = 'start post process'
POST_END          = 'end post process'
QC_START        = 'start qc'
QC_END          = 'end qc'
HTML_START        = 'start html generation'
HTML_END          = 'end html generation'
PIPE_COMPLETE     = 'complete'

STEP_ORDER = {
    PIPE_START        : 0,

    RQCFILTER_START   : 10,
    # rqcfilter.sh recorded steps:
        'clumpify start'            : 11,
        'clumpify finish'           : 12,
        'ktrim start'               : 20,
        'ktrim finish'              : 21,
        'delete temp files start'   : 30,
        'delete temp files finish'  : 31,
        'filter start'              : 40,
        'filter finish'             : 41,
        'delete temp files start'   : 50,
        'delete temp files finish'  : 51,
        'short filter start'        : 60,
        'short filter finish'       : 61,
        'delete temp files start'   : 70,
        'delete temp files finish'  : 71,
        'removeCommonMicrobes start'    : 80,
        'removeCommonMicrobes finish'   : 81,
        'delete temp files start'       : 90,
        'delete temp files finish'      : 91,
        'dehumanize start'              : 100,
        'dehumanize finish'             : 101,
        'delete temp files start'       : 110,
        'delete temp files finish'      : 111,
        'merge start'                   : 120,
        'merge finish'                  : 121,
        'khist start'                   : 130,
        'khist finish'                  : 131,
        'rqcfilter complete'            : 200,


    RQCFILTER_END     : 300,

    POST_START        : 381,
    POST_END          : 382,

    QC_START          : 400,
    QC_END            : 410,

    HTML_START        : 420,
    HTML_END          : 430,

    PIPE_COMPLETE     : 999
}

FILTER_METHODS_TXT = {"DNA": "dna.txt",
                      "FUNGAL": "fungal.txt",
                      "METAGENOME": "metagenome.txt",
                      "VIRAL-METAGENOME": "viral-metagenome.txt",
                      "ISO": "iso.txt",
                      "SAG": "sag.txt",
                      "CELL-ENRICHMENT": "cell-enrichment.txt",
                      "PLANT-2X150": "plant-2x150.txt",
                      "PLANT-2X250": "plant-2x250.txt",
                      "RNA": "rna.txt",
                      "RNAWOHUMAN": "rnawohuman.txt",
                      "3PRIMERNA": "3primerna.txt",
                      "SMRNA": "smrna.txt",
                      "METATRANSCRIPTOME": "mtaa.txt",
                      "MTF": "mtaa.txt",
                      "LFPE": "lfpe.txt",
                      "CLRS": "clrs.txt",
                      "CLIP-PE": "clip-pe.txt",
                      "NEXTERA-LMP": "nextera-lmp.txt",
                      "NEXTERA": "nextera-lmp.txt",
                      "NEXTSEQ": "nextseq.txt",
                      "ITAG": "itag.txt",
                      "MICROTRANS": "microtrans.txt",
                      "BISULPHITE": "bisulphite.txt",
                      "CHIPSEQ": "chip-seq.txt"}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## pipeline function definitions
def log_and_print(msg):
    print(msg)
    log.info(msg)


"""
Run rqcfilter.sh

@param fastq: the raw fastq input file [input of filtering]
@param outDir: the output directory
@param prodType: product type
@param status: current status
@return outFastqFile, outRrnaFastqFile

"""
def run_rqcfilter(infastq, outDir, prodType, status, enableRmoveMicrobes, enableAggressive, disableRmoveMicrobes, disableClumpify, taxList, rdb, log):
    log_and_print("\n\n%s - RUN RQCFILTER <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<%s\n" % (color['pink'], color['']))

    make_dir_p(outDir)
    opt = None

    extraOptions = ""
    optionFile = ""

    if prodType.endswith("-OLD"):
        extraOptions = " barcodefilter=f chastityfilter=f"
        prodType = prodType.replace("-OLD", "")

    if prodType in FILTER_METHODS_TXT:
        optionFile = os.path.join(SRC_ROOT, "filter_param/", FILTER_METHODS_TXT[prodType].replace(".txt", ".config"))
        log_and_print("Filter option file: %s" % optionFile)
        opt = open(optionFile, 'r').readline().rstrip()
    else:
        log_and_print("The product type, %s, is not supported yet." % prodType)
        sys.exit(2)

    assert (opt), "Null filter options."

    if prodType in ("METATRANSCRIPTOME", "MTF"):
        if infastq.endswith(".gz"):
            opt += " outribo=%s " % (os.path.basename(infastq).replace(".fastq", ".rRNA.fastq"))
        else:
            opt += " outribo=%s " % (os.path.basename(infastq).replace(".fastq", ".rRNA.fastq.gz"))

    if enableRmoveMicrobes:
        if opt.find("removemicrobes=f") != -1:
            opt = opt.replace("removemicrobes=f", "removemicrobes=t")
        opt += " removemicrobes=t "

    if disableRmoveMicrobes:
        if opt.find("removemicrobes=t") != -1:
            opt = opt.replace("removemicrobes=t", "removemicrobes=f")
        else: opt += " removemicrobes=f "

    if enableAggressive:
        opt += " aggressive=t microbebuild=3 "

    if taxList:
        opt += " taxlist=%s " % (taxList)

    ## Temp set clumpify=t for all prod types (RQC-890)
    if not disableClumpify:
        opt += " clumpify=t "
    else:
        opt += " clumpify=f "

    opt += " tmpdir=null "
    opt += extraOptions


    cmd = os.path.join(BBDIR, "rqcfilter.sh")
    filterLogFile = os.path.join(outDir, "filter.log")
    cmdStr = "%s in=%s path=%s %s usejni=f rqcfilterdata=%s > %s 2>&1" % (cmd, infastq, outDir, opt, rdb, filterLogFile)

    rtn = [None, status]
    outFastqFile = None

    shFileName = "%s/filter.sh" % outDir

    def find_filtered_fastq(adir):
        outFastqFile = None

        searchPatt = os.path.join(adir, "*.fastq.gz")
        outFastqFileList = glob.glob(searchPatt)

        assert len(outFastqFileList) >= 1, "ERROR: cannot find *.fastq.gz output file."
        for f in outFastqFileList:
            f = os.path.basename(f)
            t = f.split(".")

            if t[-3] not in ("frag", "singleton", "unknown", "rRNA", "lmp"):
                filterCode = t[-3]
            elif t[-3] == "lmp": ## nextera
                filterCode = t[-4]

            if len(t) == 7: ## ex) 12345.1.1234.ACCCC.anqdpht.fastq.gz
                fileNamePrefix = '.'.join(t[:4])
            elif len(t) == 6: ## ex) 6176.5.39297.anqrpht.fastq.gz
                fileNamePrefix = '.'.join(t[:3])
            else:
                log.warning("Unexpected filtered file name, %s", outFastqFileList)
                fileNamePrefix = '.'.join(t[:-3])
                log_and_print("Use %s as file prefix." %fileNamePrefix)

        assert filterCode and fileNamePrefix, "ERROR: unexpected filter file name: %s" % (outFastqFileList)

        of = os.path.join(adir, '.'.join([fileNamePrefix, filterCode, "fastq.gz"]))
        lof = os.path.join(adir, '.'.join([fileNamePrefix, filterCode, "lmp.fastq.gz"]))
        if os.path.isfile(of):
            outFastqFile = of
        elif os.path.isfile(lof):
            outFastqFile = lof
        else:
            log.error("Cannot find fastq.gz file.")

        # rename output file to *.filtered.fastq.gz
        f = os.path.basename(outFastqFile)
        t = f.split(".")
        if t[-3] != 'filtered':
            fto = os.path.join(adir, '.'.join(['.'.join(t[:-3]), 'filtered', "fastq.gz"]))
            shutil.move(outFastqFile, fto)
            outFastqFile = fto

        return outFastqFile

    def find_filter_number(outFastqFile):
        filteredReadNum = fastqUtil.check_fastq_format(outFastqFile) / 4

        if filteredReadNum < 0:
            log_and_print("RUN RQCFILTER - filtered fastq format error: %s." % outFastqFile)
        return filteredReadNum

    if STEP_ORDER[status] < STEP_ORDER[RQCFILTER_END]:

        create_shell(shFileName, (cmdStr,))

        log_and_print("rqcfilter cmd=[%s]" % cmdStr)
        log_and_print("sh file name=[%s]" % shFileName)

        stdOut, stdErr, exitCode = run_sh_command(shFileName, True, log, True) ## stdOut of 0 is success

        if exitCode != 0:
            log.error("Failed to run : %s, stdout : %s, stderr: %s", shFileName, stdOut, stdErr)
            return rtn

        outFastqFile = find_filtered_fastq(outDir)


        filteredReadNum = find_filter_number(outFastqFile)
        log_and_print("Read counts after RQCFILTER step = %d" % filteredReadNum)

        log_and_print("RUN RQCFILTER - completed")
        checkpoint(RQCFILTER_END, status)
        status = RQCFILTER_END

        if filteredReadNum == 0:
            log.warning("No reads left after filtering")
            checkpoint(PIPE_COMPLETE, status)
            with open(BB_STATS_LIST_FILE_NAME, 'a') as fh:
                write_stats(fh, FILTER_READ_COUNT, 0, log)
                write_stats(fh, FILTER_READ_BASE_COUNT, 0, log)
    else:
        log_and_print("No need to rerun RQCFILTER step, get filtered files and stats ... ")
        outFastqFile = find_filtered_fastq(outDir)

    rtn = [outFastqFile, status]
    return rtn

"""
gzip the final filtered fastq file, and generate STATS_LIST_FILE_NAME log files.

@param fastq: the raw fastq file
@param outDir: output dir
@param filteredFastq: the filtered fastq file
@param status: where the pipeline was at by last run
@param log

@return filteredFastq or None if error
"""
def post_process(fastq, outDir, filteredFastq, status, log):
    ## obtain read counts from input and filtered fastq files and save the values to STATS_LIST_FILE_NAME file;
    ## compress the filtered fastq file
    log_and_print("\n\n%s - RUN POST PROCESS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<%s\n" % (color['pink'], color['']))
    if STEP_ORDER[status] < STEP_ORDER[POST_END]:
        checkpoint(POST_START, status)
        rawCnt = 0
        rawBaseCnt = 0
        newCnt = 0
        newBaseCnt = 0

        stats = get_dict_obj(BB_STATS_LIST_FILE_NAME)
        rawCnt = pipeline_val('inputReads', {'type': 'int', 'vtype': 'numeric'}, stats)
        rawBaseCnt = pipeline_val('inputBases', {'type': 'int', 'vtype': 'numeric'}, stats)
        newCnt = pipeline_val('outputReads', {'type': 'int', 'vtype': 'numeric'}, stats)
        newBaseCnt = pipeline_val('outputBases', {'type': 'int', 'vtype': 'numeric'}, stats)

        readCounts = {}

        readRmPct = 100.0 * ((rawCnt - newCnt) / float(rawCnt))
        baseRmPct = 100.0 * ((rawBaseCnt - newBaseCnt) / float(rawBaseCnt))
        readCounts['readRmPct'] = '%.3f' % readRmPct
        readCounts['baseRmPct'] = '%.3f' % baseRmPct

        refStats = {}
        filterLogStat = {}

        cardinality = None
        bbdukVersion = None
        bbmapVersion = None

        if os.path.isfile("filter.log"):

            with open(os.path.join(outDir, "filter.log"), "r") as FLFH:
                isContamNumChecked = False ## Contamination will be done twice for removeribo or for MTF
                isKtrimmedTotalRemovedNumChecked = False ## for parsing "Total Removed" after ktrimming

                for l in FLFH:
                    if l.startswith("Input:"):
                        toks = re.findall("(\d+.\d*)", l.rstrip())
                        assert len(toks) == 2
                        if 'adaptertriminput' not in filterLogStat:
                            filterLogStat["adaptertriminput"] = {"numreads": toks[0], "numbases": toks[1]}
                        elif 'contamtriminput' not in filterLogStat:
                            filterLogStat["contamtriminput"] = {"numreads": toks[0], "numbases": toks[1]}

                    elif l.startswith("FTrimmed:"):
                        toks = re.findall("(\d+.\d*)", l.rstrip())
                        assert len(toks) == 4
                        filterLogStat["ftrimmed"] = {"numreads": toks[0], "percreads": toks[1], "numbases": toks[2], "percbases": toks[3]}
                    elif l.startswith("KTrimmed:"):
                        toks = re.findall("(\d+.\d*)", l.rstrip())
                        assert len(toks) == 4
                        filterLogStat["ktrimmed"] = {"numreads": toks[0], "percreads": toks[1], "numbases": toks[2], "percbases": toks[3]}
                        isKtrimmedTotalRemovedNumChecked = True

                    ## RQCSUPPORT-1987
                    elif l.startswith("Total Removed:") and isKtrimmedTotalRemovedNumChecked:
                        toks = re.findall("(\d+.\d*)", l.rstrip())
                        assert len(toks) == 4
                        filterLogStat["ktrimmed_total_removed"] = {"numreads": toks[0], "percreads": toks[1], "numbases": toks[2], "percbases": toks[3]}
                        isKtrimmedTotalRemovedNumChecked = False

                    elif l.startswith("Trimmed by overlap:"):
                        toks = re.findall("(\d+.\d*)", l.rstrip())
                        assert len(toks) == 4
                        filterLogStat["trimmedbyoverlap"] = {"numreads": toks[0], "percreads": toks[1], "numbases": toks[2], "percbases": toks[3]}

                    elif l.startswith("Result:"):
                        toks = re.findall("(\d+.\d*)", l.rstrip())
                        assert len(toks) == 4
                        if 'adaptertrimresult' not in filterLogStat:
                            filterLogStat["adaptertrimresult"] = {"numreads": toks[0], "percreads": toks[1], "numbases": toks[2], "percbases": toks[3]}
                        elif 'contamtrimresult' not in filterLogStat:
                            filterLogStat["contamtrimresult"] = {"numreads": toks[0], "percreads": toks[1], "numbases": toks[2], "percbases": toks[3]}

                    elif l.startswith("Unique 31-mers:"):
                        toks = re.findall("(\d+.\d*)", l.rstrip())
                        assert len(toks) == 2 or len(toks) == 1
                        if 'adaptertrimunique31mers' not in filterLogStat:
                            if len(toks) == 2:
                                filterLogStat["adaptertrimunique31mers"] = {"num": toks[1]}
                            else:
                                filterLogStat["adaptertrimunique31mers"] = {"num":"0"}
                        else:
                            if len(toks) == 2:
                                filterLogStat["contamtrimunique31mers"] = {"num": toks[1]}
                            else:
                                filterLogStat["contamtrimunique31mers"] = {"num":"0"}

                    elif not isContamNumChecked and l.startswith("Contaminants:"):
                        toks = re.findall("(\d+.\d*)", l.rstrip())
                        assert len(toks) == 4
                        filterLogStat["contaminants"] = {"numreads": toks[0], "percreads": toks[1], "numbases": toks[2], "percbases": toks[3]}
                        isContamNumChecked = True

                    elif l.startswith("QTrimmed:"):
                        toks = re.findall("(\d+.\d*)", l.rstrip())
                        assert len(toks) == 4
                        filterLogStat["qtrimmed"] = {"numreads": toks[0], "percreads": toks[1], "numbases": toks[2], "percbases": toks[3]}

                    elif l.startswith("Short Read Discards:"):
                        toks = re.findall("(\d+.\d*)", l.rstrip())
                        assert len(toks) == 4
                        filterLogStat["shortreaddiscards"] = {"numreads": toks[0], "percreads": toks[1], "numbases": toks[2], "percbases": toks[3]}
                    elif l.startswith("Low quality discards:"):
                        toks = re.findall("(\d+.\d*)", l.rstrip())
                        assert len(toks) == 4
                        filterLogStat["lowqualitydiscards"] = {"numreads": toks[0], "percreads": toks[1], "numbases": toks[2], "percbases": toks[3]}

                    elif l.startswith("BBDuk version"):
                        toks = re.findall("(\d+.\d*)", l.rstrip())
                        assert len(toks) == 1
                        bbdukVersion = toks[0]
                    elif l.startswith("BBMap version"):
                        toks = re.findall("(\d+.\d*)", l.rstrip())
                        assert len(toks) == 1
                        bbmapVersion = toks[0]

                    ## BBDuk 36.12 06272016
                    elif l.startswith("Adapter Sequence Removed:"):
                        toks = re.findall("(\d+.\d*)", l.rstrip())
                        assert len(toks) == 4
                        filterLogStat["adaptersequenceremoved"] = {"numreads": toks[0], "percreads": toks[1], "numbases": toks[2], "percbases": toks[3]}
                    elif l.startswith("Synthetic Contam Sequence Removed:"):
                        toks = re.findall("(\d+.\d*)", l.rstrip())
                        assert len(toks) == 4
                        filterLogStat["syntheticcontamsequenceremoved"] = {"numreads": toks[0], "percreads": toks[1], "numbases": toks[2], "percbases": toks[3]}

                    ## 08112016
                    elif l.startswith("Short Synthetic Contam Sequence Removed:"):
                        toks = re.findall("(\d+.\d*)", l.rstrip())
                        assert len(toks) == 4
                        filterLogStat["shortsyntheticcontamsequenceremoved"] = {"numreads": toks[0], "percreads": toks[1], "numbases": toks[2], "percbases": toks[3]}

                    elif l.startswith("Ribosomal Sequence Removed:"):
                        toks = re.findall("(\d+.\d*)", l.rstrip())
                        assert len(toks) == 4
                        filterLogStat["ribosomalsequenceremoved"] = {"numreads": toks[0], "percreads": toks[1], "numbases": toks[2], "percbases": toks[3]}

                    ## BBMap 36.12 06272016
                    elif l.startswith("Human Sequence Removed:"):
                        toks = re.findall("(\d+.\d*)", l.rstrip())
                        assert len(toks) == 4
                        filterLogStat["humansequenceremoved"] = {"numreads": toks[0], "percreads": toks[1], "numbases": toks[2], "percbases": toks[3]}

                    ## RQC-862, RQC-880
                    elif l.startswith("Microbial Sequence Removed:"):
                        toks = re.findall("(\d+.\d*)", l.rstrip())
                        assert len(toks) == 4
                        filterLogStat["microbialremoved"] = {"numreads": toks[0], "percreads": toks[1], "numbases": toks[2], "percbases": toks[3]}


        ##
        ## refStats.txt format
        ##
        ## name %unambiguousReads   unambiguousMB   %ambiguousReads ambiguousMB unambiguousReads    ambiguousReads
        ## human_masked 85.24693    498.92052   0.09378 0.55290 3350692 3686
        ## mouse_masked 0.03765 0.21670 0.10802 0.63690 1480    4246
        ## cat_masked   0.01862 0.09568 0.02514 0.14820 732 988
        ## dog_masked   0.00697 0.03815 0.01384 0.08160 274 544
        ##
        if os.path.isfile("refStats.txt"):
            refStatsFile = os.path.join(outDir, "refStats.txt")
            with open(refStatsFile) as RFH:
                ## Need to report 0 if nothing matched
                refStats['human'] = {"unambiguousReadsPerc":"0", "unambiguousMB":"0", "ambiguousReadsPerc":"0", "ambiguousMB":"0", "unambiguousReads":"0", "ambiguousReads":"0", "totalPerc":"0"}
                refStats['cat'] = {"unambiguousReadsPerc":"0", "unambiguousMB":"0", "ambiguousReadsPerc":"0", "ambiguousMB":"0", "unambiguousReads":"0", "ambiguousReads":"0", "totalPerc":"0"}
                refStats['dog'] = {"unambiguousReadsPerc":"0", "unambiguousMB":"0", "ambiguousReadsPerc":"0", "ambiguousMB":"0", "unambiguousReads":"0", "ambiguousReads":"0", "totalPerc":"0"}
                refStats['mouse'] = {"unambiguousReadsPerc":"0", "unambiguousMB":"0", "ambiguousReadsPerc":"0", "ambiguousMB":"0", "unambiguousReads":"0", "ambiguousReads":"0", "totalPerc":"0"}

                for l in RFH:
                    if l:
                        if l.startswith("#"):
                            continue

                        toks = l.rstrip().split()
                        assert len(toks) >= 7

                        ## the number and percent of reads that map unambiguously or ambiguously to human, cat, dog.
                        ## take the sum of the two numbers (ambiguous plus unambiguous) to use as the final percentage.
                        if l.startswith("human"):
                            refStats['human'] = {"unambiguousReadsPerc": toks[1], "unambiguousMB": toks[2], "ambiguousReadsPerc": toks[3], "ambiguousMB": toks[4], "unambiguousReads": toks[5], "ambiguousReads": toks[6], "totalPerc":float(toks[3])+float(toks[1])}
                        if l.startswith("cat"):
                            refStats['cat'] = {"unambiguousReadsPerc": toks[1], "unambiguousMB": toks[2], "ambiguousReadsPerc": toks[3], "ambiguousMB": toks[4], "unambiguousReads": toks[5], "ambiguousReads": toks[6], "totalPerc":float(toks[3])+float(toks[1])}
                        if l.startswith("dog"):
                            refStats['dog'] = {"unambiguousReadsPerc": toks[1], "unambiguousMB": toks[2], "ambiguousReadsPerc": toks[3], "ambiguousMB": toks[4], "unambiguousReads": toks[5], "ambiguousReads": toks[6], "totalPerc":float(toks[3])+float(toks[1])}
                        if l.startswith("mouse"):
                            refStats['mouse'] = {"unambiguousReadsPerc": toks[1], "unambiguousMB": toks[2], "ambiguousReadsPerc": toks[3], "ambiguousMB": toks[4], "unambiguousReads": toks[5], "ambiguousReads": toks[6], "totalPerc":float(toks[3])+float(toks[1])}

            log.debug("refStats.txt: %s", str(refStats))



        ###########################################################
        log_and_print("Write to stats file %s" % STATS_LIST_FILE_NAME)
        ###########################################################
        if os.path.isfile(STATS_LIST_FILE_NAME):
            os.remove(STATS_LIST_FILE_NAME)

        with open(BB_STATS_LIST_FILE_NAME) as bbfh:
            with open(STATS_LIST_FILE_NAME, 'a') as fh:
                for line in bbfh:
                    if not line.startswith("#") and line.strip():
                        fh.write(line)

        bbtoolsVersion = None

        stats = get_dict_obj(STATS_LIST_FILE_NAME)
        with open(STATS_LIST_FILE_NAME, 'a') as fh:
            for key in readCounts:
                if key not in stats:
                    write_stats(fh, key, readCounts[key], log)

            for key in refStats:
                for k in refStats[key]:
                    write_stats(fh, key+'_'+k, refStats[key][k], log)

            write_stats(fh, "cardinality", cardinality, log)

            ## Write refStats to filterStats.txt file
            for key in filterLogStat:
                for k in filterLogStat[key]:
                    write_stats(fh, key + '_' + k, filterLogStat[key][k], log)

            bbversionCmd = os.path.join(BBDIR, 'bbversion.sh')
            cmd = "%s" % (bbversionCmd)
            stdOut, _, exitCode = run_sh_command(cmd, True, log)
            assert stdOut is not None
            bbtoolsVersion = stdOut.strip()

            ## 05112017 Now bbtools version = bbmap version
            # bbtoolsVersion = bbmapVersion if bbmapVersion else "37.xx"
            assert bbtoolsVersion is not None
            write_stats(fh, "filter_tool", "bbtools " + bbtoolsVersion, log)
            write_stats(fh, "filter", VERSION, log)


            ## Version recording
            if bbdukVersion is None: bbdukVersion = bbtoolsVersion
            if bbmapVersion is None: bbmapVersion = bbtoolsVersion
            write_stats(fh, "bbduk_version", bbdukVersion, log)
            write_stats(fh, "bbmap_version", bbmapVersion, log)

        checkpoint(POST_END, status)
        status = POST_END
    else:
        log_and_print('No need to do post processing.')

    return filteredFastq, status


##==============================================================================
## Helper functions

def clean_up(fList, log):
    log_and_print("\n\n%s - CLEAN UP <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<%s\n" % ( color['pink'], color['']))

    for f in fList:
        if os.path.isfile(f):
            log_and_print("Removing %s ... ", f)
            os.remove(f)

    log_and_print("CLEAN UP - completed")


def create_shell(shName, cmdArray):
    with open(shName, 'w') as fh:
        fh.write("#!/bin/bash\n")
        fh.write("set -e\n")
        fh.write("set -o pipefail\n")
        #fh.write("module unload blast+; module load blast+\n")
        for item in cmdArray:
            fh.write("%s\n" % item)
        os.chmod(shName, 0755)  #-rwxr-xr-x


## checkpoint logging
def checkpoint(status, fromStatus=PIPE_START):
    if status == PIPE_START or STEP_ORDER[status] > STEP_ORDER[fromStatus]:
        checkpoint_step(STATUS_LOG_FNAME, status)

def write_stats(fh, k, v, log):
    line = "%s=%s" % (k, v)
    fh.write("%s\n" % line)

def read_qc(odir, fastq, status):
    log_and_print("\n\n%s - RUN QC <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<%s\n" % (color['pink'], color['']))
    if not os.path.isfile(fastq):
        return False, status

    qcdir = os.path.join(odir, READQC_DIR)

    if STEP_ORDER[status] < STEP_ORDER[QC_END]:
        checkpoint(QC_START, status)

        qcexe = os.path.join(os.path.dirname(__file__), 'readqc.py')
        cmd = '%s -o %s -f %s --skip-blast' % (qcexe, qcdir, fastq)
        # print('DEBUG : %s' % cmd)
        stdOut, stdErr, exitCode = run_sh_command(cmd, True)
        if exitCode != 0:
            print('ERROR : %s' % stdErr)
            return False, qcdir, status

        checkpoint(QC_END, status)
        status = QC_END
    else:
        log_and_print("No need to do qc step.")

    return True, qcdir, status

def do_html_body(odir, rawFastq, filteredFastq):

    stats = get_dict_obj(os.path.join(odir, STATS_LIST_FILE_NAME))
    tok_map = {
            'inputReads' : {'token' : '[_RAW-READ-CNT_]', 'type': 'bigint'},
            'inputBases' : {'token' : '[_RAW-BASE-CNT_]', 'type': 'bigint'},
            'outputReads' : {'token' : '[_FILTERED-READ-CNT_]', 'type': 'bigint'},
            'outputBases' : {'token' : '[_FILTERED-BASE-CNT_]', 'type': 'bigint'},
            'readRmPct' : {'token' : '[_REMOVED-READ-PCT_]', 'type': 'raw'},
            'baseRmPct' : {'token' : '[_REMOVED-BASE-PCT_]', 'type': 'raw'},

            'lowqualitydiscards_numreads' : {'token' : '[_LOW-QUAL-REMOVED-READ-CNT_]', 'type': 'bigint', 'filter': 0},
            'lowqualitydiscards_percreads' : {'token' : '[_LOW-QUAL-REMOVED-PCT_]', 'type': 'raw', 'filter': 0},

            'contaminants_numreads' : {'token' : '[_ARTI-REMOVED-READ-CNT_]', 'type': 'bigint', 'filter': 0},
            'contaminants_percreads' : {'token' : '[_ARTI-REMOVED-PCT_]', 'type': 'raw', 'filter': 0},
            'ribosomalsequenceremoved_numreads' : {'token' : '[_RRNA-REMOVED-READ-CNT_]', 'type': 'bigint', 'filter': 0},
            'ribosomalsequenceremoved_percreads' : {'token' : '[_RRNA-REMOVED-READ-PCT_]', 'type': 'raw', 'filter': 0},
            'microbialremoved_numreads' : {'token' : '[_MICROBE-REMOVED-READ-CNT_]', 'type': 'bigint', 'filter': 0},
            'microbialremoved_percreads' : {'token' : '[_MICROBE-REMOVED-READ-PCT_]', 'type': 'raw', 'filter': 0},
            'human_unambiguousreads' : {'token' : '[_HUMAN-REMOVED-READ-CNT_]', 'type': 'bigint', 'filter': 0},
            'human_unambiguousreadsperc' : {'token' : '[_HUMAN-REMOVED-READ-PCT_]', 'type': 'raw', 'filter': 0},
            'dog_unambiguousreads' : {'token' : '[_DOG-REMOVED-READ-CNT_]', 'type': 'bigint', 'filter': 0},
            'dog_unambiguousreadsperc' : {'token' : '[_DOG-REMOVED-READ-PCT_]', 'type': 'raw', 'filter': 0},
            'cat_unambiguousreads' : {'token' : '[_CAT-REMOVED-READ-CNT_]', 'type': 'bigint', 'filter': 0},
            'cat_unambiguousreadsperc' : {'token' : '[_CAT-REMOVED-READ-PCT_]', 'type': 'raw', 'filter': 0},
            'mouse_unambiguousreads' : {'token' : '[_MOUSE-REMOVED-READ-CNT_]', 'type': 'bigint', 'filter': 0},
            'mouse_unambiguousreadsperc' : {'token' : '[_MOUSE-REMOVED-READ-PCT_]', 'type': 'raw', 'filter': 0},
    }

    temp = os.path.join(PYDIR, 'template/filter_body_template.html')
    html = ''
    with open(temp, 'r') as fh:
        html = fh.read()

        ## do the place-holder replacement !!
        html = html.replace('[_RAW-FILE-LOCATION_]', rawFastq)
        html = html.replace('[_FILTERED-FILE-LOCATION_]', filteredFastq)
        fsize = format(os.stat(rawFastq).st_size / (1024*1024), ',')
        html = html.replace('[_RAW-FILE-SIZE_]', fsize)
        fsize = format(os.stat(filteredFastq).st_size / (1024*1024), ',')
        html = html.replace('[_FILTERED-FILE-SIZE_]', fsize)

        for key in tok_map:
            dat = tok_map[key]
            html = html.replace(dat['token'], pipeline_val(key, dat, stats))

        # readqc on the filter file
        if qcdir:
            hbody = do_readqc_html_body(qcdir, odir)
        else:
            hbody = ''
        html = html.replace('[_FILTERED-READ-QC_]', hbody)
    return html

def do_html(odir, qcdir, rawFastq, filteredFastq, status):
    log_and_print("\n\n%s - Create HTML file <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<%s\n" % (color['pink'], color['']))
    fname = os.path.basename(rawFastq)

    stats = get_dict_obj(os.path.join(odir, STATS_LIST_FILE_NAME))

    temp = os.path.join(PYDIR, 'template/template.html')
    with open(temp, 'r') as fh:
        html = fh.read()
        html = html.replace('[_PAGE-TITLE_]', 'Filter Report')
        html = html.replace('[_REPORT-TITLE_]', 'BBTools Filtering Report')
        html = html.replace('[_INPUT-FILE-NAME_]', fname)
        html = html.replace('[_REPORT-DATE_]', '{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()))

        hbody = do_html_body(odir, rawFastq, filteredFastq)
        html = html.replace('[_REPORT-BODY_]', hbody)

        fbasename = 'filter.log'
        fname = os.path.join(outputPath, fbasename)
        html = html.replace('[_FILTER-LOG_]', html_tag('a', fbasename, {'href': fbasename}))
        fsize = '%.1f' % (float(os.stat(fname).st_size) / 2014.0)
        html = html.replace('[_FILTER-LOG-SIZE_]', fsize)

        # fbasename = 'filter.txt'
        # fname = os.path.join(outputPath, fbasename)
        # html = html.replace('[_FILTER-REPORT_]', html_tag('a', fbasename, {'href': fbasename}))
        # fsize = '%.1f' % (float(os.stat(fname).st_size) / 2014.0)
        # html = html.replace('[_FILTER-REPORT-SIZE_]', fsize)

        ## write the html to file
        idxfile = os.path.join(odir, 'index.html')
        with open(idxfile, 'w') as fh2:
            fh2.write(html)
        print('HTML index file written to %s' % idxfile)

        # copy the css file
        cssdir = os.path.join(PYDIR, 'css')
        todir = os.path.join(odir, 'css')
        if os.path.isdir(todir):
            shutil.rmtree(todir)
        shutil.copytree(cssdir, todir, False, None)

        # copy the image file
        imgdir = os.path.join(PYDIR, 'images')
        todir = os.path.join(odir, 'images')
        if os.path.isdir(todir):
            shutil.rmtree(todir)
        shutil.copytree(imgdir, todir, False, None)

    return status

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## main program
if __name__ == "__main__":
    ## Parse options
    usage = "* Filter Pipeline, version %s\n" % (VERSION)
    origCmd = ' '.join(sys.argv)

    ## command line options
    parser = argparse.ArgumentParser(description=usage, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    PROD_TYPE = [
                    'DNA',
                    'FUNGAL',
                    'METAGENOME',
                    'VIRAL-METAGENOME',
                    'SAG',
                    'ISO',
                    'SAG',
                    'CELL-ENRICHMENT',
                    'PLANT-2x150',
                    'PLANT-2x250',
                    'RNA',
                    'SMRNA',
                    'METATRANSCRIPTOME',
                    'LFPE',
                    'CLRS',
                    'CLIP-PE',
                    'NEXTERA',
                    'ITAG',
                    'MICROTRANS',
                    'BISULPHITE',
                    '3PRIMERNA',
                    'CHIPSEQ',
                    'RNAWOHUMAN'
                    ]

    parser.add_argument("-f", "--fastq", dest="fastq", help="Set input fastq file (full path to fastq)", required=True)
    parser.add_argument("-o", "--output-path", dest="outputPath", help="Set output path to write to", required=True)
    parser.add_argument("-p", "--prod-type", dest="prodType", help="Set product type: %s" % PROD_TYPE, required=True)
    parser.add_argument("-rdb", "--ref-databases", dest="rqcfilterdata", help="Path to RQCFilterData dir", required=True)

    parser.add_argument("-t", "--taxlist", dest="taxList", help="A list of taxid(s) to exclude in CSV format", required=False)
    parser.add_argument("-ap", "--ap", dest="apNum", help="Set AP (Analysis Project) ID. Ex) -ap 123 or -ap 123,456,789", required=False)
    parser.add_argument("-at", "--at", dest="atNum", help="Set AT (Analysis Task) ID. Ex) -at 123 or -at 123,456,789", required=False)
    parser.add_argument("-v", "--version", action="version", version=VERSION)

    ## switches
    parser.add_argument("-qc", "--qc", dest="doqc", action="store_true", help="also perform readqc analysis on the filtered output", default=False)
    parser.add_argument("-a", "--aggressive", dest="enableAggressive", action="store_true", help="Enable aggressive=t and microbebuild=3", default=False)
    parser.add_argument("-b", "--skip-blast", dest="skipBlastFlag", action="store_true", help="Skip Blast search", default=False)
    parser.add_argument("-c", "--skip-cleanup", dest="skipCleanUp", action="store_true", help="Skip file clean up after pipeline complete", default=False)
    parser.add_argument("-d", "--debug", dest="doDebug", action="store_true", help="Enable debug mode")
    parser.add_argument("-l", "--disable-clumpify", dest="disableClumpify", action="store_true", help="Disable clumpify", default=False)
    parser.add_argument("-m", "--skip-microbes-removal", dest="disableRmoveMicrobes", action="store_true", help="Skip microbes removal", default=False)
    parser.add_argument("-r", "--contam", dest="enableRmoveMicrobes", action="store_true", help="Enable removemicrobes=t", default=False)
    parser.add_argument("-s", "--skip-subsampleqc", dest="skipSubsampleQc", action="store_true", help="Skip the subsample and qc step", default=False)
    parser.add_argument("-pl", "--print-log", dest="print_log", default = False, action = "store_true", help = "print log to screen")

    ## produce html when processing is done
    parser.add_argument("-html", "--html", action="store_true", help="Create html file", dest="html", default=False, required=False)

    skipCleanUp = False
    skipSubsampleQc = False
    outputPath = None ## output path, defaults to current working directory
    fastq = None ## full path to input fastq
    # logLevel = "DEBUG"
    enableRmoveMicrobes = False
    disableRmoveMicrobes = False
    disableClumpify = False
    enableAggressive = False
    skipBlastFlag = False
    apNum = None
    atNum = None

    taxList = ""

    options = parser.parse_args()
    print_log = options.print_log

    if options.outputPath:
        outputPath = options.outputPath

    if not outputPath:
        outputPath = os.getcwd()

    ## create output_directory if it doesn't exist
    if not os.path.isdir(outputPath):
        os.makedirs(outputPath)

    outputPath = os.path.realpath(outputPath)
    outputPath = outputPath.replace("/chos", "") if outputPath.startswith("/chos") else outputPath

    ## initialize my logger
    logFile = os.path.join(outputPath, "rqc_filter_pipeline.log")

    print "Started filtering pipeline with %s, writing log to: %s" % (SCRIPT_NAME, logFile)

    log = get_logger("filter", logFile, logLevel, print_log, True)

    if options.doDebug:
        DEBUG = True

    if options.skipCleanUp:
        skipCleanUp = True

    if options.skipBlastFlag:
        skipBlastFlag = True

    if options.skipSubsampleQc:
        skipSubsampleQc = True

    if options.enableRmoveMicrobes:
        if options.disableRmoveMicrobes:
            log.error("Conflict in option parameters: cannot set skip-contam with skip-microbes-removal.")
            sys.exit(0)

        enableRmoveMicrobes = True

    if options.disableRmoveMicrobes:
        if options.enableRmoveMicrobes:
            log.error("Conflict in option parameters: cannot set skip-contam with skip-microbes-removal.")
            sys.exit(0)

        disableRmoveMicrobes = True

    if options.enableAggressive:
        enableAggressive = True

    if options.disableClumpify:
        disableClumpify = True

    if options.fastq:
        fastq = options.fastq

    if options.taxList:
        taxList = options.taxList

    if options.apNum:
        apNum = options.apNum
    if options.atNum:
        atNum = options.atNum


    log_and_print("%s" % '#' * 80)
    log_and_print("  Filtering pipeline (version %s)" % VERSION)
    log_and_print("%s" % '#' * 80)


    prodType = options.prodType.upper()
    skipSubsampleQcFlag = options.skipSubsampleQc ## run subsample and qc or not

    if not os.path.isdir(outputPath):
        log.error("Cannot work with directory: %s", outputPath)

    ## check for fastq file
    if fastq:
        if not os.path.isfile(fastq):
            log.error("Input fastq file, %s not found, abort!", fastq)
    else:
        log.error("No fastq defined, abort!")

    fastq = os.path.realpath(fastq)
    fastq = fastq.replace("/chos", "") if fastq.startswith("/chos") else fastq


    ##--------------------------------
    ## init log
    log_and_print("Started pipeline, writing log to: %s" % logFile)
    log_and_print("CMD: %s" % origCmd)
    log_and_print("Fastq file: %s" % fastq)
    log_and_print("Output path: %s" % outputPath)

    os.chdir(outputPath)

    cycle = 0
    cycleMax = 1
    bIsPaired = False

    status = get_status(STATUS_LOG_FNAME, log)
    log_and_print("Starting pipeline at [%s]" % status)

    if status == PIPE_START:
        checkpoint(PIPE_START)

    ## main loop: retry upto cycleMax times
    while cycle < cycleMax:

        cycle += 1

        log_and_print("ATTEMPT [%d]" % cycle)

        filesToRemove = []  # list of intermediate files for clean up
        lastFastq = fastq   # lastFastq : fastq produced by each step, init to input
        rRnaFilterFile = None # MTF generates this 2nd filtered output file
        subsampledFastq = None
        bIsPaired = None
        filteredReadNum = -1
        filteredReadNumRrna = -1
        FragFile = SingletonFile = UnknownFile = None

        if cycle > 1:
            status = get_status(STATUS_LOG_FNAME, log)

        if status != PIPE_COMPLETE:
            ##
            ## Run rqcfilter.sh
            ##
            lastFastq, status = run_rqcfilter(lastFastq, outputPath, prodType, status, enableRmoveMicrobes, enableAggressive, disableRmoveMicrobes, disableClumpify, taxList, options.rqcfilterdata, log) ## Only MTF type generates rRnaFilterFile
            if filteredReadNum == 0:
                break

            ##--------------------------------
            ## Run post processing
            if lastFastq is not None and lastFastq != -1:
                lastFastq, status = post_process(fastq, outputPath, lastFastq, status, log)
            else:
                print "Failed @ rqcfilter"


            ##--------------------------------
            ## Clean up
            if lastFastq is not None and lastFastq != -1:

                ## run readQC on the filtered fastq
                if options.doqc:
                    rtn, qcdir, status = read_qc(outputPath, lastFastq, status)
                else:
                    qcdir = None
                    rtn = True

                ## create html file on the readQC results
                if rtn:
                    status = do_html(outputPath, qcdir, fastq, lastFastq, status)

                checkpoint(PIPE_COMPLETE)

                if not skipCleanUp:
                    clean_up(filesToRemove, log)
                else:
                    log_and_print("SKIP CLEANUP")

                log_and_print("Pipeline Completed")
                cycle = cycleMax + 1

            else:
                print "Failed @ postprocess"

        else:
            cycle = cycleMax + 1
            log_and_print("Pipeline already completed")


    print "Done."
    sys.exit(0)


## EOF

#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
stand alone read qc script based on jgi-rqc-pipeline/readqc/readqc.py v8.3.6

Command
    $ readqc.py -f FASTQ.GZ -o OUT_DIR --skip-blast -html

Outputs
  - normal QC outputs + index.html

Created: March 15, 2018

Shijie Yao (syao@lbl.gov)

Revision:

"""

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## libraries to use
import os
import sys
import argparse
import datetime
import shutil

SRC_ROOT = os.path.abspath(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(SRC_ROOT + "/lib")   # common

from readqc_constants import RQCReadQcConfig, RQCReadQc, ReadqcStats
from common import get_logger, get_status, append_rqc_stats, append_rqc_file, set_colors, get_subsample_rate,run_command
from readqc_utils import *
#from readqc_utils import checkpoint_step_wrapper, fast_subsample_fastq_sequences, write_unique_20_mers, illumina_read_gc
from rqc_fastq import get_working_read_length, read_count
from readqc_report import *
from os_utility import make_dir
from html_utility import html_tag, html_th, html_tr, html_link
from rqc_utility import get_dict_obj, pipeline_val

VERSION = "1.0.0"
LOG_LEVEL = "DEBUG"
SCRIPT_NAME = __file__

PYDIR = os.path.abspath(os.path.dirname(__file__))
BBDIR = os.path.join(PYDIR, os.path.pardir)

color = {}
color = set_colors(color, True)

"""
STEP1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"""
def do_fast_subsample_fastq_sequences(fastq, skipSubsampling, log):
    log.info("\n\n%sSTEP1 - Subsampling reads <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<%s\n", color['pink'], color[''])

    statsFile = RQCReadQcConfig.CFG["stats_file"]

    status = "1_illumina_readqc_subsampling in progress"
    checkpoint_step_wrapper(status)
    log.info("1_illumina_readqc_subsampling in progress.")

    inputReadNum = 0
    sampleRate = 0.0

    if not skipSubsampling:
        # sampleRate = RQCReadQc.ILLUMINA_SAMPLE_PCTENTAGE    ## 0.01
        inputReadNum = read_count(fastq)     ## read count from the original fastq
        assert inputReadNum > 0, "ERROR: invalid input fastq"
        sampleRate = get_subsample_rate(inputReadNum)
        log.info("Subsampling rate = %s", sampleRate)
    else:
        log.info("1_illumina_readqc_subsampling: skip subsampling. Use all the reads.")
        sampleRate = 1.0

    retCode = None
    totalReadNum = 0
    firstSubsampledFastqFileName = ""

    sequnitFileName = os.path.basename(fastq)
    sequnitFileName = sequnitFileName.replace(".fastq", "").replace(".gz", "")

    firstSubsampledFastqFileName = sequnitFileName + ".s" + str(sampleRate) + ".fastq"

    retCode, firstSubsampledFastqFileName, totalBaseCount, totalReadNum, subsampledReadNum, bIsPaired, readLength = fast_subsample_fastq_sequences(fastq, firstSubsampledFastqFileName, sampleRate, True, log)

    append_rqc_stats(statsFile, "SUBSAMPLE_RATE", sampleRate, log)

    if retCode in (RQCExitCodes.JGI_FAILURE, -2):
        log.info("1_illumina_readqc_subsampling failed.")
        status = "1_illumina_readqc_subsampling failed"
        checkpoint_step_wrapper(status)

    else:
        append_rqc_stats(statsFile, ReadqcStats.ILLUMINA_READ_BASE_COUNT, totalBaseCount, log)
        append_rqc_stats(statsFile, ReadqcStats.ILLUMINA_READ_COUNT, totalReadNum, log)

        log.info("1_illumina_readqc_subsampling complete.")
        status = "1_illumina_readqc_subsampling complete"
        checkpoint_step_wrapper(status)


    return status, firstSubsampledFastqFileName, totalReadNum, subsampledReadNum, bIsPaired, readLength


"""
STEP2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"""
def do_write_unique_20_mers(fastq, totalReadCount, log):
    log.info("\n\n%sSTEP2 - Sampling unique 25 mers <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<%s\n", color['pink'], color[''])

    if totalReadCount >= RQCReadQc.ILLUMINA_MER_SAMPLE_REPORT_FRQ * 2: ## 25000
        log.debug("read count total in step2 = %s", totalReadCount)
        filesFile = RQCReadQcConfig.CFG["files_file"]
        statsFile = RQCReadQcConfig.CFG["stats_file"]

        status = "2_unique_mers_sampling in progress"
        checkpoint_step_wrapper(status)
        log.info(status)

        retCode, newDataFile, newPngPlotFile, newHtmlPlotFile = write_unique_20_mers(fastq, log)

        if retCode != RQCExitCodes.JGI_SUCCESS:
            status = "2_unique_mers_sampling failed"
            log.error(status)
            checkpoint_step_wrapper(status)

        else:
            ## if no output files, skip this step.
            if newDataFile is not None:
                statsDict = {}

                ## in readqc_report.py
                ## 2014.07.23 read_level_mer_sampling is updated to process new output file format from bbcountunique
                log.info("2_unique_mers_sampling: post-processing the bbcountunique output file.")
                read_level_mer_sampling(statsDict, newDataFile, log)

                for k, v in statsDict.items():
                    append_rqc_stats(statsFile, k, str(v), log)

                ## outputs from bbcountunique
                append_rqc_file(filesFile, ReadqcStats.ILLUMINA_READ_20MER_UNIQUENESS_TEXT, newDataFile, log)
                append_rqc_file(filesFile, ReadqcStats.ILLUMINA_READ_20MER_UNIQUENESS_PLOT, newPngPlotFile, log)
                append_rqc_file(filesFile, ReadqcStats.ILLUMINA_READ_20MER_UNIQUENESS_D3_HTML_PLOT, newHtmlPlotFile, log)

            log.info("2_unique_mers_sampling complete.")
            status = "2_unique_mers_sampling complete"
            checkpoint_step_wrapper(status)

    else:
        ## if num reads < RQCReadQc.ILLUMINA_MER_SAMPLE_REPORT_FRQ = 25000
        ## just proceed to the next step
        log.warning("2_unique_mers_sampling can't run it because the number of reads < %s.", RQCReadQc.ILLUMINA_MER_SAMPLE_REPORT_FRQ * 2)
        status = "2_unique_mers_sampling complete"
        checkpoint_step_wrapper(status)


    return status


"""
STEP3 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"""
def do_illumina_read_gc(fastq, log):
    log.info("\n\n%sSTEP3 - Making read GC histograms <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<%s\n", color['pink'], color[''])

    filesFile = RQCReadQcConfig.CFG["files_file"]
    statsFile = RQCReadQcConfig.CFG["stats_file"]

    status = "3_illumina_read_gc in progress"
    checkpoint_step_wrapper(status)
    log.info("3_illumina_read_gc in progress.")

    reformat_gchist_file, png_file, htmlFile, mean_val, stdev_val, med_val, mode_val = illumina_read_gc(fastq, log)

    if not reformat_gchist_file:
        log.error("3_illumina_read_gc failed.")
        status = "3_illumina_read_gc failed"
        checkpoint_step_wrapper(status)

    else:
        append_rqc_stats(statsFile, ReadqcStats.ILLUMINA_READ_GC_MEAN, mean_val, log)
        append_rqc_stats(statsFile, ReadqcStats.ILLUMINA_READ_GC_STD, stdev_val, log)
        append_rqc_stats(statsFile, ReadqcStats.ILLUMINA_READ_GC_MED, med_val, log)
        append_rqc_stats(statsFile, ReadqcStats.ILLUMINA_READ_GC_MODE, mode_val, log)

        append_rqc_file(filesFile, ReadqcStats.ILLUMINA_READ_GC_TEXT, reformat_gchist_file, log)
        append_rqc_file(filesFile, ReadqcStats.ILLUMINA_READ_GC_PLOT, png_file, log)
        append_rqc_file(filesFile, ReadqcStats.ILLUMINA_READ_GC_D3_HTML_PLOT, htmlFile, log)

        log.info("3_illumina_read_gc complete.")
        status = "3_illumina_read_gc complete"
        checkpoint_step_wrapper(status)


    return status


"""
STEP4 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"""
def do_read_quality_stats(fastq, log):
    log.info("\n\n%sSTEP4 - Analyzing read quality <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<%s\n", color['pink'], color[''])

    filesFile = RQCReadQcConfig.CFG["files_file"]
    statsFile = RQCReadQcConfig.CFG["stats_file"]

    status = "4_illumina_read_quality_stats in progress"
    checkpoint_step_wrapper(status)
    log.info("4_illumina_read_quality_stats in progress.")

    readLength = 0
    readLenR1 = 0
    readLenR2 = 0
    isPairedEnd = None

    if not os.path.isfile(fastq):
        log.error("4_illumina_read_quality_stats failed. Cannot find the input fastq file")
        status = "4_illumina_read_quality_stats failed"
        checkpoint_step_wrapper(status)
        return status

    ## First figure out if it's pair-ended or not
    ## NOTE: ssize=10000 is recommended!
    readLength, readLenR1, readLenR2, isPairedEnd = get_working_read_length(fastq, log)
    log.info("Read length = %s, and is_pair_ended = %s", readLength, isPairedEnd)

    if readLength == 0:
        log.error("Failed to run get_working_read_length.")
        status = "4_illumina_read_quality_stats failed"
        checkpoint_step_wrapper(status)

    else:
        ## Pair-ended
        r1_r2_baseposqual_png = None   ## Average Base Position Quality Plot (*.qrpt.png)
        r1_r2_baseposqual_html = None  ## Average Base Position Quality D3 Plot (*.qrpt.html)
        r1_r2_baseposqual_txt = None   ## Read 1/2 Average Base Position Quality Text (*.qhist.txt)

        r1_baseposqual_box_png = None  ## Read 1 Average Base Position Quality Boxplot (*.r1.png)
        r2_baseposqual_box_png = None  ## Read 2 Average Base Position Quality Boxplot (*.r2.png)
        r1_baseposqual_box_html = None ## Read 1 Average Base Position Quality D3 Boxplot (*.r1.html)
        r2_baseposqual_box_html = None ## Read 2 Average Base Position Quality D3 Boxplot (*.r2.html)
        r1_r2_baseposqual_box_txt = None  ## Average Base Position Quality text

        r1_cyclenbase_png = None       ## Read 1 Percent N by Read Position (*.r1.fastq.base.stats.Npercent.png) --> Read 1 Cycle N Base Percent plot
        r2_cyclenbase_png = None       ## Read 2 Percent N by Read Position (*.r2.fastq.base.stats.Npercent.png) --> Read 2 Cycle N Base Percent plot

        r1_cyclenbase_txt = None       ## Read 1 Percent N by Read Position Text (*.r1.fastq.base.stats) --> Read 1 Cycle N Base Percent text
        #r2_cyclenbase_txt = None       ## Read 2 Percent N by Read Position Text (*.r2.fastq.base.stats) --> Read 2 Cycle N Base Percent text
        r1_r2_cyclenbase_txt = None    ## Merged Percent N by Read Position Text (*.r2.fastq.base.stats) --> Merged Cycle N Base Percent text

        r1_cyclenbase_html = None
        r2_cyclenbase_html = None

        r1_nuclcompfreq_png = None     ## Read 1 Nucleotide Composition Frequency Plot (*.r1.stats.png)
        r2_nuclcompfreq_png = None     ## Read 2 Nucleotide Composition Frequency Plot (*.r2.stats.png)
        r1_nuclcompfreq_html = None
        r2_nuclcompfreq_html = None

        ## Single-ended
        #se_baseposqual_txt = None      ## Average Base Position Quality Text (*.qrpt)
        #se_nuclcompfreq_png = None     ## Cycle Nucleotide Composition (*.stats.png)

        if isPairedEnd:
            append_rqc_stats(statsFile, ReadqcStats.ILLUMINA_READ_LENGTH_1, readLenR1, log)
            append_rqc_stats(statsFile, ReadqcStats.ILLUMINA_READ_LENGTH_2, readLenR2, log)

        else:
            append_rqc_stats(statsFile, ReadqcStats.ILLUMINA_READ_LENGTH_1, readLength, log)

        ## Average Base Position Quality Plot/Text using qhist.txt
        r1_r2_baseposqual_txt, r1_r2_baseposqual_png, r1_r2_baseposqual_html = gen_average_base_position_quality_plot(fastq, isPairedEnd, log) ## .reformat.qhist.txt

        log.debug("Outputs: %s %s %s", r1_r2_baseposqual_png, r1_r2_baseposqual_html, r1_r2_baseposqual_txt)

        append_rqc_file(filesFile, ReadqcStats.ILLUMINA_READ_QUAL_POS_PLOT_MERGED, r1_r2_baseposqual_png, log)
        append_rqc_file(filesFile, ReadqcStats.ILLUMINA_READ_QUAL_POS_PLOT_MERGED_D3_HTML_PLOT, r1_r2_baseposqual_html, log)
        append_rqc_file(filesFile, ReadqcStats.ILLUMINA_READ_QUAL_POS_QRPT_1, r1_r2_baseposqual_txt, log) ## for backward compatibility
        append_rqc_file(filesFile, ReadqcStats.ILLUMINA_READ_QUAL_POS_QRPT_2, r1_r2_baseposqual_txt, log) ## for backward compatibility
        append_rqc_file(filesFile, ReadqcStats.ILLUMINA_READ_QUAL_POS_QRPT_MERGED, r1_r2_baseposqual_txt, log)

        ## Average Base Position Quality Plot/Text using bqhist.txt
        r1_r2_baseposqual_box_txt, r1_baseposqual_box_png, r2_baseposqual_box_png, r1_baseposqual_box_html, r2_baseposqual_box_html = gen_average_base_position_quality_boxplot(fastq, log)

        log.debug("Read qual outputs: %s %s %s %s %s", r1_r2_baseposqual_box_txt, r1_baseposqual_box_png, r1_baseposqual_box_html, r2_baseposqual_box_png, r2_baseposqual_box_html)

        append_rqc_file(filesFile, ReadqcStats.ILLUMINA_READ_QUAL_POS_QRPT_BOXPLOT_1, r1_baseposqual_box_png, log)
        append_rqc_file(filesFile, ReadqcStats.ILLUMINA_READ_QUAL_POS_QRPT_D3_HTML_BOXPLOT_1, r1_baseposqual_box_html, log)
        append_rqc_file(filesFile, ReadqcStats.ILLUMINA_READ_QUAL_POS_QRPT_BOXPLOT_2, r2_baseposqual_box_png, log)
        append_rqc_file(filesFile, ReadqcStats.ILLUMINA_READ_QUAL_POS_QRPT_D3_HTML_BOXPLOT_2, r2_baseposqual_box_html, log)
        append_rqc_file(filesFile, ReadqcStats.ILLUMINA_READ_QUAL_POS_QRPT_BOXPLOT_TEXT, r1_r2_baseposqual_box_txt, log)

        ## ----------------------------------------------------------------------------------------------------
        ## compute Q20 of the two reads
        q20Read1 = None
        q20Read2 = None

        ## using bqhist.txt
        if r1_r2_baseposqual_box_txt:
            q20Read1 = q20_score_new(r1_r2_baseposqual_box_txt, 1, log)
            if isPairedEnd:
                q20Read2 = q20_score_new(r1_r2_baseposqual_box_txt, 2, log)

        log.debug("q20 for read 1 = %s", q20Read1)
        log.debug("q20 for read 2 = %s", q20Read2)

        if q20Read1 is not None:
            append_rqc_stats(statsFile, ReadqcStats.ILLUMINA_READ_Q20_READ1, q20Read1, log)
        else:
            log.error("Failed to get q20 read 1 from %s", r1_r2_baseposqual_box_txt)
            status = "4_illumina_read_quality_stats failed"
            checkpoint_step_wrapper(status)
            return status

        if isPairedEnd:
            if q20Read2 is not None:
                append_rqc_stats(statsFile, ReadqcStats.ILLUMINA_READ_Q20_READ2, q20Read2, log)
            else:
                log.error("Failed to get q20 read 2 from %s", r1_r2_baseposqual_box_txt)
                status = "4_illumina_read_quality_stats failed"
                checkpoint_step_wrapper(status)
                return status

        r1_r2_cyclenbase_txt, r1_nuclcompfreq_png, r1_nuclcompfreq_html, r2_nuclcompfreq_png, r2_nuclcompfreq_html = gen_cycle_nucleotide_composition_plot(fastq, readLength, isPairedEnd, log)

        log.debug("gen_cycle_nucleotide_composition_plot() ==> %s %s %s %s %s", r1_cyclenbase_txt, r1_nuclcompfreq_png, r1_nuclcompfreq_html, r2_nuclcompfreq_png, r2_nuclcompfreq_html)

        append_rqc_file(filesFile, ReadqcStats.ILLUMINA_READ_BASE_COUNT_TEXT_1, r1_r2_cyclenbase_txt, log)
        append_rqc_file(filesFile, ReadqcStats.ILLUMINA_READ_BASE_COUNT_TEXT_2, r1_r2_cyclenbase_txt, log) # reformat.sh generates a single merged output file.

        append_rqc_file(filesFile, ReadqcStats.ILLUMINA_READ_BASE_COUNT_PLOT_1, r1_nuclcompfreq_png, log)
        append_rqc_file(filesFile, ReadqcStats.ILLUMINA_READ_BASE_COUNT_D3_HTML_PLOT_1, r1_nuclcompfreq_html, log)

        if r2_nuclcompfreq_png:
            append_rqc_file(filesFile, ReadqcStats.ILLUMINA_READ_BASE_COUNT_PLOT_2, r2_nuclcompfreq_png, log)
        if r2_nuclcompfreq_html:
            append_rqc_file(filesFile, ReadqcStats.ILLUMINA_READ_BASE_COUNT_D3_HTML_PLOT_2, r2_nuclcompfreq_html, log)

        ### ---------------------------------------------------------------------------------------------------
        ## using bhist.txt
        r1_cyclenbase_txt, r1_cyclenbase_png, r1_cyclenbase_html, r2_cyclenbase_png, r2_cyclenbase_html = gen_cycle_n_base_percent_plot(fastq, readLength, isPairedEnd, log)

        log.debug("Outputs: %s %s %s", r1_cyclenbase_txt, r1_cyclenbase_png, r1_cyclenbase_html)

        append_rqc_file(filesFile, ReadqcStats.ILLUMINA_READ_BASE_PERCENTAGE_TEXT_1, r1_cyclenbase_txt, log)
        append_rqc_file(filesFile, ReadqcStats.ILLUMINA_READ_BASE_PERCENTAGE_TEXT_2, r1_cyclenbase_txt, log)

        append_rqc_file(filesFile, ReadqcStats.ILLUMINA_READ_BASE_PERCENTAGE_PLOT_1, r1_cyclenbase_png, log)
        append_rqc_file(filesFile, ReadqcStats.ILLUMINA_READ_BASE_PERCENTAGE_D3_HTML_PLOT_1, r1_cyclenbase_html, log)

        if r2_cyclenbase_png:
            append_rqc_file(filesFile, ReadqcStats.ILLUMINA_READ_BASE_PERCENTAGE_PLOT_2, r2_cyclenbase_png, log)
        if r2_cyclenbase_html:
            append_rqc_file(filesFile, ReadqcStats.ILLUMINA_READ_BASE_PERCENTAGE_D3_HTML_PLOT_2, r2_cyclenbase_html, log)

        log.info("4_illumina_read_quality_stats complete.")
        status = "4_illumina_read_quality_stats complete"
        checkpoint_step_wrapper(status)


    return status



"""
STEP5 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"""
def do_write_base_quality_stats(fastq, log):
    log.info("\n\n%sSTEP5 - Calculating base quality statistics for reads <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<%s\n", color['pink'], color[''])

    filesFile = RQCReadQcConfig.CFG["files_file"]
    statsFile = RQCReadQcConfig.CFG["stats_file"]

    status = "5_illumina_read_quality_stats in progress"
    checkpoint_step_wrapper(status)
    log.info("5_illumina_read_quality_stats in progress.")

    reformatObqhistFile = write_avg_base_quality_stats(fastq, log) ## *.reformat.obqhist.txt

    if not reformatObqhistFile:
        log.error("5_illumina_read_quality_stats failed.")
        status = "5_illumina_read_quality_stats failed"
        checkpoint_step_wrapper(status)

    else:
        ## Generate qual scores and plots of read level QC
        statsDict = {}

        retCode = base_level_qual_stats(statsDict, reformatObqhistFile, log)

        if retCode != RQCExitCodes.JGI_SUCCESS:
            log.error("5_illumina_read_quality_stats failed.")
            status = "5_illumina_read_quality_stats failed"
            checkpoint_step_wrapper(status)

        else:
            for k, v in statsDict.items():
                append_rqc_stats(statsFile, k, str(v), log)

            append_rqc_file(filesFile, ReadqcStats.ILLUMINA_READ_BASE_QUALITY_STATS, reformatObqhistFile, log)

            log.info("5_illumina_read_quality_stats complete.")
            status = "5_illumina_read_quality_stats complete"
            checkpoint_step_wrapper(status)


    return status



"""
STEP6 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"""
def do_illumina_count_q_score(fastq, log):
    log.info("\n\n%sSTEP6 - Generating quality score histogram <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<%s\n", color['pink'], color[''])

    filesFile = RQCReadQcConfig.CFG["files_file"]
    statsFile = RQCReadQcConfig.CFG["stats_file"]

    status = "6_illumina_count_q_score in progress"
    checkpoint_step_wrapper(status)
    log.info("6_illumina_count_q_score in progress.")

    qhistTxtFile, qhistPngFile, qhistHtmlPlotFile = illumina_count_q_score(fastq, log) ## *.obqhist.txt

    if not qhistTxtFile:
        log.error("6_illumina_count_q_score failed.")
        status = "6_illumina_count_q_score failed"
        checkpoint_step_wrapper(status)

    else:
        ## save qscores in statsFile
        qscore = {}
        read_level_qual_stats(qscore, qhistTxtFile, log)

        for k, v in qscore.items():
            append_rqc_stats(statsFile, k, str(v), log)

        append_rqc_file(filesFile, ReadqcStats.ILLUMINA_READ_QHIST_TEXT, qhistTxtFile, log)
        append_rqc_file(filesFile, ReadqcStats.ILLUMINA_READ_QHIST_PLOT, qhistPngFile, log)
        append_rqc_file(filesFile, ReadqcStats.ILLUMINA_READ_QHIST_D3_HTML_PLOT, qhistHtmlPlotFile, log)

        log.info("6_illumina_count_q_score complete.")
        status = "6_illumina_count_q_score complete"
        checkpoint_step_wrapper(status)


    return status



"""
STEP7 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"""
## 20140903 removed
##def do_illumina_calculate_average_quality(fastq, log):



"""
STEP8 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"""
def do_illumina_find_common_motifs(fastq, log):
    log.info("\n\n%sSTEP8 - Locating N stutter motifs in sequence reads <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<%s\n", color['pink'], color[''])

    filesFile = RQCReadQcConfig.CFG["files_file"]
    statsFile = RQCReadQcConfig.CFG["stats_file"]

    status = "8_illumina_find_common_motifs in progress"
    checkpoint_step_wrapper(status)
    log.info(status)

    retCode, statDataFile = illumina_find_common_motifs(fastq, log)

    log.info("nstutter statDataFile name = %s", statDataFile)

    if retCode != RQCExitCodes.JGI_SUCCESS:
        log.error("8_illumina_find_common_motifs failed.")
        status = "8_illumina_find_common_motifs failed"
        checkpoint_step_wrapper(status)

    else:
        ## read_level_stutter
        ## ex)
        ##688   N----------------------------------------------------------------------------------------------------------------------------
        ##-------------------------
        ##346   NNNNNNNNNNNNN----------------------------------------------------------------------------------------------------------------
        ##-------------------------
        ##53924 ------------N----------------------------------------------------------------------------------------------------------------
        ##-------------------------
        ##sum pct patterNs past 0.1 == 15.9245930330268 ( 54958 / 345114 * 100 )

        with open(statDataFile, "r") as stutFH:
            lines = stutFH.readlines()

            ## if no motifs are detected the file is empty
            if not lines:
                log.warning("The *.nstutter.stat file is not available in function read_level_stutter(). The function still returns JGI_SUCCESS.")

            else:
                assert lines[-1].find("patterNs") != -1
                t = lines[-1].split()
                ## ["sum", "pct", "patterNs", "past", "0.1", "==", "15.9245930330268", "(", "54958", "/", "345114", "*", "100", ")"]
                percent = "%.2f" % float(t[6])

                append_rqc_stats(statsFile, ReadqcStats.ILLUMINA_READ_N_FREQUENCE, percent, log)
                append_rqc_stats(statsFile, ReadqcStats.ILLUMINA_READ_N_PATTERN, str("".join(lines[:-1])), log)

        ## NOTE: ???
        append_rqc_file(filesFile, "find_common_motifs.dataFile", statDataFile, log)

        log.info("8_illumina_find_common_motifs complete.")
        status = "8_illumina_find_common_motifs complete"
        checkpoint_step_wrapper(status)


    return status


"""
STEP11 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"""
def do_illumina_detect_read_contam(fastq, bpToCut, log):
    log.info("\n\n%sSTEP11 - Detect read contam <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<%s\n", color['pink'], color[''])

    filesFile = RQCReadQcConfig.CFG["files_file"]
    statsFile = RQCReadQcConfig.CFG["stats_file"]

    status = "11_illumina_detect_read_contam in progress"
    checkpoint_step_wrapper(status)
    log.info("11_illumina_detect_read_contam in progress.")


    #########
    ## seal
    #########
    retCode2, outFileDict2, ratioResultDict2, contamStatDict = illumina_detect_read_contam3(fastq, bpToCut, log) ## seal version

    if retCode2 != RQCExitCodes.JGI_SUCCESS:
        log.error("11_illumina_detect_read_contam seal version failed.")
        status = "11_illumina_detect_read_contam failed"
        checkpoint_step_wrapper(status)

    else:
        for k, v in outFileDict2.items():
            append_rqc_file(filesFile, k + " seal", str(v), log)
            append_rqc_file(filesFile, k, str(v), log)

        for k, v in ratioResultDict2.items():
            append_rqc_stats(statsFile, k + " seal", str(v), log)
            append_rqc_stats(statsFile, k, str(v), log)

        ## contamination stat
        for k, v in contamStatDict.items():
            append_rqc_stats(statsFile, k, str(v), log)

        log.info("11_illumina_detect_read_contam seal version complete.")
        status = "11_illumina_detect_read_contam complete"
        checkpoint_step_wrapper(status)


    return status



"""
STEP13 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"""
## Removed!!
##def do_illumina_read_megablast(firstSubsampledFastqFileName, skipSubsampling, subsampledReadNum, log, blastDbPath=None):



"""
New STEP13 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"""
def do_illumina_subsampling_read_blastn(firstSubsampledFastqFileName, skipSubsampling, subsampledReadNum, log):
    log.info("\n\n%sSTEP13 - Run subsampling for Blast search <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<%s\n", color['pink'], color[''])

    subsampeldFastqFile = None
    totalReadNum = 0
    readNumToReturn = 0

    status = "13_illumina_subsampling_read_megablast in progress"
    checkpoint_step_wrapper(status)
    log.info("13_illumina_subsampling_read_megablast in progress.")

    if subsampledReadNum == 0:
        cmd = " ".join(["grep", "-c", "'^+'", firstSubsampledFastqFileName])

        stdOut, _, exitCode = run_command(cmd, True, log)

        if exitCode != 0:
            log.error("Failed to run grep cmd")
            return RQCExitCodes.JGI_FAILURE, None, None, None

        else:
            readNum = int(stdOut)

    else:
        readNum = subsampledReadNum

    log.info("Subsampled read number = %d", readNum)

    sampl_per = RQCReadQc.ILLUMINA_SAMPLE_PCTENTAGE    ## 0.01
    max_count = RQCReadQc.ILLUMINA_SAMPLE_COUNT         ## 50000

    ## TODO
    ## Use "samplereadstarget" option in reformat.sh

    if skipSubsampling:
        log.debug("No subsampling for megablast. Use fastq file, %s (readnum = %s) as query for megablast.", firstSubsampledFastqFileName, readNum)
        subsampeldFastqFile = firstSubsampledFastqFileName
        readNumToReturn = readNum

    else:
        if readNum > max_count:
            log.debug("Run the 2nd subsampling for running megablast.")
            secondSubsamplingRate = float(max_count) / readNum
            log.debug("SecondSubsamplingRate=%s, max_count=%s, readNum=%s", secondSubsamplingRate, max_count, readNum)
            log.info("Second subsampling of Reads after Percent Subsampling reads = %s with new percent subsampling %f.", readNum, secondSubsamplingRate)

            secondSubsampledFastqFile = ""
            # dataFile = ""

            sequnitFileName = os.path.basename(firstSubsampledFastqFileName)
            sequnitFileName = sequnitFileName.replace(".fastq", "").replace(".gz", "")

            secondSubsampledFastqFile = sequnitFileName + ".s" + str(sampl_per) + ".s" + str(secondSubsamplingRate) + ".n" + str(max_count) + ".fastq"

            ## ex) fq_sub_sample.pl -f .../7601.1.77813.CTTGTA.s0.01.fastq -o .../7601.1.77813.CTTGTA.s0.01.stats -r 0.0588142575171 > .../7601.1.77813.CTTGTA.s0.01.s0.01.s0.0588142575171.n50000.fastq
            ## ex) reformat.sh
            retCode, secondSubsampledFastqFile, totalBaseCount, totalReadNum, subsampledReadNum, _, _ = fast_subsample_fastq_sequences(firstSubsampledFastqFileName, secondSubsampledFastqFile, secondSubsamplingRate, False, log)

            if retCode != RQCExitCodes.JGI_SUCCESS:
                log.error("Second subsampling failed.")
                return RQCExitCodes.JGI_FAILURE, None, None, None

            else:
                log.info("Second subsampling complete.")
                log.info("Second Subsampling Total Base Count = %s.", totalBaseCount)
                log.info("Second Subsampling Total Reads = %s.", totalReadNum)
                log.info("Second Subsampling Sampled Reads = %s.", subsampledReadNum)

                if subsampledReadNum == 0:
                    log.warning("Too small first subsampled fastq file. Skip the 2nd sampling.")
                    secondSubsampledFastqFile = firstSubsampledFastqFileName

            subsampeldFastqFile = secondSubsampledFastqFile
            readNumToReturn = subsampledReadNum

        else:
            log.debug("The readNum is smaller than max_count=%s. The 2nd sampling skipped.", str(max_count))
            subsampeldFastqFile = firstSubsampledFastqFileName
            totalReadNum = readNum
            readNumToReturn = readNum

    log.info("13_illumina_subsampling_read_megablast complete.")
    status = "13_illumina_subsampling_read_megablast complete"
    checkpoint_step_wrapper(status)


    return status, subsampeldFastqFile, totalReadNum, readNumToReturn



##"""
##STEP14 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##"""
##def do_illumina_read_blastn_refseq_microbial(subsampeldFastqFile, subsampledReadNum, log, blastDbPath=None):
## 12212015 sulsj REMOVED!


"""
New STEP14 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"""
def do_illumina_read_blastn_refseq_archaea(subsampeldFastqFile, subsampledReadNum, log):
    log.info("\n\n%sSTEP14 - Run read blastn against refseq archaea <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<%s\n", color['pink'], color[''])

    statsFile = RQCReadQcConfig.CFG["stats_file"]
    retCode = None

    log.info("14_illumina_read_blastn_refseq_archaea in progress.")
    status = "14_illumina_read_blastn_refseq_archaea in progress"
    checkpoint_step_wrapper(status)

    retCode = illumina_read_blastn_refseq_archaea(subsampeldFastqFile, log)

    if retCode == RQCExitCodes.JGI_FAILURE:
        log.error("14_illumina_read_blastn_refseq_archaea failed.")
        status = "14_illumina_read_blastn_refseq_archaea failed"

    elif retCode == -143: ## timeout
        log.warning("14_illumina_read_blastn_refseq_archaea timeout.")
        status = "14_illumina_read_blastn_refseq_archaea complete"

    else:
        ## read number used in blast search?
        append_rqc_stats(statsFile, ReadqcStats.ILLUMINA_READS_NUMBER, subsampledReadNum, log)

        ret2 = read_megablast_hits("refseq.archaea", log)

        if ret2 != RQCExitCodes.JGI_SUCCESS:
            log.error("Errors in read_megablast_hits() of refseq.microbial")
            log.error("14_illumina_read_blastn_refseq_archaea reporting failed.")
            status = "14_illumina_read_blastn_refseq_archaea failed"

        else:
            log.info("14_illumina_read_blastn_refseq_archaea complete.")
            status = "14_illumina_read_blastn_refseq_archaea complete"


    return status


"""
STEP15 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"""
def do_illumina_read_blastn_refseq_bacteria(subsampeldFastqFile, log):
    log.info("\n\n%sSTEP15 - Run read blastn against refseq bacteria <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<%s\n", color['pink'], color[''])

    #statsFile = RQCReadQcConfig.CFG["stats_file"]
    #retCode = None

    log.info("15_illumina_read_blastn_refseq_bacteria in progress.")
    status = "15_illumina_read_blastn_refseq_bacteria in progress"
    checkpoint_step_wrapper(status)

    retCode = illumina_read_blastn_refseq_bacteria(subsampeldFastqFile, log)

    if retCode == RQCExitCodes.JGI_FAILURE:
        log.error("15_illumina_read_blastn_refseq_bacteria failed.")
        status = "15_illumina_read_blastn_refseq_bacteria failed"

    elif retCode == -143: ## timeout
        log.warning("15_illumina_read_blastn_refseq_bacteria timeout.")
        status = "15_illumina_read_blastn_refseq_bacteria complete"

    else:
        ret2 = read_megablast_hits("refseq.bacteria", log)

        if ret2 != RQCExitCodes.JGI_SUCCESS:
            log.error("Errors in read_megablast_hits() of refseq.microbial")
            log.error("15_illumina_read_blastn_refseq_bacteria reporting failed.")
            status = "15_illumina_read_blastn_refseq_bacteria failed"

        else:
            log.info("15_illumina_read_blastn_refseq_bacteria complete.")
            status = "15_illumina_read_blastn_refseq_bacteria complete"


    return status



"""
STEP16 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"""
def do_illumina_read_blastn_nt(subsampeldFastqFile, log):
    log.info("\n\n%sSTEP16 - Run read blastn against nt <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<%s\n", color['pink'], color[''])

    retCode = None

    log.info("16_illumina_read_blastn_nt in progress.")
    status = "16_illumina_read_blastn_nt in progress"
    checkpoint_step_wrapper(status)

    retCode = illumina_read_blastn_nt(subsampeldFastqFile, log)

    if retCode == RQCExitCodes.JGI_FAILURE:
        log.error("16_illumina_read_blastn_nt failed.")
        status = "16_illumina_read_blastn_nt failed"

    elif retCode == -143: ## timeout
        log.warning("16_illumina_read_blastn_nt timeout.")
        status = "16_illumina_read_blastn_nt complete"

    else:
        ret2 = read_megablast_hits("nt", log)

        if ret2 != RQCExitCodes.JGI_SUCCESS:
            log.error("Errors in read_megablast_hits() of nt")
            log.error("16_illumina_read_blastn_nt reporting failed.")
            status = "16_illumina_read_blastn_nt failed"

        else:
            log.info("16_illumina_read_blastn_nt complete.")
            status = "16_illumina_read_blastn_nt complete"


    return status



"""
STEP17 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

do_illumina_multiplex_statistics

Demultiplexing analysis for pooled lib

"""
def do_illumina_multiplex_statistics(fastq, log, isMultiplexed=None):
    log.info("\n\n%sSTEP17 - Run Multiplex statistics analysis <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<%s\n", color['pink'], color[''])

    filesFile = RQCReadQcConfig.CFG["files_file"]

    log.info("17_multiplex_statistics in progress.")
    status = "17_multiplex_statistics in progress"
    checkpoint_step_wrapper(status)

    log.debug("fastq file: %s", fastq)

    retCode, demultiplexStatsFile, detectionPngPlotFile, detectionHtmlPlotFile = illumina_generate_index_sequence_detection_plot(fastq, log, isMultiplexed=isMultiplexed)

    if retCode != RQCExitCodes.JGI_SUCCESS:
        log.error("17_multiplex_statistics failed.")
        status = "17_multiplex_statistics failed"
        checkpoint_step_wrapper(status)

    else:
        log.info("17_multiplex_statistics complete.")
        status = "17_multiplex_statistics complete"
        checkpoint_step_wrapper(status)

        if detectionPngPlotFile is not None:
            append_rqc_file(filesFile, ReadqcStats.ILLUMINA_READ_DEMULTIPLEX_STATS_PLOT, detectionPngPlotFile, log)

        if detectionHtmlPlotFile is not None:
            append_rqc_file(filesFile, ReadqcStats.ILLUMINA_READ_DEMULTIPLEX_STATS_D3_HTML_PLOT, detectionHtmlPlotFile, log)

        if demultiplexStatsFile is not None:
            append_rqc_file(filesFile, ReadqcStats.ILLUMINA_READ_DEMULTIPLEX_STATS, demultiplexStatsFile, log)


    return status



"""
STEP18 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

do_end_of_read_illumina_adapter_check

"""
def do_end_of_read_illumina_adapter_check(firstSubsampledFastqFileName, log):
    log.info("\n\n%sSTEP18 - Run end_of_read_illumina_adapter_check analysis <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<%s\n", color['pink'], color[''])

    filesFile = RQCReadQcConfig.CFG["files_file"]
    plotFile = None
    dataFile = None

    log.info("18_end_of_read_illumina_adapter_check in progress.")
    status = "18_end_of_read_illumina_adapter_check in progress"
    checkpoint_step_wrapper(status)

    log.debug("sampled fastq file: %s", firstSubsampledFastqFileName)

    retCode, dataFile, plotFile, htmlFile = end_of_read_illumina_adapter_check(firstSubsampledFastqFileName, log)

    if retCode != RQCExitCodes.JGI_SUCCESS:
        log.error("18_end_of_read_illumina_adapter_check failed.")
        status = "18_end_of_read_illumina_adapter_check failed"
        checkpoint_step_wrapper(status)

    else:
        log.info("18_end_of_read_illumina_adapter_check complete.")
        status = "18_end_of_read_illumina_adapter_check complete"
        checkpoint_step_wrapper(status)

        if plotFile is not None:
            append_rqc_file(filesFile, ReadqcStats.ILLUMINA_READ_END_OF_READ_ADAPTER_CHECK_PLOT, plotFile, log)
            append_rqc_file(filesFile, ReadqcStats.ILLUMINA_READ_END_OF_READ_ADAPTER_CHECK_D3_HTML_PLOT, htmlFile, log)

        if dataFile is not None:
            append_rqc_file(filesFile, ReadqcStats.ILLUMINA_READ_END_OF_READ_ADAPTER_CHECK_DATA, dataFile, log)


    return status



"""
STEP19 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

do_insert_size_analysis

"""
def do_insert_size_analysis(fastq, log):
    log.info("\n\n%sSTEP19 - Run insert size analysis <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<%s\n", color['pink'], color[''])

    filesFile = RQCReadQcConfig.CFG["files_file"]
    statsFile = RQCReadQcConfig.CFG["stats_file"]

    plotFile = None
    dataFile = None

    log.info("19_insert_size_analysis in progress.")
    status = "19_insert_size_analysis in progress"
    checkpoint_step_wrapper(status)

    log.debug("fastq file used: %s", fastq)

    retCode, dataFile, plotFile, htmlFile, statsDict = insert_size_analysis(fastq, log) ## by bbmerge.sh

    if retCode != RQCExitCodes.JGI_SUCCESS:
        log.error("19_insert_size_analysis failed.")
        status = "19_insert_size_analysis failed"
        checkpoint_step_wrapper(status)

    else:
        if plotFile is not None:
            append_rqc_file(filesFile, ReadqcStats.ILLUMINA_READ_INSERT_SIZE_HISTO_PLOT, plotFile, log)

        if dataFile is not None:
            append_rqc_file(filesFile, ReadqcStats.ILLUMINA_READ_INSERT_SIZE_HISTO_DATA, dataFile, log)

        if htmlFile is not None:
            append_rqc_file(filesFile, ReadqcStats.ILLUMINA_READ_INSERT_SIZE_HISTO_D3_HTML_PLOT, htmlFile, log)

        if statsDict:
            try:
                ## --------------------------------------------------------------------------------------------------------
                append_rqc_stats(statsFile, ReadqcStats.ILLUMINA_READ_INSERT_SIZE_AVG_INSERT, statsDict["avg_insert"], log)
                append_rqc_stats(statsFile, ReadqcStats.ILLUMINA_READ_INSERT_SIZE_STD_INSERT, statsDict["std_insert"], log)
                append_rqc_stats(statsFile, ReadqcStats.ILLUMINA_READ_INSERT_SIZE_MODE_INSERT, statsDict["mode_insert"], log)
                ## --------------------------------------------------------------------------------------------------------

                append_rqc_stats(statsFile, ReadqcStats.ILLUMINA_READ_INSERT_SIZE_TOTAL_TIME, statsDict["total_time"], log)
                append_rqc_stats(statsFile, ReadqcStats.ILLUMINA_READ_INSERT_SIZE_NUM_READS, statsDict["num_reads"], log)
                append_rqc_stats(statsFile, ReadqcStats.ILLUMINA_READ_INSERT_SIZE_JOINED_NUM, statsDict["joined_num"], log)
                append_rqc_stats(statsFile, ReadqcStats.ILLUMINA_READ_INSERT_SIZE_JOINED_PERC, statsDict["joined_perc"], log)

                append_rqc_stats(statsFile, ReadqcStats.ILLUMINA_READ_INSERT_SIZE_AMBIGUOUS_NUM, statsDict["ambiguous_num"], log)
                append_rqc_stats(statsFile, ReadqcStats.ILLUMINA_READ_INSERT_SIZE_AMBIGUOUS_PERC, statsDict["ambiguous_perc"], log)
                append_rqc_stats(statsFile, ReadqcStats.ILLUMINA_READ_INSERT_SIZE_NO_SOLUTION_NUM, statsDict["no_solution_num"], log)
                append_rqc_stats(statsFile, ReadqcStats.ILLUMINA_READ_INSERT_SIZE_NO_SOLUTION_PERC, statsDict["no_solution_perc"], log)
                append_rqc_stats(statsFile, ReadqcStats.ILLUMINA_READ_INSERT_SIZE_TOO_SHORT_NUM, statsDict["too_short_num"], log)
                append_rqc_stats(statsFile, ReadqcStats.ILLUMINA_READ_INSERT_SIZE_TOO_SHORT_PERC, statsDict["too_short_perc"], log)

                append_rqc_stats(statsFile, ReadqcStats.ILLUMINA_READ_INSERT_SIZE_INSERT_RANGE_START, statsDict["insert_range_start"], log)
                append_rqc_stats(statsFile, ReadqcStats.ILLUMINA_READ_INSERT_SIZE_INSERT_RANGE_END, statsDict["insert_range_end"], log)

                append_rqc_stats(statsFile, ReadqcStats.ILLUMINA_READ_INSERT_SIZE_90TH_PERC, statsDict["perc_90th"], log)
                append_rqc_stats(statsFile, ReadqcStats.ILLUMINA_READ_INSERT_SIZE_50TH_PERC, statsDict["perc_50th"], log)
                append_rqc_stats(statsFile, ReadqcStats.ILLUMINA_READ_INSERT_SIZE_10TH_PERC, statsDict["perc_10th"], log)
                append_rqc_stats(statsFile, ReadqcStats.ILLUMINA_READ_INSERT_SIZE_75TH_PERC, statsDict["perc_75th"], log)
                append_rqc_stats(statsFile, ReadqcStats.ILLUMINA_READ_INSERT_SIZE_25TH_PERC, statsDict["perc_25th"], log)

            except KeyError:
                log.error("19_insert_size_analysis failed (KeyError).")
                status = "19_insert_size_analysis failed"
                checkpoint_step_wrapper(status)
                return status

        log.info("19_insert_size_analysis complete.")
        status = "19_insert_size_analysis complete"
        checkpoint_step_wrapper(status)


    return status




"""
STEP21 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

do_sketch_vs_nt_refseq_silva

"""
def do_sketch_vs_nt_refseq_silva(fasta, log):
    log.info("\n\n%sSTEP21 - Run sketch vs nt, refseq, silva <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<%s\n", color['pink'], color[''])

    status = "21_sketch_vs_nt_refseq_silva in progress"
    checkpoint_step_wrapper(status)

    sketchOutDir = "sketch"
    sketchOutPath = os.path.join(outputPath, sketchOutDir)
    make_dir(sketchOutPath)
    change_mod(sketchOutPath, "0755")

    seqUnitName = os.path.basename(fasta)
    seqUnitName = file_name_trim(seqUnitName)

    filesFile = RQCReadQcConfig.CFG["files_file"]

    comOptions = "ow=t colors=f printtaxa=t depth depth2 unique2 merge"

    ## NT ##########################
    sketchOutFile = os.path.join(sketchOutPath, seqUnitName + ".sketch_vs_nt.txt")
    sendSketchShCmd = os.path.join(BBDIR, 'sendsketch.sh')
    cmd = "%s in=%s out=%s %s nt" % (sendSketchShCmd, fasta, sketchOutFile, comOptions)
    stdOut, stdErr, exitCode = run_sh_command(cmd, True, log, True)
    if exitCode != 0:
        log.error("Failed to run : %s, stdout : %s, stderr: %s", cmd, stdOut, stdErr)
        status = "21_sketch_vs_nt_refseq_silva failed"
        checkpoint_step_wrapper(status)
        return status

    append_rqc_file(filesFile, "sketch_vs_nt_output", sketchOutFile, log)

    ## Refseq ##########################
    sketchOutFile = os.path.join(sketchOutPath, seqUnitName + ".sketch_vs_refseq.txt")
    cmd = "%s in=%s out=%s %s refseq" % (sendSketchShCmd, fasta, sketchOutFile, comOptions)
    stdOut, stdErr, exitCode = run_sh_command(cmd, True, log, True)
    if exitCode != 0:
        log.error("Failed to run : %s, stdout : %s, stderr: %s", cmd, stdOut, stdErr)
        status = "21_sketch_vs_refseq failed"
        checkpoint_step_wrapper(status)
        return status

    append_rqc_file(filesFile, "sketch_vs_refseq_output", sketchOutFile, log)

    ## Silva ##########################
    sketchOutFile = os.path.join(sketchOutPath, seqUnitName + ".sketch_vs_silva.txt")
    cmd = "%s in=%s out=%s %s silva" % (sendSketchShCmd, fasta, sketchOutFile, comOptions)
    stdOut, stdErr, exitCode = run_sh_command(cmd, True, log, True)
    if exitCode != 0:
        log.error("Failed to run : %s, stdout : %s, stderr: %s", cmd, stdOut, stdErr)
        status = "21_sketch_vs_silva failed"
        checkpoint_step_wrapper(status)
        return status

    append_rqc_file(filesFile, "sketch_vs_silva_output", sketchOutFile, log)

    log.info("21_sketch_vs_nt_refseq_silva complete.")
    status = "21_sketch_vs_nt_refseq_silva complete"
    checkpoint_step_wrapper(status)

    return status



"""
STEP22 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

do_illumina_read_level_report_postprocessing

"""
## Removed!!
#def do_illumina_read_level_report(fastq, firstSubsampledFastqFileName, log):


"""
STEP23 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

do_cleanup_readqc

"""
def do_cleanup_readqc(log):
    log.info("\n\n%sSTEP23 - Cleanup <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<%s\n", color['pink'], color[''])

    status = "23_cleanup_readqc in progress"
    checkpoint_step_wrapper(status)

    ## Get option
    skipCleanup = RQCReadQcConfig.CFG["skip_cleanup"]

    if not skipCleanup:
        retCode = cleanup_readqc(log)
    else:
        log.warning("File cleaning is skipped.")
        retCode = RQCExitCodes.JGI_SUCCESS

    if retCode != RQCExitCodes.JGI_SUCCESS:
        log.error("23_cleanup_readqc failed.")
        status = "23_cleanup_readqc failed"
        checkpoint_step_wrapper(status)

    else:
        log.info("23_cleanup_readqc complete.")
        status = "23_cleanup_readqc complete"
        checkpoint_step_wrapper(status)


    return status


def stetch_section_note(table_header):
    overview = [
        'BBTools sketch.sh uses a technique called MinHash to rapidly compare large sequences.  ',
        'The result is similar to BLAST, a list of hits from a query sequence to various reference sequences, ',
        'sorted by similarity but the mechanisms are very different. ',
        'For more information, see <a href="http://bbtools.jgi.doe.gov/" target="_blank">http://bbtools.jgi.doe.gov</a>'
    ]

    legend = {
        'WKID'      : 'Weighted Kmer IDentity, which is the kmer identity compensating for differences in size.  So, comparing human chr1 to the full human genome would yield 100% WKID but approximately 10% KID.',
        'KID'       : 'Kmer IDentity, equal to matches/length; this is the fraction of shared kmers.',
        'ANI'       : 'Average Nucleotide Identity, derived from WKID and kmer length.',
        'Complt'    : 'Genome completeness (percent of the reference represented in the query).  Derived from WKID and KID.',
        'Contam'    : 'Contamination (percent of the query that does not match this reference, but matches some other reference).',
        'Depth'     : 'Per-kmer depth in sketches, indicating how many times that kmer was seen in the sequence.',
        'Depth2'    : 'Repeat-compensated depth',
        'Matches'   : 'The number of shared kmers between query and ref.',
        'Unique'    : 'The number of shared kmers between query and ref, and no other ref.',
        'Unique2'   : '??',
        'Depth2'    : 'Repeat-compensated depth',
        'noHit'     : 'Number of kmers that did not hit any reference sequence.  Though constant for a query, it will be reported differently for different references based on the relative size of the reference and query (if the reference is bigger than the query, it will report all of them).',
        'TaxID'     : 'NCBI taxonomic id, when available.',
        'gSize'     : 'Estimate of genomic size (number of unique kmers in the genome).  This is based on the smallest hash value in the list.  This is affected by blacklists or whitelists, and by using an assembly versus raw reads.',
        'gSeqs'     : 'Number of sequences used in the sketch.',
        'taxName'   : 'NCBI\'s name for that taxID.  If there is no taxID, the sequence name will be used.',
    }

    html = ''
    for name in table_header:
        if name != 'taxonomy':
            html += html_tag('li', '%s: %s' % (html_tag('span', name, {'class': 'notice'}), html_tag('span', legend.get(name), {'class': 'small-italic'})))
    html = html_tag('span', ''.join(overview), {'class': 'notice'}) + '<br /><br />Column Legend<br />' + html_tag('ul', html)

    return html


def sketch_table(fname):
    # print('DEBUG : %s' % fname)

    html = ''

    if os.path.isfile(fname):
        title = None
        header = None
        data = []

        with open(fname, 'r') as fh:
            for line in fh:
                line = line.strip()
                if line == '':
                    continue
                if line.startswith('Query:'):
                    title = line
                elif line.startswith('WKID'):
                    header = line.split('\t')
                else:
                    data.append(line)
        if title and header and data:
            html += html_th(header)
            for line in data:
                row = line.split('\t')
                html += html_tr(row)

            # html = html_tag('p', title) + html_tag('table', html, attrs={'class': 'data'})
            html = html_tag('table', html, attrs={'class': 'data'})
    return html, header

def do_html_contam_art_first_n_pb_tr(stats, files, odir, filepath_prefix):
    temp = os.path.join(PYDIR, 'template/readqc_artifacts.html')
    html = ''
    with open(temp, 'r') as fh:
        html = fh.read()

        ## the optional artifact contam
        artifact_tr = ''
        artifact_type = '50'
        artifact_val = pipeline_val('illumina read percent contamination artifact 50bp seal', {'type': 'raw'}, stats, files, filepath_prefix)
        artifact_file = pipeline_val('artifact_50bp.seal.stats seal', {'type': 'file'}, stats, files, filepath_prefix)
        if not artifact_file:
            artifact_file = pipeline_val('artifact_20bp.seal.stats seal', {'type': 'file'}, stats, files, filepath_prefix)
            artifact_type = '20'
            artifact_val = pipeline_val('illumina read percent contamination artifact 20bp seal', {'type': 'raw'}, stats, files, filepath_prefix)
        if artifact_file:
            html = html.replace('[_CONTAM-ART-FIRST-BP_]', artifact_type)
            html = html.replace('[_CONTAM-ART-FIRST-BP-SEAL_]', artifact_file)
            html = html.replace('[_CONTAM-ART-FIRST-BP-SEAL-PCT_]', artifact_val)

    return html


def do_html_body(odir, filepath_prefix):
    print('do_html_body - %s' % filepath_prefix)
    temp = os.path.join(PYDIR, 'template/readqc_body_template.html')
    statsf = os.path.join(odir, 'readqc_stats.txt')
    if not os.path.isfile(statsf):
        print('ERROR : file not found: %s' % statsf)
    stats = get_dict_obj(statsf)

    filesf = os.path.join(odir, 'readqc_files.txt')
    if not os.path.isfile(statsf):
        print('ERROR : file not found: %s' % statsf)
    files = get_dict_obj(os.path.join(odir, 'readqc_files.txt'))

    ## key (key name in readqc_stats or readqc_files) : {token (space holder in html template), type (value format)}
    tok_map = {
            ## Average Base Quality section
            'overall bases Q score mean' : {'token' : '[_BASE-QUALITY-SCORE_]', 'type': 'float', 'filter': 1},
            'overall bases Q score std' : {'token' : '[_BASE-QUALITY-SCORE-STD_]', 'type': 'float', 'filter': 1},
            'Q30 bases Q score mean' : {'token' : '[_Q30-BASE-QUALITY-SCORE_]', 'type': 'float', 'filter': 1},
            'Q30 bases Q score std' : {'token' : '[_Q30-BASE-QUALITY-SCORE-STD_]', 'type': 'float', 'filter': 1},
            'base C30' : {'token' : '[_COUNT-OF-BAESE-Q30_]', 'type': 'bigint'},
            'base C25' : {'token' : '[_COUNT-OF-BAESE-Q25_]', 'type': 'bigint'},
            'base C20' : {'token' : '[_COUNT-OF-BAESE-Q20_]', 'type': 'bigint'},
            'base C15' : {'token' : '[_COUNT-OF-BAESE-Q15_]', 'type': 'bigint'},
            'base C10' : {'token' : '[_COUNT-OF-BAESE-Q10_]', 'type': 'bigint'},
            'base C5' : {'token' : '[_COUNT-OF-BAESE-Q5_]', 'type': 'bigint'},
            'base Q30' : {'token' : '[_PCT-OF-BAESE-Q30_]', 'type': 'raw'},
            'base Q25' : {'token' : '[_PCT-OF-BAESE-Q25_]', 'type': 'raw'},
            'base Q20' : {'token' : '[_PCT-OF-BAESE-Q20_]', 'type': 'raw'},
            'base Q15' : {'token' : '[_PCT-OF-BAESE-Q15_]', 'type': 'raw'},
            'base Q10' : {'token' : '[_PCT-OF-BAESE-Q10_]', 'type': 'raw'},
            'base Q5' : {'token' : '[_PCT-OF-BAESE-Q5_]', 'type': 'raw'},

            ## Average Read Quality section
            'read Q30' : {'token' : '[_PCT-OF-READS-Q30_]', 'type': 'raw'},
            'read Q25' : {'token' : '[_PCT-OF-READS-Q25_]', 'type': 'raw'},
            'read Q20' : {'token' : '[_PCT-OF-READS-Q20_]', 'type': 'raw'},
            'read Q15' : {'token' : '[_PCT-OF-READS-Q15_]', 'type': 'raw'},
            'read Q10' : {'token' : '[_PCT-OF-READS-Q10_]', 'type': 'raw'},
            'read Q5' : {'token' : '[_PCT-OF-READS-Q5_]', 'type': 'raw'},

            'SUBSAMPLE_RATE' : {'token' : '[_SUBSAMPLE-RATE_]', 'type': 'pct', 'filter': 1},
            'read base quality stats' : {'token' : '[_AVG-READ-QUAL-HISTO-DATA_]', 'type': 'file', 'filter': 'link', 'label':'data file'},
            'ILLUMINA_READ_QHIST_D3_HTML_PLOT' : {'token' : '[_AVG-READ-QUAL-HOSTO-D3_]', 'type': 'file', 'filter': 'link', 'label':'interactive plot'},
            'read qhist plot' : {'token' : '[_AVG-READ-QUALITY-HISTOGRAM_]', 'type': 'file'},

            ## Average Base Position Quality section
            'read q20 read1' : {'token' : '[_READ_Q20_READ1_]', 'type': 'raw'},
            'read q20 read2' : {'token' : '[_READ_Q20_READ2_]', 'type': 'raw'},
            'read qual pos qrpt 1' : {'token' : '[_AVG-BASE-POS-QUAL-HISTO-DATA_]', 'type': 'file', 'filter': 'link', 'label':'data file'},
            'ILLUMINA_READ_QUAL_POS_PLOT_MERGED_D3_HTML_PLOT' : {'token' : '[_AVG-BASE-POS-QUAL-HISTO-D3_]', 'type': 'file', 'filter': 'link', 'label':'interactive plot'},
            'read qual pos plot merged' : {'token' : '[_AVG-BASE-POSITION-QUALITY_]', 'type': 'file'},

            ## Insert Size
            'ILLUMINA_READ_INSERT_SIZE_JOINED_PERC' : {'token' : '[_PCT-READS-JOINED_]', 'type': 'float', 'filter': 1},
            'ILLUMINA_READ_INSERT_SIZE_AVG_INSERT' : {'token' : '[_PCT-READS-JOINED-AVG_]', 'type': 'raw'},
            'ILLUMINA_READ_INSERT_SIZE_STD_INSERT' : {'token' : '[_PCT-READS-JOINED-STDDEV_]', 'type': 'raw'},
            'ILLUMINA_READ_INSERT_SIZE_MODE_INSERT' : {'token' : '[_PCT-READS-JOINED-MODE_]', 'type': 'raw'},

            'ILLUMINA_READ_INSERT_SIZE_HISTO_DATA' : {'token' : '[_INSERT-SIZE-HISTO-DATA_]', 'type': 'file', 'filter': 'link', 'label':'data file'},
            'ILLUMINA_READ_INSERT_SIZE_HISTO_D3_HTML_PLOT' : {'token' : '[_INSERT-SIZE-HISTO-D3_]', 'type': 'file', 'filter': 'link', 'label':'interactive plot'},
            'ILLUMINA_READ_INSERT_SIZE_HISTO_PLOT' : {'token' : '[_INSERT-SIZE-HISTOGRAM_]', 'type': 'file'},

            ## Read GC
            'read GC mean' : {'token' : '[_READ-GC-AVG_]', 'type': 'float', 'filter': 1},
            'read GC std' : {'token' : '[_READ-GC-STDDEV_]', 'type': 'float', 'filter': 1},
            'read GC text hist' : {'token' : '[_READ-QC-HISTO-DATA_]', 'type': 'file', 'filter': 'link', 'label':'data file'},
            'ILLUMINA_READ_GC_D3_HTML_PLOT' : {'token' : '[_READ-QC-HISTO-D3_]', 'type': 'file', 'filter': 'link', 'label':'interactive plot'},
            'read GC plot' : {'token' : '[_READ-GC-HIST_]', 'type': 'file'},

            ## Cycle Nucleotide Composition
            'read base count text 1' : {'token' : '[_NUC-COMP-FREQ-R1-DATA_]', 'type': 'file', 'filter': 'link', 'label':'data file'},
            'read base count text 2' : {'token' : '[_NUC-COMP-FREQ-R2-DATA_]', 'type': 'file', 'filter': 'link', 'label':'data file'},
            'ILLUMINA_READ_BASE_COUNT_D3_HTML_PLOT_1' : {'token' : '[_NUC-COMP-FREQ-R1-D3_]', 'type': 'file', 'filter': 'link', 'label':'interactive plot'},
            'ILLUMINA_READ_BASE_COUNT_D3_HTML_PLOT_2' : {'token' : '[_NUC-COMP-FREQ-R2-D3_]', 'type': 'file', 'filter': 'link', 'label':'interactive plot'},
            'read base count plot 1' : {'token' : '[_CYCLE-NUCL-COMPOSITION-READ1_]', 'type': 'file'},
            'read base count plot 2' : {'token' : '[_CYCLE-NUCL-COMPOSITION-READ2_]', 'type': 'file'},

            ## Percentage of Common Contaminants
            'illumina read percent contamination artifact seal' : {'token' : '[_CONTAM-ART-SEA-PCT_]', 'type': 'floatstr'},
            'artifact.seal.stats seal' : {'token' : '[_CONTAM-ART-SEAL_]', 'type': 'file'},

            'illumina read percent contamination DNA spikein seal' : {'token' : '[_DNA-SPIKEIN-SEAL_PCT_]', 'type': 'floatstr'},
            'DNA_spikein.seal.stats seal' : {'token' : '[_DNA-SPIKEIN-SEAL_]', 'type': 'file'},

            'illumina read percent contamination RNA spikein seal' : {'token' : '[_RNA-SPIKEIN-SEAL_PCT_]', 'type': 'floatstr'},
            'RNA_spikein.seal.stats seal' : {'token' : '[_RNA-SPIKEIN-SEAL_]', 'type': 'file'},

            'illumina read percent contamination fosmid seal' : {'token' : '[_CONTAM-FOSMID-SEAL-PCT_]', 'type': 'floatstr'},
            'fosmid.seal.stats seal' : {'token' : '[_CONTAM-FOSMID-SEAL_]', 'type': 'file'},

            'illumina read percent contamination fosmid seal' : {'token' : '[_CONTAM-FOSMID-SEAL-PCT_]', 'type': 'floatstr'},
            'fosmid.seal.stats seal' : {'token' : '[_CONTAM-FOSMID-SEAL_]', 'type': 'file'},

            'illumina read percent contamination mitochondrion seal' : {'token' : '[_CONTAM-MITO-SEAL-PCT_]', 'type': 'floatstr'},
            'mitochondrion.seal.stats seal' : {'token' : '[_CONTAM-MITO-SEAL_]', 'type': 'file'},

            'illumina read percent contamination plastid seal' : {'token' : '[_CONTAM-CHLO-SEAL-PCT_]', 'type': 'floatstr'},
            'plastid.seal.stats seal' : {'token' : '[_CONTAM-CHLO-SEAL_]', 'type': 'file'},

            'illumina read percent contamination phix seal' : {'token' : '[_CONTAM-PHIX-SEAL-PCT_]', 'type': 'floatstr'},
            'phix.seal.stats seal' : {'token' : '[_CONTAM-PHIX-SEAL_]', 'type': 'file'},

            'illumina read percent contamination rrna seal' : {'token' : '[_CONTAM-RRNA-SEAL-PCT_]', 'type': 'floatstr'},
            'rrna.seal.stats seal' : {'token' : '[_CONTAM-RRNA-SEAL_]', 'type': 'file'},

            'illumina read percent contamination microbes seal' : {'token' : '[_CONTAM-NON-SYN-SEAL-PCT_]', 'type': 'floatstr'},
            'microbes.seal.stats seal' : {'token' : '[_CONTAM-NON-SYN-SEAL_]', 'type': 'file'},

            # 'illumina read percent contamination synthetic seal' : {'token' : '[_CONTAM-SYN-SEAL-PCT_]', 'type': 'floatstr'},
            'synthetic.seal.stats seal' : {'token' : '[_CONTAM-SYN-SEAL_]', 'type': 'file'},

            'illumina read percent contamination adapters seal' : {'token' : '[_CONTAM-ADAPTER-PCT_]', 'type': 'floatstr'},
            'adapters.seal.stats seal' : {'token' : '[_CONTAM-ADAPTER-SEAL_]', 'type': 'file'},

            ## Sketch vs NT
            'sketch_vs_nt_output' : {'token' : '[_SKETCH-VS-NT_]', 'type': 'file'},
            # 'sketch_vs_nt_output' : {'token' : '[_SKETCH-VS-NT-BASE_]', 'type': 'file', 'filter': 'base'},
        }

    html = ''
    with open(temp, 'r') as fh:
        html = fh.read()

        for key in tok_map:
            dat = tok_map[key]
            val = pipeline_val(key, dat, stats, files, filepath_prefix)
            # print('key=%s; %s;         ====type=%s' % (key, val, type(val)))
            html = html.replace(dat['token'], val)

        artifact_tr = do_html_contam_art_first_n_pb_tr(stats, files, odir, filepath_prefix)
        html = html.replace('[_CONTAN-ART-SEAL-FIRST-BP_]', artifact_tr)

        ## synthetic contam
        contam_syn_file = pipeline_val('synthetic.seal.stats seal', {'type': 'file', 'filter': 'full'}, stats, files, filepath_prefix)
        if contam_syn_file and os.path.isfile(contam_syn_file):
            contam_syn_pct_key = 'illumina read percent contamination synthetic seal'
            cmd = 'grep "#Matched" %s' % contam_syn_file
            stdOut, stdErr, exitCode = run_sh_command(cmd, True)
            if exitCode == 0:
                toks = stdOut.strip().split()
                if len(toks) == 3:
                    html = html.replace('[_CONTAM-SYN-SEAL-PCT_]', toks[2][:-1])

        ###--- Sketch
        sketch_html = ''
        ## add Sketch vs NT if file exists:
        sketch_nt_file = pipeline_val('sketch_vs_nt_output', {'type': 'file', 'filter': 'full'}, stats, files, filepath_prefix)

        if os.path.isfile(sketch_nt_file):
            basename = os.path.basename(sketch_nt_file)
            # href = pipeline_val('sketch_vs_nt_output', {'type': 'file'}, stats, files, filepath_prefix)
            sketch_html = html_tag('h4', 'Sketch vs. NT') #+ '<br />' + 'Raw file:%s' % html_link(href, basename)
            sketch_tab, sketch_table_header = sketch_table(sketch_nt_file)
            sketch_html += sketch_tab


        ## add Sketch vs refseq if file exists?
        if sketch_html:
            sketch_html = html_tag('h2', 'Read Sketch', attrs={'class': 'section-title'}) \
                            + stetch_section_note(sketch_table_header)     \
                            + html_tag('div', sketch_html, {'class': 'section'})

        html = html.replace('[_SKETCH-TABLE_]', sketch_html)

        # copy the image file
        imgdir = os.path.join(PYDIR, 'images')
        todir = os.path.join(odir, 'images')
        if os.path.isdir(todir):
            shutil.rmtree(todir)
        shutil.copytree(imgdir, todir, False, None)

    return html


def do_html(odir, infile):
    # odir = os.path.abspath(odir)  ## DO NOT convert!! The relative file path in html need this original odir string matching

    screen('Output dir - %s' % odir)
    fname = os.path.basename(infile)
    screen('Create HTML page in %s for %s ..' % (odir, fname))
    temp = os.path.join(PYDIR, 'template/template.html')

    with open(temp, 'r') as fh:
        html = fh.read()
        html = html.replace('[_PAGE-TITLE_]', 'Read QC Report')
        html = html.replace('[_REPORT-TITLE_]', 'BBTools Read QC Report')
        html = html.replace('[_INPUT-FILE-NAME_]', fname)
        html = html.replace('[_REPORT-DATE_]', '{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()))

        hbody =do_html_body(odir, odir)
        html = html.replace('[_REPORT-BODY_]', hbody)

        ## write the html to file
        idxfile = os.path.join(odir, 'index.html')
        with open(idxfile, 'w') as fh2:
            fh2.write(html)
        screen('HTML index file written to %s' % idxfile)

    # copy the css file
    cssdir = os.path.join(PYDIR, 'css')
    todir = os.path.join(odir, 'css')
    if os.path.isdir(todir):
        shutil.rmtree(todir)
    shutil.copytree(cssdir, todir, False, None)

def screen(txt):
    print('.. %s' % txt)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Main Program
if __name__ == "__main__":
    desc = "RQC ReadQC Pipeline"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-f", "--fastq", dest="fastq", help="Set input fastq file (full path to fastq)", required=True)
    parser.add_argument("-o", "--output-path", dest="outputPath", help="Set output path to write to")
    parser.add_argument("-x", "--cut", dest="cutLenOfFirstBases", help="Set Read cut length (bp) for read contamination detection")
    parser.add_argument("-v", "--version", action="version", version=VERSION)

    ## For running on Crays
    parser.add_argument("-l", "--lib-name", dest="libName", help="Set library name", required=False)
    parser.add_argument("-m", "--is-multiplexed", dest="isMultiplexed", help="Set multiplexed data", required=False)
    parser.add_argument("-r", "--is-rna", dest="isRna", help="Set RNA data", required=False)

    ## produce html when processing is done
    parser.add_argument("-html", "--html", action="store_true", help="Create html file", dest="html", default=False, required=False)

    ## Toggle options
    parser.add_argument("-b", "--skip-blast", action="store_true", help="Skip the blast run", dest="skipBlast", default=False, required=False)

    parser.add_argument("-bn", "--skip-blast-nt", action="store_true", help="Skip Blast run against nt", dest="skipBlastNt", default=False, required=False)
    parser.add_argument("-br", "--skip-blast-refseq", action="store_true", help="Skip Blast run against refseq.archaea and refseq.bateria", dest="skipBlastRefseq", default=False, required=False)
    parser.add_argument("-c", "--skip-cleanup", action="store_true", help="Skip temporary file cleanup", dest="skipCleanup", default=False, required=False)
    parser.add_argument("-p", "--pooled-analysis", action="store_true", help="Enable pooled analysis (demultiplexing)", dest="pooledAnalysis", default=False, required=False)
    parser.add_argument("-s", "--skip-subsample", action="store_true", help="Skip subsampling.", dest="skipSubsample", default=False, required=False)
    parser.add_argument("-z", "--skip-localization", action="store_true", help="Skip database localization", dest="skipDbLocalization", default=False, required=False)

    options = parser.parse_args()

    outputPath = None # output path, defaults to current working directory
    fastq = None # full path to input fastq.gz

    status = "start"
    nextStepToDo = 0

    if options.outputPath:
        outputPath = options.outputPath
    else:
        outputPath = os.getcwd()

    if options.fastq:
        fastq = options.fastq

    ## create output_directory if it doesn't exist
    if not os.path.isdir(outputPath):
        os.makedirs(outputPath)

    libName = None ## for illumina_sciclone_analysis()
    isRna = None ## for illumina_sciclone_analysis()
    isMultiplexed = None ## for illumina_generate_index_sequence_detection_plot()

    bSkipBlast = False
    bSkipBlastRefseq = False ## refseq.archaea and refseq.bacteria
    bSkipBlastNt = False
    bPooledAnalysis = False

    ## RQC-743 Need to specify the first cutting bp for read contam detection (eps. for smRNA )
    firstBptoCut = 50
    if options.cutLenOfFirstBases:
        firstBptoCut = options.cutLenOfFirstBases

    if options.libName:
        libName = options.libName


    ## switches
    if options.isRna:
        isRna = options.isRna
    if options.isMultiplexed:
        isMultiplexed = options.isMultiplexed
    if options.skipBlast:
        bSkipBlast = options.skipBlast
    if options.skipBlastRefseq:
        bSkipBlastRefseq = options.skipBlastRefseq
    if options.skipBlastNt:
        bSkipBlastNt = options.skipBlastNt
    if options.pooledAnalysis:
        bPooledAnalysis = options.pooledAnalysis


    skipSubsampling = options.skipSubsample

    ## Set readqc config
    RQCReadQcConfig.CFG["status_file"] = os.path.join(outputPath, "readqc_status.log")
    RQCReadQcConfig.CFG["files_file"] = os.path.join(outputPath, "readqc_files.tmp")
    RQCReadQcConfig.CFG["stats_file"] = os.path.join(outputPath, "readqc_stats.tmp")
    RQCReadQcConfig.CFG["output_path"] = outputPath
    RQCReadQcConfig.CFG["skip_cleanup"] = options.skipCleanup
    RQCReadQcConfig.CFG["skip_localization"] = options.skipDbLocalization
    RQCReadQcConfig.CFG["log_file"] = os.path.join(outputPath, "readqc.log")

    screen("Started readqc pipeline, writing log to: %s" % (RQCReadQcConfig.CFG["log_file"]))

    log = get_logger("readqc", RQCReadQcConfig.CFG["log_file"], LOG_LEVEL, False, True)
    log.info("=================================================================")
    log.info("   Read Qc Analysis (version %s)", VERSION)
    log.info("=================================================================")

    log.info("Starting %s with %s", SCRIPT_NAME, fastq)

    if os.path.isfile(RQCReadQcConfig.CFG["status_file"]):
        status = get_status(RQCReadQcConfig.CFG["status_file"], log)
    else:
        checkpoint_step_wrapper(status)


    if not os.path.isfile(fastq):
        log.error("%s not found, aborting!", fastq)
        exit(2)
    elif status != "complete":
    ## check for fastq file
    # if os.path.isfile(fastq):
        log.info("Found %s, starting processing.", fastq)
        log.info("Latest status = %s", status)
        if status != 'start':
            nextStepToDo = int(status.split("_")[0])
            if status.find("complete") != -1:
                nextStepToDo += 1
            log.info("Next step to do = %s", nextStepToDo)

    if status != 'complete':
        bDone = False
        cycle = 0
        totalReadNum = 0
        firstSubsampledFastqFileName = ""
        secondSubsampledFastqFile = ""
        totalReadCount = 0
        subsampledReadNum = 0
        bIsPaired = False
        readLength = 0
        firstSubsampledLogFile = os.path.join(outputPath, "subsample", "first_subsampled.txt")
        secondSubsampledLogFile = os.path.join(outputPath, "subsample", "second_subsampled.txt")

        while not bDone:

            cycle += 1

            if bPooledAnalysis:
                nextStepToDo = 17
                status = "16_illumina_read_blastn_nt complete"

            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ## 1. fast_subsample_fastq_sequences
            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if nextStepToDo == 1 or status == "start":
                status, firstSubsampledFastqFileName, totalReadCount, subsampledReadNum, bIsPaired, readLength = do_fast_subsample_fastq_sequences(fastq, skipSubsampling, log)

                if status.endswith("failed"):
                    log.error("Subsampling failed.")
                    sys.exit(-1)

                ## if not skipSubsampling and subsampledReadNum == 0: ## too small input file
                if not skipSubsampling and subsampledReadNum < RQCReadQc.ILLUMINA_MIN_NUM_READS: ## min=1000
                    log.info("Too small input fastq file. Skip subsampling: total number of reads = %s, sampled number of reads = %s", totalReadCount, subsampledReadNum)
                    skipSubsampling = True
                    status, firstSubsampledFastqFileName, totalReadCount, subsampledReadNum, bIsPaired, readLength = do_fast_subsample_fastq_sequences(fastq, skipSubsampling, log)

                    if status.endswith("failed"):
                        log.error("Subsampling failed.")
                        sys.exit(-1)

                    subsampledReadNum = totalReadCount
                    if subsampledReadNum >= RQCReadQc.ILLUMINA_MIN_NUM_READS:
                        status = "1_illumina_readqc_subsampling complete"

                ## Still too low number of reads -> record ILLUMINA_TOO_SMALL_NUM_READS and quit
                ##
                if subsampledReadNum == 0 or subsampledReadNum < RQCReadQc.ILLUMINA_MIN_NUM_READS: ## min=1000
                    log.info("Too small number of reads (< %s). Stop processing.", RQCReadQc.ILLUMINA_MIN_NUM_READS)
                    print("WARNING : Too small number of reads (< %s). Stop processing." % RQCReadQc.ILLUMINA_MIN_NUM_READS)
                    log.info("Completed %s: %s", SCRIPT_NAME, fastq)

                    ## Add ILLUMINA_TOO_SMALL_NUM_READS=1 to stats to be used in reporting
                    statsFile = RQCReadQcConfig.CFG["stats_file"]
                    append_rqc_stats(statsFile, ReadqcStats.ILLUMINA_TOO_SMALL_NUM_READS, "1", log)

                    ## move rqc-files.tmp to rqc-files.txt
                    newFilesFile = os.path.join(outputPath, "readqc_files.txt")
                    newStatsFile = os.path.join(outputPath, "readqc_stats.txt")

                    cmd = "mv %s %s " % (RQCReadQcConfig.CFG["files_file"], newFilesFile)
                    log.info("mv cmd: %s", cmd)
                    run_sh_command(cmd, True, log)

                    cmd = "mv %s %s " % (RQCReadQcConfig.CFG["stats_file"], newStatsFile)
                    log.info("mv cmd: %s", cmd)
                    run_sh_command(cmd, True, log)

                    exit(0)

            if status == "1_illumina_readqc_subsampling failed":
                bDone = True


            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ## 2. write_unique_20_mers
            ## NOTE: fastq = orig input fastq,
            ##       totalReadCount = total reads in the orig input fastq
            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if nextStepToDo == 2 or status == "1_illumina_readqc_subsampling complete":
                ## Cope with restarting
                ## save subsamples file in "first_subsampled.txt" so that the file
                ## name can be read when started
                if not os.path.isfile(firstSubsampledLogFile):
                    make_dir(os.path.join(outputPath, "subsample"))
                    with open(firstSubsampledLogFile, "w") as SAMP_FH:
                        SAMP_FH.write(firstSubsampledFastqFileName + " " + str(totalReadCount) + " " + str(subsampledReadNum) + "\n")

                if firstSubsampledFastqFileName == "" and options.skipSubsample:
                    if firstSubsampledFastqFileName == "" and os.path.isfile(firstSubsampledLogFile):
                        with open(firstSubsampledLogFile, "r") as SAMP_FH:
                            l = SAMP_FH.readline().strip()
                            t = l.split()
                            assert len(t) == 2 or len(t) == 3
                            firstSubsampledFastqFileName = t[0]
                            totalReadCount = t[1]
                            subsampledReadNum = t[2]
                            log.debug("firstSubsampledFastqFileName=%s, totalReadCount=%s, subsampledReadNum=%s", firstSubsampledFastqFileName, totalReadCount, subsampledReadNum)

                    else:
                        nextStepToDo = 1
                        continue

                ## TODO: move getting totalReadNum in the func ???
                ##
                log.debug("firstSubsampledFastqFileName=%s, totalReadCount=%s, subsampledReadNum=%s", firstSubsampledFastqFileName, totalReadCount, subsampledReadNum)

                if totalReadCount is None or totalReadCount == 0:
                    nextStepToDo = 1
                    continue
                else:
                    status = do_write_unique_20_mers(fastq, totalReadCount, log)

            if status == "2_unique_mers_sampling failed":
                bDone = True

            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ## 3. generate read GC histograms: Illumina_read_gc
            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if nextStepToDo == 3 or status == "2_unique_mers_sampling complete":
                ## Read the subsampled file name
                if firstSubsampledFastqFileName == "" and os.path.isfile(firstSubsampledLogFile):
                    with open(firstSubsampledLogFile, "r") as SAMP_FH:
                        firstSubsampledFastqFileName = SAMP_FH.readline().strip().split()[0]
                status = do_illumina_read_gc(firstSubsampledFastqFileName, log)

            if status == "3_illumina_read_gc failed":
                bDone = True

            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ## 4. read_quality_stats
            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if nextStepToDo == 4 or status == "3_illumina_read_gc complete":
                if firstSubsampledFastqFileName == "" and os.path.isfile(firstSubsampledLogFile):
                    with open(firstSubsampledLogFile, "r") as SAMP_FH:
                        firstSubsampledFastqFileName = SAMP_FH.readline().strip().split()[0]
                status = do_read_quality_stats(firstSubsampledFastqFileName, log)

            if status == "4_illumina_read_quality_stats failed":
                bDone = True

            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ## 5. write_base_quality_stats
            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if nextStepToDo == 5 or status == "4_illumina_read_quality_stats complete":
                if firstSubsampledFastqFileName == "" and os.path.isfile(firstSubsampledLogFile):
                    with open(firstSubsampledLogFile, "r") as SAMP_FH:
                        firstSubsampledFastqFileName = SAMP_FH.readline().strip().split()[0]
                status = do_write_base_quality_stats(firstSubsampledFastqFileName, log)

            if status == "5_illumina_read_quality_stats failed":
                bDone = True

            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ## 6. illumina_count_q_score
            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if nextStepToDo == 6 or status == "5_illumina_read_quality_stats complete":
                if firstSubsampledFastqFileName == "" and os.path.isfile(firstSubsampledLogFile):
                    with open(firstSubsampledLogFile, "r") as SAMP_FH:
                        firstSubsampledFastqFileName = SAMP_FH.readline().strip().split()[0]
                status = do_illumina_count_q_score(firstSubsampledFastqFileName, log)

            if status == "6_illumina_count_q_score failed":
                bDone = True

            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ## 7. illumina_calculate_average_quality
            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if nextStepToDo == 7 or status == "6_illumina_count_q_score complete":
                ## Let's skip this step. (20140902)
                log.info("\n\n%sSTEP7 - Skipping 21mer analysis <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<%s\n", color['pink'], color[''])

                status = "7_illumina_calculate_average_quality in progress"
                checkpoint_step_wrapper(status)

                status = "7_illumina_calculate_average_quality complete"
                checkpoint_step_wrapper(status)

            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ## 8. illumina_find_common_motifs
            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if nextStepToDo == 8 or status == "7_illumina_calculate_average_quality complete":
                if firstSubsampledFastqFileName == "" and os.path.isfile(firstSubsampledLogFile):
                    with open(firstSubsampledLogFile, "r") as SAMP_FH:
                        firstSubsampledFastqFileName = SAMP_FH.readline().strip().split()[0]
                status = do_illumina_find_common_motifs(firstSubsampledFastqFileName, log)

            if status == "8_illumina_find_common_motifs failed":
                bDone = True

            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ## 10. illumina_run_tagdust
            ##
            ## NOTE: This step will be skipped. No need to run.
            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # if nextStepToDo == 10 or status == "9_illumina_run_dedupe complete":
            if nextStepToDo == 9 or status == "8_illumina_find_common_motifs complete":
                ## 20131023 skip this step
                log.info("\n\n%sSTEP10 - Skipping tag dust <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<%s\n", color['pink'], color[''])
                status = "10_illumina_run_tagdust in progress"
                checkpoint_step_wrapper(status)

                status = "10_illumina_run_tagdust complete"
                checkpoint_step_wrapper(status)

            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ## 11. illumina_detect_read_contam
            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if nextStepToDo == 11 or status == "10_illumina_run_tagdust complete":
                if firstSubsampledFastqFileName == "" and os.path.isfile(firstSubsampledLogFile):
                    with open(firstSubsampledLogFile, "r") as SAMP_FH:
                        firstSubsampledFastqFileName = SAMP_FH.readline().strip().split()[0]
                status = do_illumina_detect_read_contam(firstSubsampledFastqFileName, firstBptoCut, log)

            if status == "11_illumina_detect_read_contam failed":
                bDone = True


            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ## 13. illumina_subsampling_read_megablast
            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if nextStepToDo == 12 or status == "11_illumina_detect_read_contam complete":
                if not bSkipBlast:
                    if firstSubsampledFastqFileName == "" and os.path.isfile(firstSubsampledLogFile):
                        with open(firstSubsampledLogFile, "r") as SAMP_FH:
                            l = SAMP_FH.readline().strip()
                            firstSubsampledFastqFileName = l.split()[0]
                            totalReadCount = int(l.split()[1])
                            subsampledReadNum = int(l.split()[2])

                    status, secondSubsampledFastqFile, second_read_cnt_total, second_read_cnt_sampled = do_illumina_subsampling_read_blastn(firstSubsampledFastqFileName, skipSubsampling, subsampledReadNum, log)
                    log.debug("status=%s, secondSubsampledFastqFile=%s, second_read_cnt_total=%s, second_read_cnt_sampled=%s", status, secondSubsampledFastqFile, second_read_cnt_total, second_read_cnt_sampled)

                else:
                    statsFile = RQCReadQcConfig.CFG["stats_file"]
                    append_rqc_stats(statsFile, ReadqcStats.ILLUMINA_SKIP_BLAST, "1", log)
                    log.info("\n\n%sSTEP13 - Run subsampling for Blast search <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<%s\n", color['pink'], color[''])
                    log.info("13_illumina_subsampling_read_megablast skipped.\n")
                    status = "13_illumina_subsampling_read_megablast complete"
                    checkpoint_step_wrapper(status)

            if status == "13_illumina_subsampling_read_megablast failed":
                bDone = True

            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # 14. illumina_read_blastn_refseq_archaea
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if nextStepToDo == 14 or status == "13_illumina_subsampling_read_megablast complete":
                if not bSkipBlast and not bSkipBlastRefseq:
                    if secondSubsampledFastqFile == "" and not os.path.isfile(secondSubsampledLogFile):
                        nextStepToDo = 13
                        continue

                    if secondSubsampledFastqFile != "" and not os.path.isfile(secondSubsampledLogFile):
                        make_dir(os.path.join(outputPath, "subsample"))
                        with open(secondSubsampledLogFile, "w") as SAMP_FH:
                            SAMP_FH.write(secondSubsampledFastqFile + " " + str(second_read_cnt_total) + " " + str(second_read_cnt_sampled) + "\n")
                    else:
                        with open(secondSubsampledLogFile, "r") as SAMP_FH:
                            l = SAMP_FH.readline().strip()
                            secondSubsampledFastqFile = l.split()[0]
                            second_read_cnt_total = int(l.split()[1])
                            second_read_cnt_sampled = int(l.split()[2])

                    status = do_illumina_read_blastn_refseq_archaea(secondSubsampledFastqFile, second_read_cnt_sampled, log)
                    checkpoint_step_wrapper(status)

                else:
                    log.info("\n\n%sSTEP14 - Run illumina_read_blastn_refseq_archaea analysis <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<%s\n", color['pink'], color[''])
                    log.info("14_illumina_read_blastn_refseq_archaea skipped.\n")
                    status = "14_illumina_read_blastn_refseq_archaea in progress"
                    checkpoint_step_wrapper(status)
                    status = "14_illumina_read_blastn_refseq_archaea complete"
                    checkpoint_step_wrapper(status)

            if status == "14_illumina_read_blastn_refseq_archaea failed":
                bDone = True

            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # 15. illumina_read_blastn_refseq_bacteria
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if nextStepToDo == 15 or status == "14_illumina_read_blastn_refseq_archaea complete":
                if not bSkipBlast and not bSkipBlastRefseq:
                    if secondSubsampledFastqFile == "" and not os.path.isfile(secondSubsampledLogFile):
                        nextStepToDo = 13
                        continue

                    if secondSubsampledFastqFile != "" and not os.path.isfile(secondSubsampledLogFile):
                        make_dir(os.path.join(outputPath, "subsample"))
                        with open(secondSubsampledLogFile, "w") as SAMP_FH:
                            SAMP_FH.write(secondSubsampledFastqFile + " " + str(second_read_cnt_total) + " " + str(second_read_cnt_sampled) + "\n")
                    else:
                        with open(secondSubsampledLogFile, "r") as SAMP_FH:
                            l = SAMP_FH.readline().strip()
                            secondSubsampledFastqFile = l.split()[0]
                            second_read_cnt_total = int(l.split()[1])
                            second_read_cnt_sampled = int(l.split()[2])

                    status = do_illumina_read_blastn_refseq_bacteria(secondSubsampledFastqFile, log)
                    checkpoint_step_wrapper(status)

                else:
                    log.info("\n\n%sSTEP15 - Run illumina_read_blastn_refseq_bacteria analysis <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<%s\n", color['pink'], color[''])
                    log.info("15_illumina_read_blastn_refseq_bacteria skipped.\n")
                    status = "15_illumina_read_blastn_refseq_bacteria in progress"
                    checkpoint_step_wrapper(status)
                    status = "15_illumina_read_blastn_refseq_bacteria complete"
                    checkpoint_step_wrapper(status)

            if status == "15_illumina_read_blastn_refseq_bacteria failed":
                bDone = True

            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ## 16. illumina_read_blastn_nt
            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if nextStepToDo == 16 or status == "15_illumina_read_blastn_refseq_bacteria complete":
                if not bSkipBlast and not bSkipBlastNt:
                    if secondSubsampledFastqFile == "" and not os.path.isfile(secondSubsampledLogFile):
                        nextStepToDo = 13
                        continue

                    if secondSubsampledFastqFile == "" and os.path.isfile(secondSubsampledLogFile):
                        with open(secondSubsampledLogFile, "r") as SAMP_FH:
                            l = SAMP_FH.readline().strip()
                            secondSubsampledFastqFile = l.split()[0]
                            second_read_cnt_total = int(l.split()[1])
                            second_read_cnt_sampled = int(l.split()[2])

                    status = do_illumina_read_blastn_nt(secondSubsampledFastqFile, log)
                    checkpoint_step_wrapper(status)

                else:
                    log.info("\n\n%sSTEP16 - Run illumina_read_blastn_nt analysis <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<%s\n", color['pink'], color[''])
                    log.info("16_illumina_read_blastn_nt skipped.\n")
                    status = "16_illumina_read_blastn_nt in progress"
                    checkpoint_step_wrapper(status)
                    status = "16_illumina_read_blastn_nt complete"
                    checkpoint_step_wrapper(status)

            if status == "16_illumina_read_blastn_nt failed":
                bDone = True

            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ## 17. multiplex_statistics
            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if nextStepToDo == 17 or status == "16_illumina_read_blastn_nt complete":
                status = do_illumina_multiplex_statistics(fastq, log, isMultiplexed=isMultiplexed)

            if status == "17_multiplex_statistics failed":
                bDone = True

            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ## 18. end_of_read_illumina_adapter_check
            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if nextStepToDo == 18 or status == "17_multiplex_statistics complete":
                if firstSubsampledFastqFileName == "" and os.path.isfile(firstSubsampledLogFile):
                    with open(firstSubsampledLogFile, "r") as SAMP_FH:
                        firstSubsampledFastqFileName = SAMP_FH.readline().strip().split()[0]
                status = do_end_of_read_illumina_adapter_check(firstSubsampledFastqFileName, log)

            if status == "18_end_of_read_illumina_adapter_check failed":
                bDone = True

            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ## 19. insert size analysis (bbmerge.sh)
            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if nextStepToDo == 19 or status == "18_end_of_read_illumina_adapter_check complete":
                status = do_insert_size_analysis(fastq, log)

            if status == "19_insert_size_analysis failed":
                bDone = True


            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ## 21. sketch vs nt, refseq, silva
            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # if nextStepToDo == 21 or status == "20_gc_divergence_analysis complete":
            if nextStepToDo == 20 or status == "19_insert_size_analysis complete":
                status = do_sketch_vs_nt_refseq_silva(fastq, log)

            if status == "21_sketch_vs_nt_refseq_silva failed":
                bDone = True

            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ## 22. postprocessing & reporting
            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if nextStepToDo == 22 or status == "21_sketch_vs_nt_refseq_silva complete":
                ## 20131023 skip this step
                log.info("\n\n%sSTEP22 - Run illumina_readqc_report_postprocess: mv rqc-*.tmp to rqc-*.txt <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<%s\n", color['pink'], color[''])
                status = "22_illumina_readqc_report_postprocess in progress"
                checkpoint_step_wrapper(status)

                ## move rqc-files.tmp to rqc-files.txt
                newFilesFile = os.path.join(outputPath, "readqc_files.txt")
                newStatsFile = os.path.join(outputPath, "readqc_stats.txt")

                cmd = "mv %s %s " % (RQCReadQcConfig.CFG["files_file"], newFilesFile)
                log.info("mv cmd: %s", cmd)
                run_sh_command(cmd, True, log)

                cmd = "mv %s %s " % (RQCReadQcConfig.CFG["stats_file"], newStatsFile)
                log.info("mv cmd: %s", cmd)
                run_sh_command(cmd, True, log)

                log.info("22_illumina_readqc_report_postprocess complete.")
                status = "22_illumina_readqc_report_postprocess complete"
                checkpoint_step_wrapper(status)

            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ## 23. Cleanup
            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if nextStepToDo == 23 or status == "22_illumina_readqc_report_postprocess complete":
                status = do_cleanup_readqc(log)

            if status == "23_cleanup_readqc failed":
                bDone = True

            if status == "23_cleanup_readqc complete":
                status = "complete" ## FINAL COMPLETE!
                bDone = True

            ## don't cycle more than 10 times ...
            if cycle > 10:
                bDone = True

        if status != "complete":
            log.info("Status %s", status)
        else:
            log.info("\n\nCompleted %s: %s", SCRIPT_NAME, fastq)
            checkpoint_step_wrapper("complete")
            log.info("Pipeline processing Done.")
            print("Pipeline processing Done.")

    else:
        log.info("Pipeline processing is already complete, skip.")
        print("Pipeline processing is already complete, skip.")

    if status == 'complete' and options.html:
            do_html(outputPath, fastq)


    exit(0)


## EOF
